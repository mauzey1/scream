from utils import run_cmd, run_cmd_no_fail, check_minimum_python_version, get_current_head,     \
    get_current_commit, get_current_branch, expect, is_repo_clean, cleanup_repo,  \
    get_common_ancestor, merge_git_ref, checkout_git_ref, print_last_commit

from machines_specs import get_mach_compilation_resources, get_mach_testing_resources, \
    get_mach_baseline_root_dir, setup_mach_env, is_cuda_machine, \
    get_mach_cxx_compiler, get_mach_f90_compiler, get_mach_c_compiler, is_machine_supported

check_minimum_python_version(3, 4)

import os, shutil, pathlib
import concurrent.futures as threading3
import itertools
import json

###############################################################################
class TestAllScream(object):
###############################################################################

    ###########################################################################
    def __init__(self, cxx_compiler, f90_compiler, c_compiler, submit=False, parallel=False, fast_fail=False,
                 baseline_ref=None, baseline_dir=None, machine=None, no_tests=False, keep_tree=False,
                 custom_cmake_opts=(), custom_env_vars=(), preserve_env=False, tests=(),
                 integration_test="JENKINS_HOME" in os.environ, local=False, root_dir=None, work_dir=None,
                 quick_rerun=False,quick_rerun_failed=False,dry_run=False,
                 make_parallel_level=0, ctest_parallel_level=0):
    ###########################################################################

        self._cxx_compiler            = cxx_compiler
        self._f90_compiler            = f90_compiler
        self._c_compiler              = c_compiler
        self._submit                  = submit
        self._parallel                = parallel
        self._fast_fail               = fast_fail
        self._baseline_ref            = baseline_ref
        self._machine                 = machine
        self._local                   = local
        self._perform_tests           = not no_tests
        self._keep_tree               = keep_tree
        self._baseline_dir            = baseline_dir
        self._custom_cmake_opts       = custom_cmake_opts
        self._custom_env_vars         = custom_env_vars
        self._preserve_env            = preserve_env
        self._tests                   = tests
        self._root_dir                = root_dir
        self._work_dir                = work_dir
        self._integration_test        = integration_test
        self._quick_rerun             = quick_rerun
        self._quick_rerun_failed      = quick_rerun_failed
        self._dry_run                 = dry_run
        self._must_generate_baselines = False

        if self._quick_rerun_failed:
            self._quick_rerun = True

        ############################################
        #  Sanity checks and helper structs setup  #
        ############################################

        # Probe machine if none was specified
        if self._machine is None:
            # We could potentially integrate more with CIME here to do actual
            # nodename probing.
            if "CIME_MACHINE" in os.environ and is_machine_supported(os.environ["CIME_MACHINE"]):
                self._machine = os.environ["CIME_MACHINE"]
            else:
                expect(self._local,
                       "test-all-scream requires either the machine arg (-m $machine) or the -l flag,"
                       "which makes it lookf for machine specs in '~/.cime/scream_mach_specs.py'.")
                self._machine = "local"
        else:
            expect (not self._local, "Specifying a machine while passing '-l,--local' is ambiguous.")

        # Unless the user claims to know what he/she is doing, we setup the env.
        if not self._preserve_env:
            # Setup the env on this machine
            setup_mach_env(self._machine)

        # Compute root dir
        if not self._root_dir:
            self._root_dir = pathlib.Path(__file__).resolve().parent.parent
        else:
            self._root_dir = pathlib.Path(self._root_dir).resolve()
            expect(self._root_dir.is_dir() and self._root_dir.parts()[-2:] == ('scream', 'components'),
                   "Bad root-dir '{}', should be: $scream_repo/components/scream".format(self._root_dir))

        if self._work_dir is not None:
            expect(pathlib.Path(self._work_dir).absolute().is_dir(),
                   "Error! Work directory '{}' does not exist.".format(self._work_dir))
        else:
            self._work_dir = self._root_dir.absolute().joinpath("ctest-build")

        expect (not self._baseline_dir or self._work_dir != self._baseline_dir,
                "Error! For your safety, do NOT use '{}' to store baselines. Move them to a different directory (even a subdirectory of that works).".format(self._work_dir))

        expect(not (self._baseline_ref and self._baseline_dir),
               "Makes no sense to specify a baseline generation commit if using pre-existing baselines ")

        self._tests_cmake_args = {"dbg" : [("CMAKE_BUILD_TYPE", "Debug"),
                                           ("EKAT_DEFAULT_BFB", "True")],
                                  "sp"  : [("CMAKE_BUILD_TYPE", "Debug"),
                                           ("SCREAM_DOUBLE_PRECISION", "False"),
                                           ("EKAT_DEFAULT_BFB", "True")],
                                  "fpe" : [("CMAKE_BUILD_TYPE", "Debug"),
                                           ("SCREAM_PACK_SIZE", "1"),
                                           ("SCREAM_SMALL_PACK_SIZE", "1"),
                                           ("EKAT_DEFAULT_BFB", "True")]}

        self._test_full_names = { "dbg" : "full_debug",
                                  "sp"  : "full_sp_debug",
                                  "fpe" : "debug_nopack_fpe"}

        if not self._tests:
            # always do dbg and sp tests, do not do fpe test on CUDA
            self._tests = ["dbg", "sp"]
            if not is_cuda_machine(self._machine):
                self._tests.append("fpe")
        else:
            for t in self._tests:
                expect(t in self._test_full_names,
                       "Requested test '{}' is not supported by test-all-scream, please choose from: {}".\
                           format(t, ", ".join(self._test_full_names.keys())))

        os.chdir(str(self._root_dir)) # needed, or else every git command will need repo=root_dir
        expect(get_current_commit(), "Root dir: {}, does not appear to be a git repo".format(self._root_dir))

        self._original_branch = get_current_branch()
        self._original_commit = get_current_commit()

        print_last_commit(git_ref=self._original_branch, dry_run=self._dry_run)

        ############################################
        #    Deduce compilers if needed/possible   #
        ############################################

        if self._cxx_compiler is None:
            self._cxx_compiler = get_mach_cxx_compiler(self._machine)
        if self._f90_compiler is None:
            self._f90_compiler = get_mach_f90_compiler(self._machine)
        if self._c_compiler is None:
            self._c_compiler = get_mach_c_compiler(self._machine)

        if not self._dry_run:
            self._f90_compiler = run_cmd_no_fail("which {}".format(self._f90_compiler))
            self._cxx_compiler = run_cmd_no_fail("which {}".format(self._cxx_compiler))
            self._c_compiler   = run_cmd_no_fail("which {}".format(self._c_compiler))

        ###################################
        #      Compute baseline info      #
        ###################################

        default_baselines_root_dir = pathlib.Path(self._work_dir,"baselines")
        if self._baseline_dir is None:
            if self._baseline_ref is None:
                # Compute baseline ref
                if self._keep_tree:
                    self._baseline_ref = "HEAD"
                elif self._integration_test:
                    self._baseline_ref = "origin/master"
                    merge_git_ref(git_ref="origin/master", verbose=True, dry_run=self._dry_run)
                else:
                    self._baseline_ref = get_common_ancestor("origin/master")
                    # Prefer a symbolic ref if possible
                    if self._baseline_ref is None or self._baseline_ref == get_current_commit(commit="origin/master"):
                        self._baseline_ref = "origin/master"
            self._must_generate_baselines = True

            self._baseline_dir = pathlib.Path(default_baselines_root_dir).absolute()

        else:
            # We treat the "AUTO" string as a request for automatic baseline dir.
            if self._baseline_dir == "AUTO":
                self._baseline_dir = get_mach_baseline_root_dir(self._machine)

            self._baseline_dir = pathlib.Path(self._baseline_dir).absolute()

            # Make sure the baseline root directory exists
            expect(self._baseline_dir.is_dir(), "Baseline_dir {} is not a dir".format(self._baseline_dir))

            if self._integration_test:
                self._baseline_ref = "origin/master"
                merge_git_ref(git_ref=self._baseline_ref, verbose=True, dry_run=self._dry_run)
            else:
                for test in self._tests:
                    test_baseline_dir = self.get_preexisting_baseline(test)
                    expect(test_baseline_dir.is_dir(), "Missing baseline {}".format(test_baseline_dir))

        # Name of the file used to store/check the git sha of the repo used to generate baselines,
        # and name of the file used to store/check the builds for which baselines are available
        # Store it once to avoid typos-like bugs
        self._baseline_sha_file   = pathlib.Path(self._baseline_dir, "baseline_git_sha")
        self._baseline_names_file = pathlib.Path(self._baseline_dir, "baseline_names")

        if self._integration_test:
            master_sha = get_current_commit(commit=self._baseline_ref)
            if not self.baselines_are_present():
                print ("Some baselines were not found. Rebuilding them.")
                self._must_generate_baselines = True
            elif self.baselines_are_expired(expected_baseline_sha=master_sha):
                print ("Baselines expired. Rebuilding them.")
                self._must_generate_baselines = True
            else:
                print ("Baselines found and not expired. Skipping baselines generation.")

        if self._must_generate_baselines:
            print("Using commit {} to generate baselines".format(self._baseline_ref))

        ##################################################
        #   Deduce how many testing resources per test   #
        ##################################################

        if ctest_parallel_level > 0:
            ctest_max_jobs = ctest_parallel_level
            print("Note: honoring requested value for ctest parallel level: {}".format(ctest_max_jobs))
        elif "CTEST_PARALLEL_LEVEL" in os.environ:
            ctest_max_jobs = int(os.environ["CTEST_PARALLEL_LEVEL"])
            print("Note: honoring environment value for ctest parallel level: {}".format(ctest_max_jobs))
        else:
            ctest_max_jobs = get_mach_testing_resources(self._machine)
            print("Note: no value passed for --ctest-parallel-level. Using the default for this machine: {}".format(ctest_max_jobs))

        self._testing_res_count = {"dbg" : ctest_max_jobs,
                                   "sp"  : ctest_max_jobs,
                                   "fpe" : ctest_max_jobs}

        # Deduce how many compilation resources per test
        if make_parallel_level > 0:
            make_max_jobs = make_parallel_level
            print("Note: honoring requested value for make parallel level: {}".format(make_max_jobs))
        else:
            make_max_jobs = get_mach_compilation_resources(self._machine)
            print("Note: no value passed for --make-parallel-level. Using the default for this machine: {}".format(make_max_jobs))

        self._compile_res_count = {"dbg" : make_max_jobs,
                                   "sp"  : make_max_jobs,
                                   "fpe" : make_max_jobs}

        if self._parallel:
            # We need to be aware that other builds may be running too.
            # (Do not oversubscribe the machine)
            make_remainder = make_max_jobs % len(self._tests)
            make_count     = make_max_jobs // len(self._tests)
            ctest_remainder = ctest_max_jobs % len(self._tests)
            ctest_count     = ctest_max_jobs // len(self._tests)

            # In case we have more items in self._tests than cores/gpus (unlikely)
            if make_count == 0:
                make_count = 1
            if ctest_count == 0:
                ctest_count = 1

            for test in self._tests:
                self._compile_res_count[test] = make_count
                if self._tests.index(test)<make_remainder:
                    self._compile_res_count[test] = make_count + 1

                self._testing_res_count[test] = ctest_count
                if self._tests.index(test)<ctest_remainder:
                    self._testing_res_count[test] = ctest_count + 1

                print("test {} can use {} jobs to compile, and {} jobs for testing".format(test,self._compile_res_count[test],self._testing_res_count[test]))

        if self._keep_tree:
            expect(not self._integration_test, "Should not be doing keep-tree with integration testing")
            print("WARNING! You have uncommitted changes in your repo.",
                  "         The PASS/FAIL status may depend on these changes",
                  "         so if you want to keep them, don't forget to create a commit.",sep="\n")
            if self._baseline_dir is None:
                # Make sure the baseline ref is HEAD
                expect(self._baseline_ref == "HEAD",
                       "The option --keep-tree is only available when testing against pre-built baselines "
                       "(--baseline-dir) or HEAD (-b HEAD)")
        else:
            expect(self._dry_run or is_repo_clean(),
                   "Repo must be clean before running. If testing against HEAD or pre-built baselines, "
                   "you can pass `--keep-tree` to allow non-clean repo.")

    ###############################################################################
    def get_baseline_dir(self, test):
    ###############################################################################
        return pathlib.Path(self._baseline_dir, self._test_full_names[test])

    ###############################################################################
    def get_test_dir(self, test):
    ###############################################################################
        return pathlib.Path(self._work_dir, self._test_full_names[test])

    ###############################################################################
    def get_preexisting_baseline(self, test):
    ###############################################################################
        expect(self._baseline_dir is not None, "Cannot supply preexisting baseline without baseline dir")
        return pathlib.Path(self._baseline_dir, self._test_full_names[test], "data")

    ###############################################################################
    def baselines_are_present (self):
    ###############################################################################
        # Check that all baselines are present (one subdir for all values of self._tests)

        # Sanity check
        expect(self._baseline_dir is not None, "Error! This routine should only be called when testing against pre-existing baselines.")

        # Even if a single baseline is missing, we consider all the baselines not present
        for test in self._tests:
            test_baseline_dir = pathlib.Path(self._baseline_dir, self._test_full_names[test], "data")
            if not test_baseline_dir.is_dir():
                return False

        # Note: inside this script we don't know what kind of file should be in the baseline dirs.
        #       If the actual files are missing, some other part of the testing will crash.
        return True

    ###############################################################################
    def baselines_are_expired (self, expected_baseline_sha):
    ###############################################################################
        # Baselines are expired if either:
        #  2) there is no file in baseline_dir containing the sha of the baselines
        #  3) the baselines sha does not match the one passed to this function

        # Sanity check
        expect(self._baseline_dir is not None, "Error! This routine should only be called when testing against pre-existing baselines.")

        # The file specifying what baselines were built during last baselines generation msut be there
        if not self._baseline_names_file.exists():
            return True

        # It might happen that we generate baselines for all build types, then later on
        # for some reason we manually generate baselines for only one build type. The other
        # baselines will still be there, but may be expired. Therefore, we check the
        # baselines_names file, to see what baselines were built last time. If all the
        # baselines we need are there, then we're good
        valid_baselines = run_cmd_no_fail("cat {}".format(self._baseline_names_file.resolve()))
        for test in self._tests:
            if not test in valid_baselines:
                return True

        # No sha file => baselines expired
        if not self._baseline_sha_file.exists():
            return True

        # Different sha => baselines expired
        baseline_sha = run_cmd_no_fail("cat {}".format(self._baseline_sha_file))
        return expected_baseline_sha != baseline_sha

    ###############################################################################
    def get_machine_file(self):
    ###############################################################################
        if self._local:
            return pathlib.Path("~/.cime/scream_mach_file.cmake").expanduser()
        else:
            return pathlib.Path(self._root_dir, "cmake", "machine-files", "{}.cmake".format(self._machine))

    ###############################################################################
    def generate_cmake_config(self, extra_configs, for_ctest=False):
    ###############################################################################

        # Ctest only needs config options, and doesn't need the leading 'cmake '
        result  = "{}-C {}".format("" if for_ctest else "cmake ", self.get_machine_file())

        # Netcdf should be available. But if the user is doing a testing session
        # where all netcdf-related code is disabled, he/she should be able to run
        # even if no netcdf is available
        stat, f_path, err = run_cmd("nf-config --prefix")
        if stat == 0:
            result += " -DNetCDF_Fortran_PATHS={}".format(f_path)
        stat, c_path, err = run_cmd("nc-config --prefix")
        if stat == 0:
            result += " -DNetCDF_C_PATHS={}".format(c_path)

        # Test-specific cmake options
        for key, value in extra_configs:
            result += " -D{}={}".format(key, value)

        # User-requested config options
        custom_opts_keys = []
        for custom_opt in self._custom_cmake_opts:
            expect ("=" in custom_opt, "Error! Syntax error in custom cmake options. Should be `VAR_NAME=VALUE`.")
            if "=" in custom_opt:
                name, value = custom_opt.split("=", 1)
                # Some effort is needed to ensure quotes are perserved
                result += " -D{}='{}'".format(name, value)
                custom_opts_keys.append(name)

        # Common config options (unless already specified by the user)
        if "CMAKE_CXX_COMPILER" not in custom_opts_keys:
            result += " -DCMAKE_CXX_COMPILER={}".format(self._cxx_compiler)
        if "CMAKE_C_COMPILER" not in custom_opts_keys:
            result += " -DCMAKE_C_COMPILER={}".format(self._c_compiler)
        if "CMAKE_Fortran_COMPILER" not in custom_opts_keys:
            result += " -DCMAKE_Fortran_COMPILER={}".format(self._f90_compiler)

        if "SCREAM_DYNAMICS_DYCORE" not in custom_opts_keys:
            result += " -DSCREAM_DYNAMICS_DYCORE=HOMME"

        return result

    ###############################################################################
    def get_taskset_id(self, test):
    ###############################################################################
        # Note: we need to loop through the whole list, since the compile_res_count
        #       might not be the same for all test.

        it = itertools.takewhile(lambda name: name!=test, self._tests)
        offset = sum(self._compile_res_count[prevs] for prevs in it)

        start = offset
        end   = offset + self._compile_res_count[test] - 1

        return start, end

    ###############################################################################
    def create_ctest_resource_file(self, test, build_dir):
    ###############################################################################
        # Create a json file in the test build dir, which ctest will then use
        # to schedule tests in parallel.
        # In the resource file, we have N res groups with 1 slot, with N being
        # what's in self._testing_res_count[test]. On CPU machines, res groups
        # are cores, on GPU machines, res groups are GPUs. In other words, a
        # res group is where we usually bind an individual MPI rank.
        # The id of the res groups on is offset so that it is unique across all builds

        # Note: on CPU, the actual ids are pointless, since ctest job is already places
        # on correct cores with taskset. On GPU, however, these ids are used to select
        # the kokkos device where the tests are run on.

        if self._parallel:
            it = itertools.takewhile(lambda name: name!=test, self._tests)
            start = sum(self._testing_res_count[prevs] for prevs in it)
            end   = start + self._testing_res_count[test]
        else:
            start = 0
            end   = self._testing_res_count[test]

        data = {}

        # This is the only version numbering supported by ctest, so far
        data['version'] = {"major":1,"minor":0}

        devices = []
        for res_id in range(start,end):
            devices.append({"id":str(res_id)})

        # Add resource groups
        data['local'] = [{"devices":devices}]

        with open("{}/ctest_resource_file.json".format(build_dir),'w') as outfile:
            json.dump(data,outfile,indent=2)

    ###############################################################################
    def generate_ctest_config(self, cmake_config, extra_configs, test):
    ###############################################################################
        name = self._test_full_names[test]
        result = ""
        if self._submit:
            result += "CIME_MACHINE={} ".format(self._machine)

        test_dir = self.get_test_dir(test)
        self.create_ctest_resource_file(test,test_dir)

        result += "SCREAM_BUILD_PARALLEL_LEVEL={} CTEST_PARALLEL_LEVEL={} ctest -V --output-on-failure ".format(self._compile_res_count[test], self._testing_res_count[test])
        result += "--resource-spec-file {}/ctest_resource_file.json ".format(test_dir)

        if self._baseline_dir is not None:
            cmake_config += " -DSCREAM_TEST_DATA_DIR={}".format(self.get_preexisting_baseline(test))

        if not self._submit:
            result += "-DNO_SUBMIT=True "

        for key, value in extra_configs:
            result += "-D{}={} ".format(key, value)

        work_dir = pathlib.Path(self._work_dir).joinpath(name)
        result += "-DBUILD_WORK_DIR={} ".format(work_dir)
        result += "-DBUILD_NAME_MOD={} ".format(name)
        result += '-S {}/cmake/ctest_script.cmake -DCMAKE_COMMAND="{}" '.format(self._root_dir, cmake_config)

        # Ctest can only competently manage test pinning across a single instance of ctest. For
        # multiple concurrent instances of ctest, we have to help it.
        if self._parallel:
            start, end = self.get_taskset_id(test)
            result = "taskset -c {}-{} sh -c '{}'".format(start,end,result)

        return result

    ###############################################################################
    def generate_baselines(self, test):
    ###############################################################################
        test_dir = self.get_baseline_dir(test)

        cmake_config = self.generate_cmake_config(self._tests_cmake_args[test])
        cmake_config += " -DSCREAM_BASELINES_ONLY=ON"

        print("===============================================================================")
        print("Generating baseline for test {} with config '{}'".format(self._test_full_names[test], cmake_config))
        print("===============================================================================")

        # We cannot just crash if we fail to generate baselines, since we would
        # not get a dashboard report if we did that. Instead, just ensure there is
        # no baseline file to compare against if there's a problem.
        stat, _, err = run_cmd("{} {}".format(cmake_config, self._root_dir), from_dir=test_dir, verbose=True, dry_run=self._dry_run)
        if stat!= 0:
            print ("WARNING: Failed to configure baselines:\n{}".format(err))
            return False

        cmd = "make -j{} && make -j{} baseline".format(self._compile_res_count[test],self._testing_res_count[test])
        if self._parallel:
            start, end = self.get_taskset_id(test)
            cmd = "taskset -c {}-{} sh -c '{}'".format(start,end,cmd)

        stat, _, err = run_cmd(cmd, from_dir=test_dir, verbose=True, dry_run=self._dry_run)

        if stat != 0:
            print("WARNING: Failed to create baselines:\n{}".format(err))
            return False
        else:
            # Clean up the directory, by removing everything but the 'data' subfolder
            run_cmd_no_fail(r"find -maxdepth 1 -not -name data ! -path . -exec rm -rf {} \;",
                            from_dir=test_dir, verbose=True, dry_run=self._dry_run)

        return True

    ###############################################################################
    def generate_all_baselines(self):
    ###############################################################################
        git_head_ref = get_current_head()

        print("###############################################################################")
        print("Generating baselines for ref {}".format(self._baseline_ref))
        print("###############################################################################")

        # First, create build directories (one per test)
        for test in self._tests:
            test_dir = self.get_baseline_dir(test)

            # Create this test's build dir
            if test_dir.exists():
                shutil.rmtree(str(test_dir))

            test_dir.mkdir(parents=True)

        checkout_git_ref(self._baseline_ref, verbose=True, dry_run=self._dry_run)

        success = True
        num_workers = len(self._tests) if self._parallel else 1
        with threading3.ProcessPoolExecutor(max_workers=num_workers) as executor:

            future_to_test = {
                executor.submit(self.generate_baselines, test) : test
                for test in self._tests}

            for future in threading3.as_completed(future_to_test):
                test = future_to_test[future]
                success &= future.result()

                if not success and self._fast_fail:
                    print('Generation of baselines for build {} failed'.format(self._test_full_names[test]))
                    return False

        if success:
            # Store the sha used for baselines generation
            run_cmd_no_fail("echo '{}' > {}".format(get_current_commit(commit=self._baseline_ref),self._baseline_sha_file))
            # Store the name of the builds for which we created a baseline
            tmp_string = ""
            for test in self._tests:
                tmp_string += " {}".format(test)
            run_cmd_no_fail("echo '{}' > {}".format(tmp_string,self._baseline_names_file))

        checkout_git_ref(git_head_ref, verbose=True, dry_run=self._dry_run)

        return success

    ###############################################################################
    def run_test(self, test):
    ###############################################################################
        git_head = get_current_head()

        print("===============================================================================")
        print("Testing '{}' for test '{}'".format(git_head, self._test_full_names[test]))
        print("===============================================================================")

        test_dir = self.get_test_dir(test)
        cmake_config = self.generate_cmake_config(self._tests_cmake_args[test], for_ctest=True)
        ctest_config = self.generate_ctest_config(cmake_config, [], test)

        if self._quick_rerun and pathlib.Path("{}/CMakeCache.txt".format(test_dir)).is_file():
            # Do not purge bld dir, and do not rerun config step.
            # Note: make will still rerun cmake if some cmake file has changed
            ctest_config += "-DSKIP_CONFIG_STEP=TRUE "
            if self._quick_rerun_failed:
                ctest_config += "--rerun-failed "
        else:
            # This directory might have been used also to build the model to generate baselines.
            # Although it's ok to build in the same dir, we MUST make sure to erase cmake's cache
            # and internal files from the previous build (CMakeCache.txt and CMakeFiles folder)
            run_cmd_no_fail("rm -rf CMake*", from_dir=test_dir, dry_run=self._dry_run)

        success = run_cmd(ctest_config, from_dir=test_dir, arg_stdout=None, arg_stderr=None, verbose=True, dry_run=self._dry_run)[0] == 0

        return success

    ###############################################################################
    def run_all_tests(self):
    ###############################################################################
        print("###############################################################################")
        print("Running tests!")
        print("###############################################################################")

        # First, create build directories (one per test)
        for test in self._tests:
            test_dir = self.get_test_dir(test)

            # Create this test's build dir
            if test_dir.exists():
                if not self._quick_rerun:
                    shutil.rmtree(str(test_dir))
                    test_dir.mkdir(parents=True)
            else:
                test_dir.mkdir(parents=True)

        success = True
        tests_success = {
            test : False
            for test in self._tests}

        num_workers = len(self._tests) if self._parallel else 1
        with threading3.ProcessPoolExecutor(max_workers=num_workers) as executor:
            future_to_test = {
                executor.submit(self.run_test,test) : test
                for test in self._tests}

            for future in threading3.as_completed(future_to_test):
                test = future_to_test[future]
                tests_success[test] = future.result()
                success &= tests_success[test]
                # If failed, and fast fail is requested, return immediately
                # Note: this is effective only if num_worksers=1
                if not success and self._fast_fail:
                    break

        for t,s in tests_success.items():
            if not s:
                print("Build type {} failed. Here's a list of failed tests:".format(self._test_full_names[t]))
                test_dir = self.get_test_dir(t)
                test_results_dir = pathlib.Path(test_dir, "Testing", "Temporary")
                for result in test_results_dir.glob("LastTestsFailed*"):
                    print(result.read_text())

        return success

    ###############################################################################
    def test_all_scream(self):
    ###############################################################################

        # Add any override the user may have requested
        for env_var in self._custom_env_vars:
            key,val = env_var.split("=",2)
            os.environ.update( { key : val } )

        success = True
        try:
            # If needed, generate baselines first
            if self._must_generate_baselines:
                expect(self._baseline_ref is not None, "Missing baseline ref")

                success = self.generate_all_baselines()
                if not success:
                    print ("Error(s) occurred during baselines generation phase")
                    return False

            # If requested, run tests
            if self._perform_tests:
                success &= self.run_all_tests()
                if not success:
                    print ("Error(s) occurred during test phase")

        finally:
            if not self._keep_tree:
                # Cleanup the repo if needed
                cleanup_repo(self._original_branch, self._original_commit, dry_run=self._dry_run)

        return success
