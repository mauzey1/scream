#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"

#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestVdShocSolve {

  static void run_bfb()
  {
    VdShocDecompData f90_decomp_data[] = {
      VdShocDecompData(1, 3, 4, 1),
    };

    VdShocSolveData f90_solve_data[] = {
      VdShocSolveData(f90_decomp_data[0].shcol(), f90_decomp_data[0].nlev()),
    };

    static constexpr Int num_runs = sizeof(f90_decomp_data) / sizeof(VdShocDecompData);

    // Generate random input data
    for (Int i = 0; i < num_runs; ++i) {
      VdShocDecompData& d_f90_decomp = f90_decomp_data[i];
      VdShocSolveData& d_f90_solve = f90_solve_data[i];
      d_f90_decomp.randomize();
      d_f90_solve.randomize();
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    VdShocDecompData cxx_decomp_data[] = {
      VdShocDecompData(f90_decomp_data[0])
    };
    VdShocSolveData cxx_solve_data[] = {
      VdShocSolveData(f90_solve_data[0])
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (Int i = 0; i < num_runs; ++i) {
      VdShocDecompData& d_f90_decomp = f90_decomp_data[i];
      VdShocSolveData& d_f90_solve = f90_solve_data[i];

      // expects data in C layout
      vd_shoc_decomp(d_f90_decomp);

      // Copy decomp data to solver
      for (Int k = 0; k < d_f90_decomp.total(d_f90_decomp.dl); ++k) {
        d_f90_solve.dl[k] = d_f90_decomp.dl[k];
        d_f90_solve.du[k] = d_f90_decomp.du[k];
        d_f90_solve.d [k] = d_f90_decomp.d [k];
      }
      vd_shoc_solve(d_f90_solve);
    }

    // Get data from cxx
    for (Int i = 0; i < num_runs; ++i) {
      VdShocDecompData& d_cxx_decomp = cxx_decomp_data[i];
      VdShocSolveData& d_cxx_solve = cxx_solve_data[i];

      d_cxx_decomp.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      vd_shoc_decomp_f(d_cxx_decomp.shcol(), d_cxx_decomp.nlev(), d_cxx_decomp.nlevi(),
                       d_cxx_decomp.kv_term, d_cxx_decomp.tmpi, d_cxx_decomp.rdp_zt,
                       d_cxx_decomp.dtime, d_cxx_decomp.flux, d_cxx_decomp.du, d_cxx_decomp.dl,
                       d_cxx_decomp.d);
      d_cxx_decomp.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout

      // Copy decomp data to solver
      for (Int k = 0; k < d_cxx_decomp.total(d_cxx_decomp.dl); ++k) {
        d_cxx_solve.dl[k] = d_cxx_decomp.dl[k];
        d_cxx_solve.du[k] = d_cxx_decomp.du[k];
        d_cxx_solve.d [k] = d_cxx_decomp.d [k];
      }

      d_cxx_solve.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      vd_shoc_solve_f(d_cxx_solve.shcol(), d_cxx_solve.nlev(), d_cxx_solve.du, d_cxx_solve.dl,
                      d_cxx_solve.d, d_cxx_solve.var);
      d_cxx_solve.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      VdShocSolveData& d_f90 = f90_solve_data[i];
      VdShocSolveData& d_cxx = cxx_solve_data[i];
      for (Int k = 0; k < d_f90.total(d_f90.var); ++k) {
        REQUIRE(d_f90.total(d_f90.var) == d_cxx.total(d_cxx.var));
        REQUIRE(d_f90.var[k] == d_cxx.var[k]);
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("vd_shoc_solve_bfb", "[shoc]")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestVdShocSolve;

  TestStruct::run_bfb();
}

} // empty namespace
