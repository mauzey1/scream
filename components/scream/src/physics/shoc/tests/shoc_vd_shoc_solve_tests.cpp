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
//    VdShocSolveData f90_data[] = {
//      // TODO
//    };

//    static constexpr Int num_runs = sizeof(f90_data) / sizeof(VdShocSolveData);

//    // Generate random input data
//    for (auto& d : f90_data) {
//      d.randomize();
//    }

//    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
//    // inout data is in original state
//    VdShocSolveData cxx_data[] = {
//      // TODO
//    };

//    // Assume all data is in C layout

//    // Get data from fortran
//    for (auto& d : f90_data) {
//      // expects data in C layout
//      vd_shoc_solve(d);
//    }

//    // Get data from cxx
//    for (auto& d : cxx_data) {
//      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
//      vd_shoc_solve_f(d.shcol(), d.nlev(), d.du, d.dl, d.d, d.var);
//      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
//    }

//    // Verify BFB results, all data should be in C layout
//    for (Int i = 0; i < num_runs; ++i) {
//      VdShocSolveData& d_f90 = f90_data[i];
//      VdShocSolveData& d_cxx = cxx_data[i];
//      for (Int k = 0; k < d_f90.total(d_f90.var); ++k) {
//        REQUIRE(d_f90.total(d_f90.var) == d_cxx.total(d_cxx.var));
//        REQUIRE(d_f90.var[k] == d_cxx.var[k]);
//      }

//    }
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
