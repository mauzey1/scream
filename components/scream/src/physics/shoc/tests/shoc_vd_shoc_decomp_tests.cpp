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
struct UnitWrap::UnitTest<D>::TestVdShocDecomp {

  static void run_bfb()
  {
    VdShocDecompData f90_data[] = {
      VdShocDecompData(10, 71, 72, 1),
      VdShocDecompData(10, 12, 13, 0),
      VdShocDecompData(7,  16, 17, 100),
      VdShocDecompData(2, 7, 8, 7.5),
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(VdShocDecompData);

    // Generate random input data
    for (auto& d : f90_data) {
      d.randomize();
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    VdShocDecompData cxx_data[] = {
      VdShocDecompData(f90_data[0]),
      VdShocDecompData(f90_data[1]),
      VdShocDecompData(f90_data[2]),
      VdShocDecompData(f90_data[3])
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      vd_shoc_decomp(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      vd_shoc_decomp_f(d.shcol(), d.nlev(), d.nlevi(), d.kv_term, d.tmpi, d.rdp_zt, d.dtime, d.flux, d.du, d.dl, d.d);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      VdShocDecompData& d_f90 = f90_data[i];
      VdShocDecompData& d_cxx = cxx_data[i];
      for (Int k = 0; k < d_f90.total(d_f90.du); ++k) {
        REQUIRE(d_f90.total(d_f90.du) == d_cxx.total(d_cxx.du));
        REQUIRE(d_f90.du[k] == d_cxx.du[k]);
        REQUIRE(d_f90.total(d_f90.du) == d_cxx.total(d_cxx.dl));
        REQUIRE(d_f90.dl[k] == d_cxx.dl[k]);
        REQUIRE(d_f90.total(d_f90.du) == d_cxx.total(d_cxx.d));
        REQUIRE(d_f90.d[k] == d_cxx.d[k]);
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("vd_shoc_decomp_bfb", "[shoc]")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestVdShocDecomp;

  TestStruct::run_bfb();
}

} // empty namespace
