#ifndef SHOC_VD_SHOC_SOLVE_IMPL_HPP
#define SHOC_VD_SHOC_SOLVE_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

//struct ShocSolveImpl {

//  template <typename TridiagArray>
//  Kokkos::View<typename TridiagArray::value_type**, Kokkos::LayoutRight>
//  KOKKOS_INLINE_FUNCTION
//  get_diags (const TridiagArray& A, const int diag_idx) {
//    return Kokkos::View<typename TridiagArray::value_type**, Kokkos::LayoutRight>(
//      &A.impl_map().reference(diag_idx, 0, 0),
//      A.extent_int(1), A.extent_int(2));
//  }

//  template <typename DataArray>
//  KOKKOS_INLINE_FUNCTION
//  Kokkos::View<typename DataArray::value_type*>
//  get_x (const DataArray& X) {
//    assert(X.extent_int(1) == 1);
//    return Kokkos::View<typename DataArray::value_type*>(
//      &X.impl_map().reference(0, 0), X.extent_int(0));
//  }

//  template <typename Scalar>
//  using TridiagArray = Kokkos::View<Scalar***, Kokkos::LayoutRight>;
//  template <typename Scalar>
//  using DataArray = Kokkos::View<Scalar**, Kokkos::LayoutRight>;

//  using SHF = Functions<Real, DefaultDevice>;
//  using KT         = typename SHF::KT;
//  using ExeSpace   = typename KT::ExeSpace;
//  using MemberType = typename SHF::MemberType;
//  using Scalar     = typename SHF::Scalar;
//  using Spack      = typename SHF::Spack;

//  TridiagArray A;
//  DataArray X;

//  template <typename ViewType>
//  void init (const ViewType dl, const ViewType d, const ViewType du, const ViewType rhs)
//  {
//    A = TridiagArray("A", 3, dl.extent_int(1), ekat::npack<Spack>(dl.extent_int(0)));
//    X = DataArray("X", dl.extent_int(1), ekat::npack<Spack>(dl.extent_int(0)));
//  }

//  void shoc_solve (const int nrows, const int nprob, TridiagArray A, DataArray X)
//  {


//      const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1);

//      const auto f = KOKKOS_LAMBDA (const MemberType& team) {
//        const auto dl = get_diags(A, 0);
//        const auto d  = get_diags(A, 1);
//        const auto du = get_diags(A, 2);
//        ekat::bfb(team, dl, d, du, X);
//      };
//      Kokkos::parallel_for(policy, f);
//  }

  //template<typename S, typename D>
  //KOKKOS_FUNCTION
  //void Functions<S,D>::vd_shoc_solve(TridiagArray)
  //{

  //}

//  }; // struct ShocSolveImpl

} // namespace shoc
} // namespace scream

#endif
