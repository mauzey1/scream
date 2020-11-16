#ifndef SHOC_VD_SHOC_DECOMP_IMPL_HPP
#define SHOC_VD_SHOC_DECOMP_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc vd_shoc_decomp. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::vd_shoc_decomp(
  const MemberType&            team,
  const Int&                   nlev,
  const uview_1d<const Spack>& kv_term,
  const uview_1d<const Spack>& tmpi,
  const uview_1d<const Spack>& rdp_zt,
  const Scalar&                dtime,
  const Scalar&                flux,
  const uview_1d<Spack>&       du,
  const uview_1d<Spack>&       dl,
  const uview_1d<Spack>&       d)
{
  const auto ggr = C::gravit;

  const auto skv_term = scalarize(kv_term);
  const auto stmpi = scalarize(tmpi);

  const Int nlev_pack = ekat::npack<Spack>(nlev);

  // Compute entries of the tridiagonal system
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_pack), [&] (const Int& k) {
    // Compute shift of kv_term and tmpi
    Spack kv_term_k, kv_term_kp1, tmpi_k, tmpi_kp1;
    auto range_pack = ekat::range<IntSmallPack>(k*Spack::n);
    auto shift_range = range_pack;
    shift_range.set(range_pack > nlev-1, 1); // don't calculate shift above nlev-1
    ekat::index_and_shift<1>(skv_term, shift_range, kv_term_k, kv_term_kp1);
    ekat::index_and_shift<1>(stmpi, shift_range, tmpi_k, tmpi_kp1);

    // Determine superdiagonal (du) and subdiagonal (dl) coeffs of the
    // tridiagonal diffusion matrix.
    du(k) = -kv_term_kp1*tmpi_kp1*rdp_zt(k);
    dl(k) = -kv_term(k)*tmpi(k)*rdp_zt(k);

    // The bottom element of the superdiagonal (du) and the top element of
    // the subdiagonal (dl) is set to zero (not included in linear system).
    du(k).set(range_pack == nlev-1, 0);
    dl(k).set(range_pack == 0, 0);

    // The diagonal elements are a combination of du and dl (d=1-du-dl). Surface
    // fluxes are applied explicitly in the diagonal at the top level.
    d(k).set(range_pack == 0, 1 - du(k));
    d(k).set(0 < range_pack && range_pack < nlev-1, 1 - du(k) - dl(k));
    d(k).set(range_pack == nlev-1, 1 - dl(k) + flux*dtime*ggr*rdp_zt(k));
  });
}

} // namespace shoc
} // namespace scream

#endif
