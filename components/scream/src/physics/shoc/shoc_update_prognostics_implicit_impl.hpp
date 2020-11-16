#ifndef SHOC_UPDATE_PROGNOSTICS_IMPLICIT_IMPL_HPP
#define SHOC_UPDATE_PROGNOSTICS_IMPLICIT_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc update_prognostics_implicit. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
KOKKOS_FUNCTION
void Functions<S,D>::update_prognostics_implicit(
  const MemberType&            team,
  const Int&                   nlev,
  const Int&                   nlevi,
  const Int&                   num_tracer,
  const Scalar&                dtime,
  const uview_1d<const Spack>& dz_zt,
  const uview_1d<const Spack>& dz_zi,
  const uview_1d<const Spack>& rho_zt,
  const uview_1d<const Spack>& zt_grid,
  const uview_1d<const Spack>& zi_grid,
  const uview_1d<const Spack>& tk,
  const uview_1d<const Spack>& tkh,
  const Scalar&                uw_sfc,
  const Scalar&                vw_sfc,
  const Scalar&                wthl_sfc,
  const Scalar&                wqw_sfc,
  const uview_1d<const Spack>& wtracer_sfc,
  const uview_1d<Spack>&       rdp_zt,
  const uview_1d<Spack>&       tmpi,
  const uview_1d<Spack>&       tkh_zi,
  const uview_1d<Spack>&       tk_zi,
  const uview_1d<Spack>&       rho_zi,
  const uview_1d<Spack>&       thetal,
  const uview_1d<Spack>&       qw,
  const uview_2d<Spack>&       tracer,
  const uview_1d<Spack>&       tke,
  const uview_1d<Spack>&       u_wind,
  const uview_1d<Spack>&       v_wind)
{
  const auto nlev_packs  = ekat::npack<Spack>(nlev);
  const auto nlevi_packs = ekat::npack<Spack>(nlevi);

  const auto last_nlev_pack = nlev-1;
  const auto last_nlev_indx = (nlev-1)%Spack::n;
  const auto last_nlevi_pack = nlevi-1;
  const auto last_nlevi_indx = (nlevi-1)%Spack::n;

  // linearly interpolate tkh, tk, and air density onto the interface grids
  linear_interp(team,zt_grid,zi_grid,tkh,tkh_zi,nlev,nlevi,0);
  linear_interp(team,zt_grid,zi_grid,tk,tk_zi,nlev,nlevi,0);
  linear_interp(team,zt_grid,zi_grid,rho_zt,rho_zi,nlev,nlevi,0);

  // Define the tmpi variable, which is really dt*(g*rho)**2/dp
  // at interfaces. Substitue dp = g*rho*dz in the above equation
  compute_tmpi(team, nlevi, dtime, rho_zi, dz_zi, tmpi);

  // compute 1/dp term, needed in diffusion solver
  dp_inverse(team, nlev, rho_zt, dz_zt, rdp_zt);

  // compute terms needed for the implicit surface stress (ksrf)
  // and tke flux calc (wtke_sfc)
  Scalar ksrf, wtke_sfc;
  {
    const auto wsmin = 1;
    const auto ksrfmin = 1e-4;
    const auto ustarmin = 0.01;

    const auto rho = rho_zi(last_nlevi_pack)[last_nlevi_indx];
    const auto uw = uw_sfc;
    const auto vw = vw_sfc;

    const auto taux = rho*uw;
    const auto tauy = rho*vw;

    const auto u_wind_sfc = u_wind(last_nlev_pack)[last_nlev_indx];
    const auto v_wind_sfc = v_wind(last_nlev_pack)[last_nlev_indx];

    const auto ws = std::max(std::sqrt((u_wind_sfc*u_wind_sfc) + v_wind_sfc*v_wind_sfc), sp(wsmin));
    const auto tau = std::sqrt(taux*taux + tauy*tauy);
    ksrf = std::max(tau/ws, sp(ksrfmin));

    const auto ustar = std::max(std::sqrt(std::sqrt(uw*uw + vw*vw)), sp(ustarmin));
    wtke_sfc = ustar*ustar*ustar;
  }

  // compute surface fluxes for liq. potential temp, water and tke
  {
    const auto rho_zi_sfc = rho_zi(last_nlevi_pack)[last_nlevi_indx];
    const auto rdp_zt_sfc = rdp_zt(last_nlev_pack)[last_nlev_indx];
    auto thetal_sfc = thetal(last_nlev_pack)[last_nlev_indx];
    auto qw_sfc = qw(last_nlev_pack)[last_nlev_indx];
    auto tke_sfc = tke(last_nlev_pack)[last_nlev_indx];
    auto tracer_sfc = Kokkos::subview(tracer, nlev-1, Kokkos::ALL());

    sfc_fluxes(team, num_tracer, dtime, rho_zi_sfc, rdp_zt_sfc,
               wthl_sfc, wqw_sfc, wtke_sfc, wtracer_sfc,
               thetal_sfc, qw_sfc, tke_sfc, tracer_sfc);

    thetal(last_nlev_pack)[last_nlev_indx] = thetal_sfc;
    qw(last_nlev_pack)[last_nlev_indx] = qw_sfc;
    tke(last_nlev_pack)[last_nlev_indx] = tke_sfc;
  }



}

} // namespace shoc
} // namespace scream

#endif
