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
  Scalar ksrf;
  {
    // Minimum wind speed for ksrfturb computation [ m/s ]
    const auto wsmin = 1;

    // Minimum surface drag coefficient  [ kg/s/m^2 ]
    const auto ksrfmin = 1e-4;

    const auto rho = rho_zi(nlevi_packs-1)[(nlevi-1)%Spack::n];
    const auto uw = uw_sfc;
    const auto vw = vw_sfc;

    const auto taux = rho*uw;
    const auto tauy = rho*vw;

    const auto u_wind_sfc = u_wind(nlev_packs-1)[(nlev-1)%Spack::n];
    const auto v_wind_sfc = v_wind(nlev_packs-1)[(nlev-1)%Spack::n];

    // compute the wind speed
    const auto ws = std::max(std::sqrt((u_wind_sfc*u_wind_sfc) + v_wind_sfc*v_wind_sfc), wsmin);
    const auto tau = std::sqrt(taux*taux + tauy*tauy);
    ksrf = std::max(tau/ws, ksrfmin);
  }

  // compute term needed for tke flux calc (wtke_sfc)
//  wtke_sfc(1:shcol) = tke_srf_flux_term(shcol, uw_sfc, vw_sfc)

//  ! compute surface fluxes for liq. potential temp, water and tke
//  call sfc_fluxes(shcol, num_tracer, dtime, rho_zi(:,nlevi), rdp_zt(:,nlev), &
//                  wthl_sfc, wqw_sfc, wtke_sfc, wtracer_sfc, &
//                  thetal(:,nlev), qw(:,nlev), tke(:,nlev), tracer(:,nlev,:))


}

} // namespace shoc
} // namespace scream

#endif
