#ifndef SHOC_VD_SHOC_SOLVE_IMPL_HPP
#define SHOC_VD_SHOC_SOLVE_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc vd_shoc_solve. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::vd_shoc_solve(const Int& shcol, const Int& nlev, const uview_1d<const Spack>& du, const uview_1d<const Spack>& dl, const uview_1d<const Spack>& d, const uview_1d<Spack>& var)
{

}

} // namespace shoc
} // namespace scream

#endif
