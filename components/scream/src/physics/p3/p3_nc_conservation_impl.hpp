#ifndef P3_NC_CONSERVATION_IMPL_HPP
#define P3_NC_CONSERVATION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 nc_conservation. Clients should NOT
 * #include this file, but include p3_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::nc_conservation(const Spack& nc, const Spack& nc_selfcollect_tend, const Spack& dt, Spack& nc_collect_tend, Spack& nc2ni_immers_freeze_tend, Spack& nc_accret_tend, Spack& nc2nr_autoconv_tend)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace p3
} // namespace scream

#endif