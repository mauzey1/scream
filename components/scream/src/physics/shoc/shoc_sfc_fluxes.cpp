#include "shoc_sfc_fluxes_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing sfc_fluxes on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
