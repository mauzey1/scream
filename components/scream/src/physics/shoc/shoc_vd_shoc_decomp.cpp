#include "shoc_vd_shoc_decomp_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing vd_shoc_decomp on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
