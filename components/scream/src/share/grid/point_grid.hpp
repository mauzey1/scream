#ifndef SCREAM_POINT_GRID_HPP
#define SCREAM_POINT_GRID_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/scream_types.hpp"

#include "ekat/mpi/ekat_comm.hpp"

namespace scream
{

/*
 * A grid consisting of a bunch of points
 *
 * A point grid simply stores a bunch of dofs gids. Unlike SEGrid,
 * which also stores info on the SE element each dof belongs to,
 * this class does not store any additional info regarding its dofs.
 * In particular, the map lid->idx is an identity, and the native
 * layout has only one tag: Column.
 *
 * This grid is typical of Physics parametrizations, and is also used
 * to interface to the component coupler.
 */

class PointGrid : public AbstractGrid
{
public:

  // The dofs lat-lon coordinates, and the area associated with a dof
  using geo_view_2d = kokkos_types::view_2d<double>;
  using geo_view_1d = kokkos_types::view_1d<double>;

  PointGrid (const std::string& grid_name,
             const dofs_list_type& dofs_gids,
             const int num_global_cols,
             const int num_vertical_levels);
  virtual ~PointGrid () = default;

  // Native layout of a dof. This is the natural way to index a dof in the grid.
  // E.g., for a 2d structured grid, this could be a set of 2 indices.
  FieldLayout get_native_dof_layout () const override;

  // Get dofs geo info
  const geo_view_1d& get_dofs_area   () const { return m_dofs_area; }
  const geo_view_2d& get_dofs_coords () const { return m_dofs_coords; }

  void set_geometry_data (const geo_view_2d& dofs_coords,
                          const geo_view_1d& dofs_area);

private:

  // Geometric info
  geo_view_2d      m_dofs_coords;
  geo_view_1d      m_dofs_area;
};

// Create a point grid, with linear range of gids, evenly partitioned
// among the ranks in the given communicator.
PointGrid
create_point_grid (const std::string& name,
                   const int num_global_cols,
                   const int num_vertical_lev,
                   const ekat::Comm& comm);

} // namespace scream

#endif // SCREAM_POINT_GRID_HPP
