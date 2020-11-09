#include "share/grid/point_grid.hpp"

#include <numeric>

namespace scream {

PointGrid::
PointGrid (const std::string&     grid_name,
           const dofs_list_type&  dofs_gids,
           const int              num_global_cols,
           const int              num_vertical_levels)
 : AbstractGrid(GridType::Point,grid_name)
{
  m_num_global_dofs = num_global_cols;
  m_num_local_dofs  = dofs_gids.size();
  m_num_vert_levs   = num_vertical_levels;

  // The lid->idx map is the identity map.
  m_lid_to_idx = decltype(m_lid_to_idx)("lid to idx",m_num_local_dofs,1);

  auto h_lid_to_idx = Kokkos::create_mirror_view(m_lid_to_idx);
  std::iota(h_lid_to_idx.data(),h_lid_to_idx.data()+m_num_local_dofs,0);
  Kokkos::deep_copy(m_lid_to_idx,h_lid_to_idx);
}

FieldLayout
PointGrid::get_native_dof_layout () const
{
  using namespace ShortFieldTagsNames;
  return FieldLayout({COL},{m_num_local_dofs});
}

PointGrid
create_point_grid (const std::string& grid_name,
                   const int num_global_cols,
                   const int num_vertical_lev,
                   const ekat::Comm& comm)
{
  // Compute how many columns are owned by this rank
  const int num_procs = comm.size();

  auto num_my_cols = num_global_cols / num_procs;
  int remainder   = num_global_cols % num_procs;
  int dof_offset  = num_my_cols*comm.rank();
  if (comm.rank() < remainder) {
    ++num_my_cols;
    dof_offset += comm.rank();
  } else {
    dof_offset += remainder;
  }

  PointGrid::dofs_list_type dofs_gids ("phys dofs",num_my_cols);
  auto h_dofs_gids = Kokkos::create_mirror_view(dofs_gids);
  std::iota(h_dofs_gids.data(),h_dofs_gids.data()+num_my_cols,dof_offset);
  Kokkos::deep_copy(dofs_gids,h_dofs_gids);

  return PointGrid(grid_name,dofs_gids,num_global_cols,num_vertical_lev);
}

} // namespace scream
