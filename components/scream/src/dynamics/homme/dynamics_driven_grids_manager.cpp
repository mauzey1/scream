#include "dynamics/homme/dynamics_driven_grids_manager.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"
#include "dynamics/homme/physics_dynamics_remapper.hpp"

#include "share/grid/se_grid.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/remap/inverse_remapper.hpp"

// Get all Homme's compile-time dims
#include "homme_dimensions.hpp"

namespace scream
{

DynamicsDrivenGridsManager::DynamicsDrivenGridsManager (const ekat::Comm& /* comm */, const ekat::ParameterList& p)
 : m_params(p)
{
  // Nothing to do here
}

DynamicsDrivenGridsManager::remapper_ptr_type
DynamicsDrivenGridsManager::do_create_remapper (const grid_ptr_type from_grid,
                                                const grid_ptr_type to_grid) const {
  using PDR = PhysicsDynamicsRemapper<remapper_type::real_type>;
  auto pd_remapper = std::make_shared<PDR>(m_grids.at("SE Physics"),m_grids.at("SE Dynamics"));
  if (from_grid->name()=="SE Physics" &&
      to_grid->name()=="SE Dynamics") {
    return pd_remapper;
  } else if (from_grid->name()=="SE Dynamics" &&
             to_grid->name()=="SE Physics") {
    return std::make_shared<InverseRemapper<Real>>(pd_remapper);
  }
  return nullptr;
}

void DynamicsDrivenGridsManager::build_grid (const std::string& grid_name)
{
  // We cannot start to build SCREAM grid structures until homme has inited the geometry,
  // so do that if not already done.
  // NOTE: this *requires* that init_parallel_f90 and init_params_f90 have already been
  //       called. But we can't do it from here, so simply call init_geometry_f90.
  //       if parallel/params have not been called yet, the code will abort.
  if (!is_geometry_inited_f90()) {
    init_geometry_f90();
  }

  if (grid_name=="SE Physics") {
    build_physics_grid();
  } else if (grid_name=="SE Dynamics") {
    build_dynamics_grid();
  }

  if (grid_name==m_params.get<std::string>("Reference Grid")) {
    m_grids["Reference"] = get_grid(grid_name);
  }
}

void DynamicsDrivenGridsManager::build_dynamics_grid () {
  if (m_grids.find("SE Dynamics")==m_grids.end()) {

    // Initialize the dyn grid
    const int nelemd = get_num_owned_elems_f90();
    const int nlev   = get_nlev_f90();
    auto dyn_grid = std::make_shared<SEGrid>("SE Dynamics",nelemd,NP,nlev);

    // Create dynamics dofs map
    AbstractGrid::dofs_list_type      dofs("dyn dofs",nelemd*NP*NP);
    AbstractGrid::lid_to_idx_map_type lids_to_elgpgp("dyn lid to elgpgp",nelemd*NP*NP,3);

    auto h_lids_to_elgpgp = Kokkos::create_mirror_view(lids_to_elgpgp);
    auto h_dofs = Kokkos::create_mirror_view(dofs);

    // Get (ie,igp,jgp,gid) data for each dof
    get_cols_indices_f90(h_dofs.data(),h_lids_to_elgpgp.data(),false);

    Kokkos::deep_copy(dofs,h_dofs);
    Kokkos::deep_copy(lids_to_elgpgp,h_lids_to_elgpgp);

    dyn_grid->set_dofs (dofs, lids_to_elgpgp);

    // Set the grid in the map
    m_grids["SE Dynamics"] = m_grids["Dynamics"] = dyn_grid;
  }
}

void DynamicsDrivenGridsManager::build_physics_grid () {
  if (m_grids.find("Physics")==m_grids.end()) {

    // Initialize the phys grid
    const int nlev = get_nlev_f90();

    // Create the physics dofs map
    const int num_cols = get_num_owned_columns_f90 ();
    AbstractGrid::dofs_list_type dofs("phys dofs",num_cols);
    auto h_dofs = Kokkos::create_mirror_view(dofs);

    // Get gid for each dof
    get_cols_gids_f90(h_dofs.data(), true);

    Kokkos::deep_copy(dofs,h_dofs);

    auto phys_grid = std::make_shared<PointGrid>("Physics",dofs,nlev);

    // Get dofs area/coords
    PointGrid::geo_view_1d d_area("",num_cols);
    PointGrid::geo_view_2d d_coords("",2,num_cols);
    auto h_area   = Kokkos::create_mirror_view(d_area);
    auto h_coords = Kokkos::create_mirror_view(d_coords);

    get_cols_geo_specs_f90(h_coords.data(), h_area.data());
    Kokkos::deep_copy(d_area, h_area);
    Kokkos::deep_copy(d_coords, h_coords);

    phys_grid->set_geometry_data(d_coords, d_area);

    // Set the grid in the map
    m_grids["Physics"] = phys_grid;
  }
}

} // namespace scream
