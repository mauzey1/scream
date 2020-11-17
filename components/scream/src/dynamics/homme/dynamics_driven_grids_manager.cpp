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
  // Valid names for the dyn grid
  auto& dgn = m_dyn_grid_aliases;

  dgn.insert("Dynamics");
  dgn.insert("SE Dynamics");
  dgn.insert("GLL");
  dgn.insert("Dyn Grid");

  // Valid names for the different phys grids

  // No redistribution of columns
  auto& pg0 = m_phys_grid_aliases[0];
  pg0.insert("Physics");
  pg0.insert("Physics GLL"); // Phys columns are SE gll points

  auto& pg2 = m_phys_grid_aliases[2];
  pg2.insert("Phys Grid");
  pg2.insert("Physics PG2"); // Phys columns are FV points in 2x2 subcells of SE cell

  auto& pg3 = m_phys_grid_aliases[3];
  pg3.insert("Physics PG3"); // Phys columns are FV points in 3x3 subcells of SE cell

  // Twin columns
  auto& pg10 = m_phys_grid_aliases[10];
  pg10.insert("Physics Twin");     // Same as GLL
  pg10.insert("Physics GLL Twin"); // Phys columns are SE gll points

  auto& pg12 = m_phys_grid_aliases[10];
  pg12.insert("Phys Grid Twin");
  pg12.insert("Physics PG2 Twin"); // Phys columns are FV points in 2x2 subcells of SE cell

  auto& pg13 = m_phys_grid_aliases[10];
  pg13.insert("Physics PG3 Twin"); // Phys columns are FV points in 2x2 subcells of SE cell

  // TODO: add other rebalancing?

  for (const auto& gn : dgn) {
    m_valid_grid_names.insert(gn);
  }
  for (const auto& it : m_phys_grid_aliases) {
    for (const auto& gn : it.second) {
      m_valid_grid_names.insert(gn);
    }
  }
}

DynamicsDrivenGridsManager::remapper_ptr_type
DynamicsDrivenGridsManager::do_create_remapper (const grid_ptr_type from_grid,
                                                const grid_ptr_type to_grid) const {
  using PDR = PhysicsDynamicsRemapper<remapper_type::real_type>;
  auto pd_remapper = std::make_shared<PDR>(m_grids.at("Physics"),m_grids.at("Dynamics"));
  if (from_grid->name()=="SE Physics" &&
      to_grid->name()=="SE Dynamics") {
    return pd_remapper;
  } else if (from_grid->name()=="SE Dynamics" &&
             to_grid->name()=="SE Physics") {
    return std::make_shared<InverseRemapper<Real>>(pd_remapper);
  }
  return nullptr;
}

void DynamicsDrivenGridsManager::
build_grids (const std::set<std::string>& grid_names,
             const std::string& reference_grid) {
  // Sanity check first
  for (const auto& gn : grid_names) {
    EKAT_REQUIRE_MSG (supported_grids().count(gn)==1,
                      "Error! Grid '" + gn + "' is not supported by this grid manager.\n");
  }
  // We cannot start to build SCREAM grid structures until homme has inited the geometry,
  // so do that if not already done.
  // NOTE: this *requires* that init_parallel_f90 and init_params_f90 have already been
  //       called. But we can't do it from here, so simply call init_geometry_f90.
  //       if parallel/params have not been called yet, the code will abort.
  if (!is_geometry_inited_f90()) {
    init_geometry_f90();
  }

  // We know we need the dyn grid, so build it
  build_dynamics_grid ();

  for (const auto& gn : grid_names) {
    build_physics_grid(gn);
  }

  // Set the ptr to the ref grid
  m_grids["Reference"] = get_grid(reference_grid); 
}

void DynamicsDrivenGridsManager::build_dynamics_grid () {
  if (m_grids.find("Dynamics")==m_grids.end()) {

    // Initialize the dyn grid
    const int nelemd = get_num_owned_elems_f90();
    const int nlev   = get_nlev_f90();
    auto dyn_grid = std::make_shared<SEGrid>("SE Dynamics",nelemd,NP,nlev);

    const int ndofs = nelemd*NP*NP;

    // Create dynamics dofs map
    AbstractGrid::dofs_list_type      dofs("dyn dofs",ndofs);
    AbstractGrid::lid_to_idx_map_type lids_to_elgpgp("dyn lid to elgpgp",ndofs,3);
    AbstractGrid::geo_view_type       lat("lat",ndofs);

    auto h_lids_to_elgpgp = Kokkos::create_mirror_view(lids_to_elgpgp);
    auto h_dofs = Kokkos::create_mirror_view(dofs);

    // Get (ie,igp,jgp,gid) data for each dof
    get_dyn_grid_data_f90 (h_dofs.data(),h_lids_to_elgpgp.data(), h_lat.data(), h_lon.data());

    Kokkos::deep_copy(dofs,h_dofs);
    Kokkos::deep_copy(lids_to_elgpgp,h_lids_to_elgpgp);

    dyn_grid->set_dofs (dofs, lids_to_elgpgp);

    // Set the grid in the map
    for (const auto& gn : m_dyn_grid_aliases) {
      m_grids[gn] = dyn_grid;
    }
  }
}

void DynamicsDrivenGridsManager::
build_physics_grid (const std::string& name) {

  // Codes for the physics grids to build
  constexpr int gll = 0;    // Physics GLL
  constexpr int pg2 = 2;    // Physics PG2
  constexpr int pg3 = 3;    // Physics PG3
  constexpr int gll_t = 10;  // Physics GLL Twin
  constexpr int pg2_t = 12;  // Physics PG2 Twin
  constexpr int pg3_t = 13;  // Physics PG3 Twin

  int pg_type;

  if (name=="Physics" || name=="Physics GLL") {
    pg_type = gll;
  } else if (name=="Phys Grid" || name=="Physics PG2") {
    pg_type = pg2;
  } else if (name=="Physics Twin" || name=="Physics GLL Twin") {
    pg_type = gll_t;
  } else if (name=="Phys Grid Twin" || name=="Physics PG2 Twin") {
    pg_type = pg2_t;
  } else if (name=="Physics PG3") {
    pg_type = pg3;
  } else if (name=="Physics PG3 Twin") {
    pg_type = pg3_t;
  }
  
  if (m_grids.find(name)==m_grids.end()) {

    // Initialize the phys grid
    const int nlev = get_nlev_f90();

    // Create the physics dofs map
    const int num_cols = get_num_local_columns_f90 ();
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
