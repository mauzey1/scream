module homme_grid_mod
  use iso_c_binding,     only: c_ptr, c_f_pointer, c_int, c_bool
  use shr_kind_mod,      only: r8 => shr_kind_r8
  use parallel_mod,      only: abortmp
  use kinds,             only: iulog
  use homme_context_mod, only: par, masterproc

  implicit none
  private

  ! Routines that modify state
  public :: init_geometry_f90
  public :: finalize_geometry_f90

  ! Routines to get information
  public :: get_cols_gids_f90, get_cols_indices_f90, get_cols_geo_specs_f90
  public :: get_num_local_columns_f90, get_num_global_columns_f90
  public :: get_num_owned_elems_f90, get_np_f90, get_nlev_f90

  integer :: fv_nphys = 0

  ! Global geo info
  real (r8), allocatable :: pg_lat (:)
  real (r8), allocatable :: pg_lon (:)
  real (r8), allocatable :: pg_area (:)

contains

  subroutine init_geometry_f90 () bind(c)
    use prim_driver_base,     only: prim_init1_geometry, prim_init1_cleanup, MetaVertex, GridEdge
    use prim_cxx_driver_base, only: init_cxx_connectivity
    use dimensions_mod,       only: nelemd
    use homme_context_mod,    only: is_geometry_inited, is_parallel_inited, &
                                    elem, par, dom_mt

    if (.not. is_parallel_inited) then
      call abortmp ("Error! 'homme_init_parallel_f90' must be called *before* init_geometry_f90.\n")
    endif

    if (masterproc)  write(iulog,*) "homme_grid_mod: Initing geometry..."

    call prim_init1_geometry(elem,par,dom_mt)

    call init_cxx_connectivity(nelemd,GridEdge,MetaVertex,par)

    ! Cleanup the tmp stuff used in prim_init1_geometry
    call prim_init1_cleanup()

    is_geometry_inited = .true.
  end subroutine init_geometry_f90

  subroutine finalize_geometry_f90 () bind(c)
    use homme_context_mod,    only: is_geometry_inited, masterproc

    if (is_geometry_inited) then
      if (fv_nphys>0) then
        if (masterproc)  write(iulog,*) "homme_grid_mod: phygrid finalization..."
        call gfr_finish()
      endif
    endif

    is_geometry_inited = .false.
  end subroutine finalize_geometry_f90

  subroutine get_dyn_grid_data_f90 (gids_ptr, elgp_ptr, lat_ptr, lon_ptr) bind(c)
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: gids_ptr, elgp_ptr, lat_ptr, lon_ptr

  end subroutine get_dyn_grid_data_f90

  subroutine get_phys_grid_data_f90 (pg_type, gids_ptr, lat_ptr, lon_ptr, area_ptr) bind(c)
    use parallel_mod, only: abortmp
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: gids_ptr, elgp_ptr, lat_ptr, lon_ptr
    integer (kind=c_int), intent(in) :: pg_type
    !
    ! Local(s)
    !
    integer :: N, ngcols

    ! Possible values for pg_type:
    !   0 : physics grid on GLL nodes
    !   2 : physics grid on FV points in 2x2 subcell (PG2)
    !   3 : physics grid on FV points in 3x3 subcell (PG3)
    !  10 : GLL points, use twin columns
    !  12 : PG2 points, use twin columns
    !  13 : PG3 points, use twin columns
    ! In general:
    !   *0 => GLL nodes, *N,N>0 => FV points
    !   1* => redistribute using twin columns

    N = mod(pg_type,10)

    if (N .gt. 0) then
      if (fv_nphys==0) then
        if (masterproc)  write(iulog,*) "homme_grid_mod: phygrid initialization..."
        fv_nphys = N
        call gfr_init(par, elem, fv_nphys)
      elseif (fv_nphys .ne. N) then
        call abortmp ("Error! FV phys grid was already inited with a different order.")
      endif
    endif

    if (pg_type .ge. 10) then
      call abortmp ("Error! Twin columns phys grid not yet implemented.")
    endif

    if (.not. associated(pg_lat)) then
      ngcols = get_num_global_columns_f90()

      allocate(pg_area(ngcols))
      allocate(pg_lat(ngcols))
      allocate(pg_lon(ngcols))

      call compute_global_area (pg_area)
      call compute_global_coords (pg_lat, pg_lon, fv_nphys)
    endif

    call get_my_phys_data (N, gids_ptr, lat_ptr, lon_ptr, area_ptr)
  end subroutine get_phys_grid_data_f90

  subroutine get_cols_gids_f90 (gids_ptr, phys_grid) bind(c)
    use homme_context_mod, only: is_geometry_inited
    use dimensions_mod,    only: np, nelemd
    use element_mod,       only: index_t
    use kinds,             only: long_kind
    use homme_context_mod, only: elem
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: gids_ptr
    logical (kind=c_bool), intent(in) :: phys_grid
    !
    ! Local(s)
    !
    type (index_t) :: idxP
    integer (kind=long_kind), pointer :: gids(:)
    integer :: i,j,icol,ie,ip,jp,fv_cnt,ncols

    if (.not. is_geometry_inited) then
      call abortmp ("Error! Geometry was not inited yet.\n")
    endif

    if (phys_grid) then
      ncols = get_num_local_columns_f90()
    else
      ncols = nelemd*np*np
    endif

    call c_f_pointer(gids_ptr, gids, [ncols])

    icol = 1
    do ie=1,nelemd
      if (phys_grid .and. fv_nphys>0) then
        fv_cnt = 1
        do j = 1,fv_nphys
          do i = 1,fv_nphys
            gids(icol) = (elem(ie)%GlobalId-1)*fv_nphys*fv_nphys + fv_cnt
            fv_cnt = fv_cnt+1
            icol = icol+1
          end do
        end do
      else
        idxP = elem(ie)%idxP
        if (phys_grid) then
          ncols = idxP%NumUniquePts
        else
          ncols = np*np
        endif
        do i=1,ncols
          ip = idxP%ia(i)
          jp = idxP%ja(i)

          gids(icol) = elem(ie)%gdofP(ip,jp)-1
          icol = icol+1
        enddo
      endif
    enddo
  end subroutine get_cols_gids_f90

  subroutine get_cols_indices_f90 (gids_ptr, elgp_ptr, phys_grid) bind(c)
    use homme_context_mod, only: is_geometry_inited, elem
    use dimensions_mod,    only: np, nelemd
    use element_mod,       only: index_t
    use kinds,             only: long_kind
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: gids_ptr, elgp_ptr
    logical (kind=c_bool), intent(in) :: phys_grid
    !
    ! Local(s)
    !
    type (index_t) :: idxP
    integer (kind=long_kind), pointer :: gids(:)
    integer (kind=c_int), pointer :: elgp(:,:)
    integer :: i,j,icol,ie,ip,jp,fv_cnt,ncols

    if (.not. is_geometry_inited) then
      call abortmp ("Error! Geometry was not inited yet.\n")
    endif

    if (phys_grid) then
      ncols = get_num_local_columns_f90()
    else
      ncols = nelemd*np*np
    endif

    call c_f_pointer(gids_ptr, gids, [ncols])
    call c_f_pointer(elgp_ptr, elgp, [3,ncols])

    if (fv_nphys>0) then
      do ie = 1,nelemd
        fv_cnt = 1
        do j = 1,fv_nphys
          do i = 1,fv_nphys
            gids(icol) = (elem(ie)%GlobalId-1)*fv_nphys*fv_nphys + fv_cnt
            fv_cnt = fv_cnt+1
            icol = icol+1
          end do
        end do
      end do
    else
      icol = 1
      do ie=1,nelemd
        idxP = elem(ie)%idxP
        if (phys_grid) then
          ncols = idxP%NumUniquePts
          do i=1,ncols
            ip = idxP%ia(i)
            jp = idxP%ja(i)
            elgp(1,icol) = ie-1
            elgp(2,icol) = jp-1
            elgp(3,icol) = ip-1

            gids(icol) = elem(ie)%gdofP(ip,jp)-1
            icol = icol + 1
          enddo
        else
          do jp=1,4
            do ip=1,4
              elgp(1,icol) = ie-1
              elgp(2,icol) = jp-1
              elgp(3,icol) = ip-1

              gids(icol) = elem(ie)%gdofP(ip,jp)-1
              icol = icol + 1
            enddo
          enddo
        endif
      enddo
    endif
  end subroutine get_cols_indices_f90

  subroutine get_cols_geo_specs_f90 (coords_ptr, area_ptr, which_grid) bind(c)
    use iso_c_binding,     only: c_double
    use element_mod,       only: index_t
    use homme_context_mod, only: is_geometry_inited, elem
    use dimensions_mod,    only: nelemd
    use dyn_grid,          only: compute_global_area, compute_global_coords
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: coords_ptr, area_ptr
    logical (kind=c_bool), intent(in) :: phys_grid
    !
    ! Local(s)
    !
    real (kind=c_double), pointer :: coords(:,:)
    real (kind=c_double), pointer :: area(:)
    type (index_t) :: idxP
    integer :: icol, ie, ncols, lcols
    ! integer :: i,ip,jp,ie,ncols,icol
    real (kind=c_double), pointer :: lat_g(:)
    real (kind=c_double), pointer :: lon_g(:)
    real (kind=c_double), pointer :: area_g(:)

    if (.not. is_geometry_inited) then
      call abortmp ("Error! Geometry was not inited yet.\n")
    endif

    lcols = get_num_local_columns_f90()

    call c_f_pointer(coords_ptr, coords, [lcols,2])
    call c_f_pointer(area_ptr,   area,   [lcols])

    ncols = get_num_global_columns_f90()

    allocate(lat_g(ncols))
    allocate(lon_g(ncols))
    allocate(area_g(ncols))

    call compute_global_area(area_g)
    call compute_global_coords(clat=lat_g,clon=lon_g)

    deallocate(lat_g)
    deallocate(lon_g)
    deallocate(area_g)
  end subroutine get_cols_geo_specs_f90

  function get_num_local_columns_f90 () result (num_cols) bind(c)
    use homme_context_mod, only: is_geometry_inited, elem
    use dimensions_mod,    only: nelemd, np
    !
    ! Local(s)
    !
    integer (kind=c_int) :: num_cols
    integer :: ie

    if (.not. is_geometry_inited) then
      call abortmp ("Error! Geometry was not inited yet.\n")
    endif

    if (fv_nphys>0) then
      num_cols = nelemd*fv_nphys*fv_nphys
    else
      num_cols = 0
      do ie=1,nelemd
        num_cols = num_cols + elem(ie)%idxP%NumUniquePts
      enddo
    endif
  end function get_num_local_columns_f90

  function get_num_global_columns_f90 () result (num_cols) bind(c)
    use homme_context_mod, only: is_geometry_inited, elem
    use dimensions_mod,    only: nelem, np
    !
    ! Local(s)
    !
    integer (kind=c_int) :: num_cols
    integer :: ie

    if (.not. is_geometry_inited) then
      call abortmp ("Error! Geometry was not inited yet.\n")
    endif

    if (fv_nphys>0) then
      num_cols = nelem*fv_nphys*fv_nphys
    else
      num_cols = nelem*(np-1)*(np-1)+2
    endif
  end function get_num_global_columns_f90

  function get_num_owned_elems_f90 () result(num_elems) bind(c)
    use homme_context_mod, only: is_geometry_inited
    use dimensions_mod,    only: nelemd
    !
    ! Local(s)
    !
    integer (kind=c_int) :: num_elems

    if (.not. is_geometry_inited) then
      call abortmp ("Error! Geometry was not inited yet.\n")
    endif

    num_elems = nelemd
  end function get_num_owned_elems_f90

  function get_np_f90 () result (np_out) bind(c)
    use dimensions_mod, only: np
    !
    ! Local(s)
    !
    integer (kind=c_int) :: np_out

    np_out = np
  end function get_np_f90

  function get_nlev_f90 () result (nlev_out) bind(c)
    use dimensions_mod, only: nlev
    !
    ! Local(s)
    !
    integer (kind=c_int) :: nlev_out

    nlev_out = nlev
  end function get_nlev_f90

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !           Private subroutines below this line             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_my_phys_data (pg_type, gids_ptr, lat_ptr, lon_ptr, area_ptr)
    use parallel_mod, only: abortmp
    !
    ! Input(s)
    !
    type (c_ptr), intent(in) :: gids_ptr, elgp_ptr, lat_ptr, lon_ptr
    integer, intent(in) :: pg_type
    !
    ! Local(s)
    !
    integer :: N, nlcols
    real (kind=c_double), pointer :: lat(:), lon(:), area(:)
    real (kind=c_int),    pointer :: gids(:)

    nlcols = get_num_local_columns_f90()

    call c_f_pointer (lat_ptr,  lat,  [nlcols])
    call c_f_pointer (lon_ptr,  lon,  [nlcols])
    call c_f_pointer (area_ptr, area, [nlcols])
    call c_f_pointer (gids_ptr, gids, [nlcols])

  end subroutine get_my_phys_data

end module homme_grid_mod
