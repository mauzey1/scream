module dyn_grid

  use element_mod,        only: element_t
  use kinds,              only: iulog
  use shr_kind_mod,       only: r8 => shr_kind_r8
  use shr_const_mod,      only: SHR_CONST_PI
  use dimensions_mod,     only: nelem, nelemd, nelemdmax, ne, np, npsq
  use homme_context_mod,  only: par, npes, iam, masterproc

  implicit none
  private

! We need MPI in here, so include it
#include <mpif.h>

  real(r8),        parameter :: rad2deg = 180.0_r8 / SHR_CONST_PI

  public :: compute_global_area
  public :: compute_global_coords

!===================================================================================================
contains
!===================================================================================================

  !
  !=================================================================================================
  !
  subroutine compute_global_area(area_d,fv_nphys)
    use dof_mod,                only: UniqueCoords, UniquePoints
    use gllfvremap_mod,         only: gfr_f_get_area
    !------------------------------Arguments------------------------------------
    real(r8), pointer :: area_d(:)
    integer, intent(in) :: fv_nphys
    !----------------------------Local-Variables--------------------------------
    real(r8), dimension(np,np)  :: areaw
    real(r8), allocatable       :: area_fv(:,:)
    real(r8), allocatable       :: rbuf(:)
    integer,  dimension(npes)   :: displace  ! MPI data displacement for gathering
    integer,  dimension(npes)   :: recvcnts  ! MPI send buffer count for gathering
    integer,  allocatable       :: col_id_global(:)
    integer,  allocatable       :: col_id_local(:)
    real(r8), allocatable       :: area_local(:)
    integer  :: ncol_fv_gbl, ncol_fv_lcl
    integer  :: ie, sb, eb, i, j, ip, fv_cnt, icol
    integer  :: ierr
    integer  :: ibuf
    !---------------------------------------------------------------------------
    if (masterproc) then
      write(iulog,*) 'INFO: Non-scalable action: Computing global area in SE dycore.'
    end if

    if (fv_nphys > 0) then
      
      ! ncol_fv_gbl = fv_nphys*fv_nphys*nelem
      ! ncol_fv_lcl = fv_nphys*fv_nphys*nelemd
      ! allocate(rbuf(ncol_fv_gbl))
      ! allocate(col_id_local(ncol_fv_lcl))
      ! allocate(col_id_global(ncol_fv_gbl))
      ! allocate(area_local(ncol_fv_gbl))

      ! ! Get area for local blocks
      ! icol = 1
      ! do ie = 1,nelemd
      !   do j = 1,fv_nphys
      !     do i = 1,fv_nphys
      !       area_local(icol) = gfr_f_get_area(ie, i, j)
      !       icol = icol+1
      !     end do ! i
      !   end do ! j
      ! end do ! ie

      ! ! gather send buffer count as local cell count
      ! recvcnts(:) = 0
      ! call mpi_allgather(ncol_fv_lcl, 1, mpi_integer, recvcnts, 1, mpi_integer, par%comm, ierr)

      ! ! determine displacement for MPI gather
      ! displace(1) = 0
      ! do ip = 2,npes
      !   displace(ip) = displace(ip-1) + recvcnts(ip-1)
      ! end do
      ! ! Check to make sure we counted correctly
      ! if (masterproc) then
      !   if ( displace(npes) + recvcnts(npes) /= ncol_fv_gbl ) then
      !     call endrun('compute_global_area: bad MPI displace array size')
      !   end if
      ! end if ! masterproc

      ! ! gather element IDs for sorting
      ! col_id_global(:) = -1
      ! icol = 1
      ! do ie = 1,nelemd
      !   fv_cnt = 1
      !   do j = 1,fv_nphys
      !     do i = 1,fv_nphys
      !       col_id_local(icol) = (elem(ie)%GlobalId-1)*fv_nphys*fv_nphys + fv_cnt
      !       fv_cnt = fv_cnt+1
      !       icol = icol+1
      !     end do
      !   end do
      ! end do
      ! call mpi_allgatherv( col_id_local(1:ncol_fv_lcl), recvcnts(iam+1), mpi_integer, col_id_global, &
      !                      recvcnts(:), displace(:), mpi_integer, par%comm, ierr)


      ! call mpi_allgatherv( area_local(:), recvcnts(iam+1), mpi_real8, rbuf, &
      !                      recvcnts(:), displace(:), mpi_real8, par%comm, ierr)

      ! ! sort according to global element ID
      ! do icol = 1,ncol_fv_gbl
      !   area_d( col_id_global(icol) ) = rbuf(icol)
      ! end do

      ! deallocate(rbuf)
      ! deallocate(area_local)
      ! deallocate(col_id_local)
      ! deallocate(col_id_global)

    else ! physics is on GLL grid
    
      ! allocate(rbuf(ngcols_d))
      ! do ie = 1, nelemdmax        
      !   if(ie <= nelemd) then
      !     displace(iam+1) = elem(ie)%idxp%UniquePtOffset-1
      !     recvcnts(iam+1) = elem(ie)%idxP%NumUniquePts
      !     eb = displace(iam+1) + elem(ie)%idxp%NumUniquePts
      !     areaw = 1.0_r8 / elem(ie)%rspheremp(:,:)         
      !     call UniquePoints(elem(ie)%idxP, areaw, area_d(displace(iam+1)+1:eb))
      !   else
      !     displace(iam+1) = 0
      !     recvcnts(iam+1) = 0
      !   end if
      !   ibuf = displace(iam+1)
      !   call mpi_allgather(ibuf, 1, mpi_integer, displace, 1, mpi_integer, par%comm, ierr)
      !   ibuf = recvcnts(iam+1)
      !   call mpi_allgather(ibuf, 1, mpi_integer, recvcnts, 1, mpi_integer, par%comm, ierr)
      !   sb = displace(iam+1) + 1
      !   eb = displace(iam+1) + recvcnts(iam+1)
      !   rbuf(1:recvcnts(iam+1)) = area_d(sb:eb)
      !   call mpi_allgatherv(rbuf, recvcnts(iam+1), mpi_real8, area_d,       &
      !                       recvcnts(:), displace(:), mpi_real8, par%comm, ierr)
      ! end do ! ie
      ! deallocate(rbuf)

    end if ! fv_nphys > 0

  end subroutine compute_global_area
  !
  !=================================================================================================
  !
  subroutine compute_global_coords(clat, clon, lat_out, lon_out, corner_lat_out, corner_lon_out, fv_nphys)
    use dof_mod,                only: UniqueCoords, UniquePoints
    use gllfvremap_mod,         only: gfr_f_get_latlon, gfr_f_get_corner_latlon
    !------------------------------Arguments------------------------------------
    integer, intent(in) :: fv_nphys
    real(r8),           intent(out) :: clat(:)                ! radians
    real(r8),           intent(out) :: clon(:)                ! radians
    real(r8), optional, intent(out) :: lat_out(:)             ! degrees
    real(r8), optional, intent(out) :: lon_out(:)             ! degrees
    real(r8), optional, intent(out) :: corner_lat_out(:,:)    ! degrees
    real(r8), optional, intent(out) :: corner_lon_out(:,:)    ! degrees
    !----------------------------Local-Variables--------------------------------
    real(r8), allocatable     :: rbuf(:)
    integer,  dimension(npes) :: displace  ! MPI data displacement for gathering
    integer,  dimension(npes) :: recvcnts  ! MPI send buffer count for gathering
    integer,  allocatable     :: col_id_global(:)
    integer,  allocatable     :: col_id_local(:)
    real(r8), allocatable     :: lat_rad_local(:)
    real(r8), allocatable     :: lon_rad_local(:)
    real(r8), allocatable     :: corner_lat_rad_local(:,:)
    real(r8), allocatable     :: corner_lon_rad_local(:,:)
    integer  :: ncol_fv_gbl, ncol_fv_lcl
    integer  :: ie, sb, eb, j, i, ip, fv_cnt, icol, c
    integer  :: ierr
    integer  :: ibuf
    !---------------------------------------------------------------------------
    if (masterproc) then
      write(iulog,*) 'INFO: Non-scalable action: Computing global coords in SE dycore.'
    end if

    if (fv_nphys > 0) then

      ! ncol_fv_gbl = fv_nphys*fv_nphys*nelem
      ! ncol_fv_lcl = fv_nphys*fv_nphys*nelemd
      ! allocate(rbuf(ncol_fv_gbl))
      ! allocate(col_id_local(ncol_fv_lcl))
      ! allocate(col_id_global(ncol_fv_gbl))
      ! allocate(lat_rad_local(ncol_fv_lcl))
      ! allocate(lon_rad_local(ncol_fv_lcl))

      ! ! calculate coordinates for local blocks
      ! icol = 1
      ! do ie = 1,nelemd
      !   do j = 1,fv_nphys
      !     do i = 1,fv_nphys
      !       call gfr_f_get_latlon(ie, i, j, lat_rad_local(icol), lon_rad_local(icol))
      !       icol = icol+1
      !     end do ! i
      !   end do ! j
      ! end do ! ie

      ! ! gather send buffer count as local cell count
      ! recvcnts(:) = 0
      ! call mpi_allgather(ncol_fv_lcl, 1, mpi_integer, recvcnts, 1, mpi_integer, par%comm, ierr)

      ! ! determine displacement for MPI gather
      ! displace(1) = 0
      ! do ip = 2,npes
      !   displace(ip) = displace(ip-1) + recvcnts(ip-1)
      ! end do

      ! ! Check to make sure we counted correctly
      ! if (masterproc) then
      !   if ( displace(npes) + recvcnts(npes) /= ncol_fv_gbl ) then
      !     call endrun('compute_global_coords: bad MPI displace array size')
      !   end if
      ! end if

      ! ! gather element IDs for sorting
      ! col_id_global(:) = -1
      ! icol = 1
      ! do ie = 1,nelemd
      !   fv_cnt = 1
      !   do j = 1,fv_nphys
      !     do i = 1,fv_nphys
      !       col_id_local(icol) = (elem(ie)%GlobalId-1)*fv_nphys*fv_nphys + fv_cnt
      !       fv_cnt = fv_cnt+1
      !       icol = icol+1
      !     end do
      !   end do
      ! end do
      ! call mpi_allgatherv( col_id_local(1:ncol_fv_lcl), recvcnts(iam+1), mpi_integer, col_id_global, &
      !                      recvcnts(:), displace(:), mpi_integer, par%comm, ierr)

      ! ! Gather global latitudes
      ! call mpi_allgatherv( lat_rad_local(:), recvcnts(iam+1), mpi_real8, rbuf, &
      !                      recvcnts(:), displace(:), mpi_real8, par%comm, ierr)
      
      ! ! sort latitude according to global element ID
      ! do icol = 1,ncol_fv_gbl
      !   clat( col_id_global(icol) ) = rbuf(icol)
      ! end do

      ! ! Gather global longitudes
      ! call mpi_allgatherv( lon_rad_local(:), recvcnts(iam+1), mpi_real8, rbuf, &
      !                      recvcnts(:), displace(:), mpi_real8, par%comm, ierr)
      
      ! ! sort longitude according to global element ID
      ! do icol = 1,ncol_fv_gbl
      !   clon( col_id_global(icol) ) = rbuf(icol)
      ! end do

      ! ! Create version in degrees if requested
      ! if (present(lat_out)) lat_out(:) = clat(:) * rad2deg
      ! if (present(lon_out)) lon_out(:) = clon(:) * rad2deg

      ! !----------------------------------------------------
      ! ! Get cell corners (only needed for writing scrip file)
      ! !----------------------------------------------------
      ! if ( present(corner_lat_out) .or. present(corner_lon_out) ) then
      !   allocate(corner_lat_rad_local(ncol_fv_lcl,4))
      !   allocate(corner_lon_rad_local(ncol_fv_lcl,4))
      !   icol = 1
      !   do ie = 1,nelemd
      !     do j = 1,fv_nphys
      !       do i = 1,fv_nphys
      !         do c = 1,4
      !           call gfr_f_get_corner_latlon(ie, i, j, c, &
      !                corner_lat_rad_local(icol,c), corner_lon_rad_local(icol,c))
      !         end do
      !         icol = icol+1
      !       end do ! i
      !     end do ! j
      !   end do ! ie
      !   ! Gather corner coordinates (one corner at a time)
      !   do c = 1,4
      !     ! Gather latitude
      !     call mpi_allgatherv( corner_lat_rad_local(:,c), recvcnts(iam+1), mpi_real8, rbuf, &
      !                          recvcnts(:), displace(:), mpi_real8, par%comm, ierr)
      !     ! sort latitude 
      !     do icol = 1,ncol_fv_gbl
      !       corner_lat_out( col_id_global(icol), c) = rbuf(icol) * rad2deg
      !     end do
      !     ! Gather longitude
      !     call mpi_allgatherv( corner_lon_rad_local(:,c), recvcnts(iam+1), mpi_real8, rbuf, &
      !                          recvcnts(:), displace(:), mpi_real8, par%comm, ierr)
      !     ! sort longitude 
      !     do icol = 1,ncol_fv_gbl
      !       corner_lon_out( col_id_global(icol), c) = rbuf(icol) * rad2deg
      !     end do
      !   end do ! c
      !   ! Deallocate stuff for corners
      !   deallocate(corner_lat_rad_local)
      !   deallocate(corner_lon_rad_local)
      ! end if ! present(corner_lat_out) .or. present(corner_lon_out)
      ! !----------------------------------------------------
      ! !----------------------------------------------------

      ! ! Deallocate stuff
      ! deallocate(rbuf)
      ! deallocate(lat_rad_local)
      ! deallocate(lon_rad_local)
      ! deallocate(col_id_local)
      ! deallocate(col_id_global)

    else ! physics is on GLL grid

      ! allocate(rbuf(ngcols_d))

      ! clat(:) = -iam
      ! clon(:) = -iam
      ! if (present(lon_out)) then
      !   lon_out(:) = -iam
      ! end if
      ! if (present(lat_out)) then
      !   lat_out(:) = -iam
      ! end if

      ! do ie = 1, nelemdmax
      !   if(ie <= nelemd) then
      !     displace(iam+1) = elem(ie)%idxp%UniquePtOffset-1
      !     eb = displace(iam+1) + elem(ie)%idxp%NumUniquePts
      !     recvcnts(iam+1) = elem(ie)%idxP%NumUniquePts
      !     call UniqueCoords(elem(ie)%idxP, elem(ie)%spherep,                    &
      !          clat(displace(iam+1)+1:eb),                                       &
      !          clon(displace(iam+1)+1:eb))
      !     if (present(lat_out)) then
      !       lat_out(displace(iam+1)+1:eb) = clat(displace(iam+1)+1:eb) * rad2deg
      !     end if
      !     if (present(lon_out)) then
      !       lon_out(displace(iam+1)+1:eb) = clon(displace(iam+1)+1:eb) * rad2deg
      !     end if
      !   else
      !     displace(iam+1) = 0
      !     recvcnts(iam+1) = 0
      !   end if
      !   ibuf = displace(iam+1)
      !   call mpi_allgather(ibuf, 1, mpi_integer, displace, &
      !        1, mpi_integer, par%comm, ierr)

      !   ibuf = recvcnts(iam+1)
      !   call mpi_allgather(ibuf, 1, mpi_integer, recvcnts, &
      !        1, mpi_integer, par%comm, ierr)

      !   sb = displace(iam+1) + 1
      !   eb = displace(iam+1) + recvcnts(iam+1)

      !   rbuf(1:recvcnts(iam+1)) = clat(sb:eb)
      !   call mpi_allgatherv(rbuf, recvcnts(iam+1), mpi_real8, clat,            &
      !          recvcnts(:), displace(:), mpi_real8, par%comm, ierr)
      !   if (present(lat_out)) then
      !     rbuf(1:recvcnts(iam+1)) = lat_out(sb:eb) 
      !     call mpi_allgatherv(rbuf, recvcnts(iam+1), mpi_real8, lat_out,       &
      !          recvcnts(:), displace(:), mpi_real8, par%comm, ierr)
      !   end if

      !   rbuf(1:recvcnts(iam+1)) = clon(sb:eb)
      !   call mpi_allgatherv(rbuf, recvcnts(iam+1), mpi_real8, clon,            &
      !          recvcnts(:), displace(:), mpi_real8, par%comm, ierr)
      !   if (present(lon_out)) then
      !     rbuf(1:recvcnts(iam+1)) = lon_out(sb:eb) 
      !     call mpi_allgatherv(rbuf, recvcnts(iam+1), mpi_real8, lon_out,       &
      !          recvcnts(:), displace(:), mpi_real8, par%comm, ierr)
      !   end if
      ! end do
      
      ! deallocate(rbuf)

    end if ! fv_nphys > 0

  end subroutine compute_global_coords
  !
  !=================================================================================================
  !
end module dyn_grid
