!===============================================================================
! Bacteria for Modal Aerosol Model
!===============================================================================
module bacteria_model 
  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,   only: masterproc
  use abortutils,   only: endrun

  implicit none
  private

  public :: bacteria_names
  public :: bacteria_nbin
  public :: bacteria_nnum
  public :: bacteria_indices
  public :: bacteria_emis
  public :: bacteria_readnl
  public :: bacteria_init
  public :: bacteria_active

  integer, parameter :: bacteria_nbin = 2
  integer, parameter :: bacteria_nnum = 2

#if  ( defined MODAL_AERO_3MODE )
  character(len=6), parameter :: bacteria_names(bacteria_nbin+bacteria_nnum) = (/ 'dst_a1', 'dst_a3', 'num_a1', 'num_a3' /)
  real(r8),         parameter :: bacteria_dmt_grd(bacteria_nbin+1) = (/ 0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8/)
  real(r8),         parameter :: bacteria_emis_sclfctr(bacteria_nbin) = (/ 0.032_r8,0.968_r8 /)
#elif ( defined MODAL_AERO_7MODE )
  character(len=6), parameter :: bacteria_names(bacteria_nbin+bacteria_nnum) = (/ 'dst_a5', 'dst_a7', 'num_a5', 'num_a7' /)
  real(r8),         parameter :: bacteria_dmt_grd(bacteria_nbin+1) = (/ 0.1e-6_r8, 2.0e-6_r8, 10.0e-6_r8/)
  real(r8),         parameter :: bacteria_emis_sclfctr(bacteria_nbin) = (/ 0.13_r8, 0.87_r8 /)
#endif

  integer  :: bacteria_indices(bacteria_nbin+bacteria_nnum)
  real(r8) :: bacteria_dmt_vwr(bacteria_nbin)
  real(r8) :: bacteria_stk_crc(bacteria_nbin)

  real(r8)          :: bacteria_emis_fact = -1.e36_r8        ! tuning parameter for bacteria emissions
  character(len=cl) :: soil_erod_file = 'soil_erod_file' ! full pathname for soil erodibility dataset

  logical :: bacteria_active = .false.

 contains

  !=============================================================================
  ! reads bacteria namelist options
  !=============================================================================
  subroutine bacteria_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'bacteria_readnl'

    namelist /bacteria_nl/ bacteria_emis_fact, soil_erod_file

    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'bacteria_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, bacteria_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(bacteria_emis_fact, 1,                   mpir8,   0, mpicom)
    call mpibcast(soil_erod_file, len(soil_erod_file), mpichar, 0, mpicom)
#endif

  end subroutine bacteria_readnl

  !=============================================================================
  !=============================================================================
  subroutine bacteria_init()
    use soil_erod_mod, only: soil_erod_init
    use constituents,  only: cnst_get_ind
    use bacteria_common,   only: bacteria_set_params

    integer :: n

    do n = 1, bacteria_nbin
       call cnst_get_ind(bacteria_names(n), bacteria_indices(n),abort=.false.)
    end do
    do n = 1, bacteria_nnum
       call cnst_get_ind(bacteria_names(bacteria_nbin+n), bacteria_indices(bacteria_nbin+n),abort=.false.)
    enddo 
    bacteria_active = any(bacteria_indices(:) > 0)
    if (.not.bacteria_active) return
   
    call  soil_erod_init( bacteria_emis_fact, soil_erod_file )

    call bacteria_set_params( bacteria_nbin, bacteria_dmt_grd, bacteria_dmt_vwr, bacteria_stk_crc )

  end subroutine bacteria_init

  !===============================================================================
  !===============================================================================
  subroutine bacteria_emis( ncol, lchnk, bacteria_flux_in, cflx, soil_erod )
    use soil_erod_mod, only : soil_erod_fact
    use soil_erod_mod, only : soil_erodibility
    use mo_constants,  only : bacteria_density
    use physconst,     only : pi

  ! args
    integer,  intent(in)    :: ncol, lchnk
    real(r8), intent(in)    :: bacteria_flux_in(:,:)
    real(r8), intent(inout) :: cflx(:,:)
    real(r8), intent(out)   :: soil_erod(:)

  ! local vars
    integer :: i, m, idst, inum
    real(r8) :: x_mton
    real(r8),parameter :: soil_erod_threshold = 0.1_r8

    ! set bacteria emissions

    col_loop: do i =1,ncol

       soil_erod(i) = soil_erodibility( i, lchnk )

       if( soil_erod(i) .lt. soil_erod_threshold ) soil_erod(i) = 0._r8

       ! rebin and adjust bacteria emissons..
       do m = 1,bacteria_nbin

          idst = bacteria_indices(m)

          cflx(i,idst) = sum( -bacteria_flux_in(i,:) ) &
               * bacteria_emis_sclfctr(m)*soil_erod(i)/soil_erod_fact*1.15_r8

          x_mton = 6._r8 / (pi * bacteria_density * (bacteria_dmt_vwr(m)**3._r8))                

          inum = bacteria_indices(m+bacteria_nbin)

          cflx(i,inum) = cflx(i,idst)*x_mton

       enddo

    end do col_loop

  end subroutine bacteria_emis

end module bacteria_model
