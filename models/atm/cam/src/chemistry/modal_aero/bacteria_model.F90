!===============================================================================
! Bacteria for Modal Aerosol Model
!===============================================================================
module bacteria_model 
  use shr_kind_mod, only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,   only: masterproc
  use abortutils,   only: endrun
  use camsrfexch,   only: cam_in_t
  use tracer_data,  only: trfld, trfile
  use cam_logfile,  only: iulog

  implicit none
  private

  public :: bacteria_emis
  public :: bacteria_names
  public :: bacteria_nbin
  public :: bacteria_nnum
  public :: bacteria_indices
  public :: bacteria_readnl
  public :: bacteria_init
  public :: bacteria_active

  integer, parameter :: bacteria_nbin = 1
  integer, parameter :: bacteria_nnum = 1

! Tracer fields
  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

#if  ( defined MODAL_AERO_3MODE )
  character(len=6), parameter :: bacteria_names(bacteria_nbin+bacteria_nnum) = (/ 'bac_a3', 'num_a3' /)

  ! Set bacteria diameter per mode (m).
  real(r8),         parameter :: bacteria_diameter(bacteria_nbin) = (/ 4.0e-6 /)

#elif ( defined MODAL_AERO_7MODE )
  character(len=6), parameter :: bacteria_names(bacteria_nbin+bacteria_nnum) = (/ 'bac_a7', 'num_a7' /)

  ! Set bacteria diameter per mode (m).
  real(r8),         parameter :: bacteria_diameter(bacteria_nbin) = (/ 4.0e-6 /)
#endif

  integer  :: bacteria_indices(bacteria_nbin+bacteria_nnum)

  real(r8)          :: bacteria_emis_fact = 1.0_r8     ! overall scaling parameter for bacteria emissions

  logical :: bacteria_active = .true.

  integer, parameter :: & ! number of bacteria data fields
         n_bacteria_data = 1

  character(len=32)   :: specifier(n_bacteria_data) = ''
  character(len=256)  :: filename = ' '
  character(len=256)  :: filelist = ' '
  character(len=256)  :: datapath = ' '
  character(len=32)   :: datatype = 'CYCLICAL'
  integer             :: data_cycle_yr = 0
  logical             :: rmv_file = .false.
  integer             :: fixed_ymd = 0
  integer             :: fixed_tod = 0

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
   character(len=*), parameter :: subname = 'bacteria_data_readnl'

   character(len=32)   :: bacteria_specifier(n_bacteria_data)
   character(len=256)  :: bacteria_filename
   character(len=256)  :: bacteria_filelist
   character(len=256)  :: bacteria_datapath
   character(len=32)   :: bacteria_datatype
   integer             :: bacteria_cycle_yr
   logical             :: bacteria_rmv_file
   integer             :: bacteria_fixed_ymd
   integer             :: bacteria_fixed_tod

   namelist /bacteria_nl/ &
      bacteria_specifier, &
      bacteria_filename,  &
      bacteria_filelist,  &
      bacteria_datapath,  &
      bacteria_datatype,  &
      bacteria_rmv_file,  &
      bacteria_cycle_yr,  &
      bacteria_fixed_ymd, &
      bacteria_fixed_tod
   !-----------------------------------------------------------------------------                             
   ! Initialize namelist variables from local module variables.                                               
   bacteria_specifier= specifier
   bacteria_filename = filename
   bacteria_filelist = filelist
   bacteria_datapath = datapath
   bacteria_datatype = datatype
   bacteria_rmv_file = rmv_file
   bacteria_cycle_yr = data_cycle_yr
   bacteria_fixed_ymd= fixed_ymd
   bacteria_fixed_tod= fixed_tod

   ! Read aerosol namelist                                                                                    
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
      call freeunit(unitn)
      close(unitn)
   endif

!      bacteria_specifier, & ! Names of variables containing aerosol data in the prescribed aerosol datasets.  
!      bacteria_filename,  & ! Filename of dataset for prescribed marine organic matter emissions.             
!      bacteria_filelist,  & ! Filename of file that contains a sequence of filenames for prescribed           
!                         & ! aerosols.  The filenames in this file are relative to the directory specied     
!                         & ! by bacteria_datapath.                                                            
!      bacteria_datapath,  & ! Full pathname of the directory that contains the files specified in bacteria_filelist.
!      bacteria_datatype,      & ! Type of time interpolation for data in bacteria files.                       
!                         & ! Can be set to 'CYCLICAL', 'SERIAL', 'INTERP_MISSING_MONTHS', or 'FIXED'.        
!      bacteria_rmv_file,  & ! Remove the file containing prescribed aerosol deposition fluxes from local disk when no longer needed.
!      bacteria_cycle_yr,  & ! The  cycle year of the prescribed aerosol flux data                             
!                         & ! if bacteria_datatype  is 'CYCLICAL'.
!      bacteria_fixed_ymd, & ! The date at which the prescribed aerosol flux data is fixed                     
!                         & ! if bacteria_datatype is 'FIXED'.                                                 
!      bacteria_fixed_tod    ! The time of day (seconds) corresponding to bacteria_fixed_ymd                    
!                           ! at which the prescribed aerosol flux data is fixed                              
!                           ! if bacteria_datatype is 'FIXED'.                                                 
#ifdef SPMD
   ! Broadcast namelist variables                                                                             
   call mpibcast(bacteria_specifier,len(bacteria_specifier)*n_bacteria_data,   mpichar, 0, mpicom)
   call mpibcast(bacteria_filename, len(bacteria_filename),   mpichar, 0, mpicom)
   call mpibcast(bacteria_filelist, len(bacteria_filelist),   mpichar, 0, mpicom)
   call mpibcast(bacteria_datapath, len(bacteria_datapath),   mpichar, 0, mpicom)
   call mpibcast(bacteria_datatype, len(bacteria_datatype),   mpichar, 0, mpicom)
   call mpibcast(bacteria_rmv_file, 1, mpilog, 0, mpicom)
   call mpibcast(bacteria_cycle_yr, 1, mpiint, 0, mpicom)
   call mpibcast(bacteria_fixed_ymd,1, mpiint, 0, mpicom)
   call mpibcast(bacteria_fixed_tod,1, mpiint, 0, mpicom)
#endif

   ! Update module variables with user settings.                                                              
   specifier     = bacteria_specifier
   filename      = bacteria_filename
   filelist      = bacteria_filelist
   datapath      = bacteria_datapath
   datatype      = bacteria_datatype
   rmv_file      = bacteria_rmv_file
   data_cycle_yr = bacteria_cycle_yr
   fixed_ymd     = bacteria_fixed_ymd
   fixed_tod     = bacteria_fixed_tod

  end subroutine bacteria_readnl

  !=============================================================================

  !===============================================================================
  !===============================================================================
  subroutine bacteria_emis( state, cam_in )
    use mo_constants,  only : bacteria_density
    use physconst,     only : pi
    use physics_types, only: physics_state
  ! args
    ! Arguments:
    type(physics_state),    intent(in )   :: state   ! Physics state variables
    type(cam_in_t), target, intent(inout) :: cam_in  ! import state

  ! local vars
    integer    :: ncol, lchnk
    integer :: i, m, ibac, inum
    real(r8) :: x_mton
    real(r8), pointer :: cflx(:,:), bacteria_flux_in(:)

  ! Get model state
    lchnk = state%lchnk
    ncol = state%ncol

  ! Pointer to state emissions
    cflx   => cam_in%cflx

  ! Pointer to input emissions
    nullify(bacteria_flux_in)

    fldloop: do i=1,n_bacteria_data
      select case (trim(fields(i)%fldnam))
      case( "bac_flx" )
         bacteria_flux_in => fields(i)%data(1:ncol,1,lchnk)
      case default
         if ( masterproc ) then
            write(iulog,*) 'Unknown field name '//fields%fldnam//' in progseasalts field ...'
         endif
      end select
    end do fldloop

    ! set bacteria emissions

    ! rebin and adjust bacteria emissons..
    do m = 1,bacteria_nbin

       ! If bacteria will be emitted into more than one bin, this code
       ! will need to be modified.

       ibac = bacteria_indices(m)
       cflx(:ncol,ibac) = bacteria_emis_fact * bacteria_flux_in(:)

       ! Mass-to-number conversion factor
       x_mton = 6._r8 / (pi * bacteria_density * (bacteria_diameter(m)**3._r8))

       inum = bacteria_indices(m+bacteria_nbin)
       cflx(:ncol,ibac) = cflx(:ncol,ibac)*x_mton

    enddo

  end subroutine bacteria_emis


subroutine bacteria_init()
    !-----------------------------------------------------------------------                                     
    !                                            
    ! Purpose: READ INPUT FILES, CREATE FIELDS, and horizontally interpolate bacteria data                       
    !
    ! Method:                                                                                                    
    !                                                                                                            
    ! Author: S. M. Burrows, adapted from dust_initialize                                                        
    !                                                                                                            
    !-----------------------------------------------------------------------                                  
    use tracer_data,      only : trcdata_init
    use cam_history,      only : addfld, add_default, phys_decomp
    use ppgrid,           only : pcols
    use ppgrid,           only : begchunk, endchunk

    !-----------------------------------------------------------------------                                     
    !   ... local variables                                                                                      
    !-----------------------------------------------------------------------                                  
    integer :: i,j,c
    integer :: ierr,did,vid,nlon,nlat,ntimes,ncols
    integer :: number_flds

    if ( masterproc ) then
       write(iulog,*) 'bacteria organics are prescribed in :'//trim(filename)
    endif

    allocate (file%in_pbuf(n_bacteria_data))
    file%in_pbuf(:) = .false.

    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, data_cycle_yr, fixed_ymd, fixed_tod, datatype)

    number_flds = 0
    if (associated(fields)) number_flds = size( fields )

    if ( number_flds .eq. n_bacteria_data ) then
       if ( masterproc ) then
          write(iulog,"(A21,I3,A20)") 'Successfully read in ',number_flds,' bacteria data fields'
       endif
    else if( number_flds < 1 ) then
       if ( masterproc ) then
          write(iulog,*) 'Failed to read in any bacteria data'
          write(iulog,*) ' '
       endif
       return
    else if ( number_flds .ne. n_bacteria_data ) then
       if ( masterproc ) then
          write(iulog,"(A8,I3,A20)") 'Read in ',number_flds,' bacteria data fields'
          write(iulog,"(A9,I3,A20)") 'Expected ',n_bacteria_data,' bacteria data fields'
          write(iulog,*) ' '
          return
       endif
    end if
    ! Following loop adds fields for output.                                                                  
    !   Note that the units are given in fields(i)%units, avgflag='A' indicates output mean                   
    fldloop:do i = 1,n_bacteria_data

       if ( masterproc ) then
          write(iulog,*) 'adding field '//fields(i)%fldnam//' ...'
       endif

!!$       ! Set names of variable tendencies and declare them as history variables
!!$       !    addfld(fname,                 unite,              numlev, avgflag, long_name, decomp_type, ...)

       if ( trim(fields(i)%fldnam) == "bacteria_flux" ) then
          call addfld(trim(fields(i)%fldnam), 'm^-2 s^-1 ', 1, 'A', 'bacteria input data: '//fields(i)%fldnam, phys_decomp )
          call add_default (fields(i)%fldnam, 1, ' ')
       endif

    enddo fldloop
    if ( masterproc ) then
       write(iulog,*) 'Done initializing bacteria emissions data'
    endif

    bacteria_active = .true.

  end subroutine bacteria_init

end module bacteria_model
