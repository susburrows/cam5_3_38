module restFileMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Reads from or writes to/ the CLM restart file.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use decompMod   , only : bounds_type
  use spmdMod     , only : masterproc
  use clm_varpar  , only : crop_prog
  use abortutils  , only : endrun
  use shr_log_mod , only : errMsg => shr_log_errMsg
  use clm_varctl  
  use ncdio_pio   , only : file_desc_t, ncd_pio_createfile, ncd_pio_openfile, ncd_global, &
                           ncd_pio_closefile, ncd_defdim, ncd_putatt, ncd_enddef, check_dim
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: restFile_read
  public :: restFile_write
  public :: restFile_open
  public :: restFile_close
  public :: restFile_getfile
  public :: restFile_filename        ! Sets restart filename
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: restFile_read_pfile     
  private :: restFile_write_pfile    ! Writes restart pointer file
  private :: restFile_closeRestart   ! Close restart file and write restart pointer file
  private :: restFile_dimset
  private :: restFile_dimcheck
  private :: restFile_enddef
  !
  ! !PRIVATE TYPES: None
  private
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine restFile_write( bounds, file, rdate, noptr )
    !
    ! !DESCRIPTION:
    ! Read/write CLM restart file.
    !
    ! !USES:
    use clm_time_manager , only : timemgr_restart_io, get_nstep
    use subgridRestMod   , only : SubgridRest
    use BiogeophysRestMod, only : BiogeophysRest
    use CNRestMod        , only : CNRest
    use CropRestMod      , only : CropRest
    use accumulMod       , only : accumulRest
    use SLakeRestMod     , only : SLakeRest
    use ch4RestMod       , only : ch4Rest
    use histFileMod      , only : hist_restart_ncd
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds          ! bounds
    character(len=*) , intent(in) :: file            ! output netcdf restart file
    character(len=*) , intent(in), optional :: rdate ! restart file time stamp for name
    logical,           intent(in), optional :: noptr ! if should NOT write to the restart pointer file
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid ! netcdf id
    integer :: i       ! index
    logical :: ptrfile ! write out the restart pointer file
    !-----------------------------------------------------------------------

    if ( present(noptr) )then
       ptrfile = .not. noptr
    else
       ptrfile = .true.
    end if

    ! --------------------------------------------
    ! Open restart file
    ! --------------------------------------------

    call restFile_open( flag='write', file=file, ncid=ncid )

    ! --------------------------------------------
    ! Define dimensions and variables
    ! --------------------------------------------

    call restFile_dimset ( ncid )

    ! Define restart file variables

    call timemgr_restart_io(ncid, flag='define')

    call SubgridRest(bounds, ncid, flag='define' )

    call BiogeophysRest(bounds, ncid, flag='define' )

    if (use_cn) then
       call CNRest(bounds,  ncid, flag='define')
       if ( crop_prog ) call CropRest( bounds, ncid, flag='define' )
    end if
    
    call accumulRest( ncid, flag='define' )

    call SLakeRest( ncid, flag='define' )

    if (use_lch4) then
       call ch4Rest (bounds, ncid, flag='define')
    end if

    if (present(rdate)) then 
       call hist_restart_ncd (bounds, ncid, flag='define', rdate=rdate )
    end if

    call restFile_enddef( ncid )

    ! --------------------------------------------
    ! Write restart file variables
    ! --------------------------------------------
    
    call timemgr_restart_io( ncid, flag='write' )

    call SubgridRest( bounds, ncid, flag='write' )

    call BiogeophysRest( bounds, ncid, flag='write' )

    if (use_cn) then
       call CNRest( bounds, ncid, flag='write' )
       if ( crop_prog ) call CropRest( bounds, ncid, flag='write' )
    end if

    call SLakeRest( ncid, flag='write' )
    if (use_lch4) then
       call ch4Rest ( bounds, ncid, flag='write' )
    end if

    call accumulRest( ncid, flag='write' )
    
    call hist_restart_ncd ( bounds, ncid, flag='write' )

    ! --------------------------------------------
    ! Close restart file and write restart pointer file
    ! --------------------------------------------
    
    call restFile_close( ncid )
    call restFile_closeRestart( file )
    
    ! Write restart pointer file
    
    if ( ptrfile ) call restFile_write_pfile( file )
    
    ! Write out diagnostic info

    if (masterproc) then
       write(iulog,*) 'Successfully wrote out restart data at nstep = ',get_nstep()
       write(iulog,'(72a1)') ("-",i=1,60)
    end if
    
  end subroutine restFile_write

  !-----------------------------------------------------------------------
  subroutine restFile_read( bounds, file )
    !
    ! !DESCRIPTION:
    ! Read a CLM restart file.
    !
    ! !USES:
    use BiogeophysRestMod, only : BiogeophysRest
    use CNRestMod        , only : CNRest
    use CropRestMod      , only : CropRest
    use SLakeRestMod     , only : SLakeRest
    use ch4RestMod       , only : ch4Rest
    use accumulMod       , only : accumulRest
    use histFileMod      , only : hist_restart_ncd
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: file  ! output netcdf restart file
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid ! netcdf id
    integer :: i              ! index
    !-----------------------------------------------------------------------

    ! Open file

    call restFile_open( flag='read', file=file, ncid=ncid )

    ! Read file

    call restFile_dimcheck( ncid )

    call SLakeRest( ncid, flag='read' )

    call BiogeophysRest( bounds, ncid, flag='read' )

    if (use_cn) then
       call CNRest( bounds, ncid, flag='read' )
       if ( crop_prog ) call CropRest( bounds, ncid, flag='read' )
    end if

    if (use_lch4) then
       call ch4Rest( bounds, ncid, flag='read' )
    end if

    call accumulRest( ncid, flag='read' )

    call hist_restart_ncd (bounds, ncid, flag='read')

    ! Close file 

    call restFile_close( ncid )

    ! Write out diagnostic info

    if (masterproc) then
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*) 'Successfully read restart data for restart run'
       write(iulog,*)
    end if

  end subroutine restFile_read

  !-----------------------------------------------------------------------
  subroutine restFile_getfile( file, path )
    !
    ! !DESCRIPTION:
    ! Determine and obtain netcdf restart file
    !
    ! !USES:
    use clm_varctl, only : caseid, nrevsn, nsrest, brnch_retain_casename
    use clm_varctl, only : nsrContinue, nsrBranch
    use fileutils , only : getfil
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(out) :: file  ! name of netcdf restart file
    character(len=*), intent(out) :: path  ! full pathname of netcdf restart file
    !
    ! !LOCAL VARIABLES:
    integer :: status                      ! return status
    integer :: length                      ! temporary          
    character(len=256) :: ftest,ctest      ! temporaries
    !-----------------------------------------------------------------------

    ! Continue run:
    ! Restart file pathname is read restart pointer file 

    if (nsrest==nsrContinue) then
       call restFile_read_pfile( path )
       call getfil( path, file, 0 )
    end if

    ! Branch run: 
    ! Restart file pathname is obtained from namelist "nrevsn"
    ! Check case name consistency (case name must be different for branch run, 
    ! unless namelist specification states otherwise)

    if (nsrest==nsrBranch) then
       length = len_trim(nrevsn)
       if (nrevsn(length-2:length) == '.nc') then
          path = trim(nrevsn) 
       else
          path = trim(nrevsn) // '.nc'
       end if
       call getfil( path, file, 0 )

       ! tcraig, adding xx. and .clm2 makes this more robust
       ctest = 'xx.'//trim(caseid)//'.clm2'
       ftest = 'xx.'//trim(file)
       status = index(trim(ftest),trim(ctest))
       if (status /= 0 .and. .not.(brnch_retain_casename)) then
          if (masterproc) then
             write(iulog,*) 'Must change case name on branch run if ',&
                  'brnch_retain_casename namelist is not set'
             write(iulog,*) 'previous case filename= ',trim(file),&
                  ' current case = ',trim(caseid), &
                  ' ctest = ',trim(ctest), &
                  ' ftest = ',trim(ftest)
          end if
          call endrun(msg=errMsg(__FILE__, __LINE__)) 
       end if
    end if

  end subroutine restFile_getfile

  !-----------------------------------------------------------------------
  subroutine restFile_read_pfile( pnamer )
    !
    ! !DESCRIPTION:
    ! Setup restart file and perform necessary consistency checks
    !
    ! !USES:
    use fileutils , only : opnfil, getavu, relavu
    use clm_varctl, only : rpntfil, rpntdir, inst_suffix
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(out) :: pnamer ! full path of restart file
    !
    ! !LOCAL VARIABLES:
    !EOP
    integer :: i                  ! indices
    integer :: nio                ! restart unit
    integer :: status             ! substring check status
    character(len=256) :: locfn   ! Restart pointer file name
    !-----------------------------------------------------------------------

    ! Obtain the restart file from the restart pointer file. 
    ! For restart runs, the restart pointer file contains the full pathname 
    ! of the restart file. For branch runs, the namelist variable 
    ! [nrevsn] contains the full pathname of the restart file. 
    ! New history files are always created for branch runs.

    if (masterproc) then
       write(iulog,*) 'Reading restart pointer file....'
    endif

    nio = getavu()
    locfn = trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
    call opnfil (locfn, nio, 'f')
    read (nio,'(a256)') pnamer
    call relavu (nio)

    if (masterproc) then
       write(iulog,*) 'Reading restart data.....'
       write(iulog,'(72a1)') ("-",i=1,60)
    end if

  end subroutine restFile_read_pfile

  !-----------------------------------------------------------------------
  subroutine restFile_closeRestart( file )
    !
    ! !DESCRIPTION:
    ! Close restart file and write restart pointer file if
    ! in write mode, otherwise just close restart file if in read mode
    !
    ! !USES:
    use clm_time_manager, only : is_last_step
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*) , intent(in) :: file  ! local output filename
    !
    ! !CALLED FROM:
    ! subroutine restart in this module
    !
    ! !REVISION HISTORY:
    ! Author: Mariana Vertenstein
    !
    !
    ! !LOCAL VARIABLES:
    !EOP
    integer :: i                   !index
    !-----------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Successfully wrote local restart file ',trim(file)
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*)
    end if

  end subroutine restFile_closeRestart

  !-----------------------------------------------------------------------
  subroutine restFile_write_pfile( fnamer )
    !
    ! !DESCRIPTION:
    ! Open restart pointer file. Write names of current netcdf restart file.
    !
    ! !USES:
    use clm_varctl, only : rpntdir, rpntfil, inst_suffix
    use fileutils , only : relavu
    use fileutils , only : getavu, opnfil
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fnamer
    !
    ! !LOCAL VARIABLES:
    integer :: m                    ! index
    integer :: nio                  ! restart pointer file
    character(len=256) :: filename  ! local file name
    !-----------------------------------------------------------------------

    if (masterproc) then
       nio = getavu()
       filename= trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
       call opnfil( filename, nio, 'f' )

       write(nio,'(a)') fnamer
       call relavu( nio )
       write(iulog,*)'Successfully wrote local restart pointer file'
    end if

  end subroutine restFile_write_pfile

  !-----------------------------------------------------------------------
  subroutine restFile_open( flag, file, ncid )

    use clm_time_manager, only : get_nstep

    implicit none
    character(len=*),  intent(in) :: flag ! flag to specify read or write
    character(len=*),  intent(in) :: file ! filename
    type(file_desc_t), intent(out):: ncid ! netcdf id

    integer :: omode                              ! netCDF dummy variable
    character(len= 32) :: subname='restFile_open' ! subroutine name

    if (flag == 'write') then

       ! Create new netCDF file (in define mode) and set fill mode
       ! to "no fill" to optimize performance

       if (masterproc) then	
          write(iulog,*)
          write(iulog,*)'restFile_open: writing restart dataset at ',&
               trim(file), ' at nstep = ',get_nstep()
          write(iulog,*)
       end if
       call ncd_pio_createfile(ncid, trim(file))

    else if (flag == 'read') then

       ! Open netcdf restart file

       if (masterproc) then
          write(iulog,*) 'Reading restart dataset'
       end if
       call ncd_pio_openfile (ncid, trim(file), 0)

    end if

  end subroutine restFile_open

  !-----------------------------------------------------------------------
  character(len=256) function restFile_filename( rdate )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use clm_varctl, only : caseid, inst_suffix
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: rdate   ! input date for restart file name 
    !-----------------------------------------------------------------------

    restFile_filename = "./"//trim(caseid)//".clm2"//trim(inst_suffix)//&
         ".r."//trim(rdate)//".nc"
    if (masterproc) then
       write(iulog,*)'writing restart file ',trim(restFile_filename),' for model date = ',rdate
    end if

  end function restFile_filename

  !------------------------------------------------------------------------
  subroutine restFile_dimset( ncid )
    !
    ! !DESCRIPTION:
    ! Read/Write initial data from/to netCDF instantaneous initial data file
    !
    ! !USES:
    use clm_time_manager, only : get_nstep, get_curr_date
    use spmdMod     , only : mpicom, MPI_LOGICAL
    use clm_varctl  , only : caseid, ctitle, version, username, hostname, fsurdat, &
         conventions, source
    use clm_varpar  , only : numrad, nlevlak, nlevsno, nlevgrnd, nlevurb, nlevcan
    use clm_varpar  , only : cft_lb, cft_ub, maxpatch_glcmec
    use decompMod   , only : get_proc_global
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    !
    ! !LOCAL VARIABLES:
    integer :: dimid               ! netCDF dimension id
    integer :: numg                ! total number of gridcells across all processors
    integer :: numl                ! total number of landunits across all processors
    integer :: numc                ! total number of columns across all processors
    integer :: nump                ! total number of pfts across all processors
    integer :: ier                 ! error status
    integer :: strlen_dimid        ! string dimension id
    character(len=  8) :: curdate  ! current date
    character(len=  8) :: curtime  ! current time
    character(len=256) :: str
    character(len= 32) :: subname='restFile_dimset' ! subroutine name
    !------------------------------------------------------------------------

    call get_proc_global(numg, numl, numc, nump)

    ! Define dimensions

    call ncd_defdim(ncid , 'gridcell', numg           ,  dimid)
    call ncd_defdim(ncid , 'landunit', numl           ,  dimid)
    call ncd_defdim(ncid , 'column'  , numc           ,  dimid)
    call ncd_defdim(ncid , 'pft'     , nump           ,  dimid)

    call ncd_defdim(ncid , 'levgrnd' , nlevgrnd       ,  dimid)
    call ncd_defdim(ncid , 'levurb'  , nlevurb        ,  dimid)
    call ncd_defdim(ncid , 'levlak'  , nlevlak        ,  dimid)
    call ncd_defdim(ncid , 'levsno'  , nlevsno        ,  dimid)
    call ncd_defdim(ncid , 'levsno1' , nlevsno+1      ,  dimid)
    call ncd_defdim(ncid , 'levtot'  , nlevsno+nlevgrnd, dimid)
    call ncd_defdim(ncid , 'numrad'  , numrad         ,  dimid)
    call ncd_defdim(ncid , 'levcan'  , nlevcan        ,  dimid)
    call ncd_defdim(ncid , 'string_length', 64        ,  dimid)
    if (create_glacier_mec_landunit) then
       call ncd_defdim(ncid , 'glc_nec', maxpatch_glcmec, dimid)
    end if

    ! Define global attributes

    call ncd_putatt(ncid, NCD_GLOBAL, 'Conventions', trim(conventions))
    call getdatetime(curdate, curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call ncd_putatt(ncid, NCD_GLOBAL, 'history' , trim(str))
    call ncd_putatt(ncid, NCD_GLOBAL, 'username', trim(username))
    call ncd_putatt(ncid, NCD_GLOBAL, 'host'    , trim(hostname))
    call ncd_putatt(ncid, NCD_GLOBAL, 'version' , trim(version))
    call ncd_putatt(ncid, NCD_GLOBAL, 'source'  , trim(source))
    str = '$Id: restFileMod.F90 41292 2012-10-26 13:51:45Z erik $'
    call ncd_putatt(ncid, NCD_GLOBAL, 'revision_id'    , trim(str))
    call ncd_putatt(ncid, NCD_GLOBAL, 'case_title'     , trim(ctitle))
    call ncd_putatt(ncid, NCD_GLOBAL, 'case_id'        , trim(caseid))
    call ncd_putatt(ncid, NCD_GLOBAL, 'surface_dataset', trim(fsurdat))
    call ncd_putatt(ncid, NCD_GLOBAL, 'title', 'CLM Restart information')
    if (create_glacier_mec_landunit) then
       call ncd_putatt(ncid, ncd_global, 'created_glacier_mec_landunits', 'true')
    else
       call ncd_putatt(ncid, ncd_global, 'created_glacier_mec_landunits', 'false')
    end if

    call ncd_putatt(ncid, ncd_global, 'ipft_not_vegetated'                       ,0 ) 
    call ncd_putatt(ncid, ncd_global, 'ipft_needleleaf_evergreen_temperate_tree' ,1 ) 
    call ncd_putatt(ncid, ncd_global, 'ipft_needleleaf_evergreen_boreal_tree'    ,2 ) 
    call ncd_putatt(ncid, ncd_global, 'ipft_needleleaf_deciduous_boreal_tree'    ,3 ) 
    call ncd_putatt(ncid, ncd_global, 'ipft_broadleaf_evergreen_tropical_tree'   ,4 ) 
    call ncd_putatt(ncid, ncd_global, 'ipft_broadleaf_evergreen_temperate_tree'  ,5 ) 
    call ncd_putatt(ncid, ncd_global, 'ipft_broadleaf_deciduous_tropical_tree'   ,6 ) 
    call ncd_putatt(ncid, ncd_global, 'ipft_broadleaf_deciduous_temperate_tree'  ,7 ) 
    call ncd_putatt(ncid, ncd_global, 'ipft_broadleaf_deciduous_boreal_tree'     ,8 ) 
    call ncd_putatt(ncid, ncd_global, 'ipft_broadleaf_evergreen_shrub'           ,9 ) 
    call ncd_putatt(ncid, ncd_global, 'ipft_broadleaf_deciduous_temperate_shrub' ,10) 
    call ncd_putatt(ncid, ncd_global, 'ipft_broadleaf_deciduous_boreal_shrub'    ,11) 
    call ncd_putatt(ncid, ncd_global, 'ipft_c3_arctic_grass'                     ,12) 
    call ncd_putatt(ncid, ncd_global, 'ipft_c3_non-arctic_grass'                 ,13) 
    call ncd_putatt(ncid, ncd_global, 'ipft_c4_grass'                            ,14) 
    call ncd_putatt(ncid, ncd_global, 'ipft_c3_crop'                             ,15) 
    call ncd_putatt(ncid, ncd_global, 'ipft_c3_irrigated'                        ,16) 
    call ncd_putatt(ncid, ncd_global, 'ipft_corn'                                ,17) 
    call ncd_putatt(ncid, ncd_global, 'ipft_irrigated_corn'                      ,18) 
    call ncd_putatt(ncid, ncd_global, 'ipft_spring_temperate_cereal'             ,19) 
    call ncd_putatt(ncid, ncd_global, 'ipft_irrigated_spring_temperate_cereal'   ,20) 
    call ncd_putatt(ncid, ncd_global, 'ipft_winter_temperate_cereal'             ,21) 
    call ncd_putatt(ncid, ncd_global, 'ipft_irrigated_winter_temperate_cereal'   ,22) 
    call ncd_putatt(ncid, ncd_global, 'ipft_soybean'                             ,23) 
    call ncd_putatt(ncid, ncd_global, 'ipft_irrigated_soybean'                   ,24) 

    call ncd_putatt(ncid, ncd_global, 'cft_lb'                                 , cft_lb)
    call ncd_putatt(ncid, ncd_global, 'cft_ub'                                 , cft_ub)

    call ncd_putatt(ncid, ncd_global, 'icol_vegetated_or_bare_soil'            , 1) 
    call ncd_putatt(ncid, ncd_global, 'icol_crop'                              , 2) 
    call ncd_putatt(ncid, ncd_global, 'icol_crop_noncompete'                   , '2*100+m, m=cft_lb,cft_ub')
    call ncd_putatt(ncid, ncd_global, 'icol_landice'                           , 3) 
    call ncd_putatt(ncid, ncd_global, 'icol_landice_multiple_elevation_classes', '4*100+m, m=1,glcnec')  
    call ncd_putatt(ncid, ncd_global, 'icol_deep_lake'                         , 5) 
    call ncd_putatt(ncid, ncd_global, 'icol_wetland'                           , 6) 
    call ncd_putatt(ncid, ncd_global, 'icol_urban_roof'                        , 71)
    call ncd_putatt(ncid, ncd_global, 'icol_urban_sunwall'                     , 72)
    call ncd_putatt(ncid, ncd_global, 'icol_urban_shadewall'                   , 73)
    call ncd_putatt(ncid, ncd_global, 'icol_urban_impervious_road'             , 74)
    call ncd_putatt(ncid, ncd_global, 'icol_urban_pervious_road'               , 75)

    call ncd_putatt(ncid, ncd_global, 'ilun_vegetated_or_bare_soil'             , 1)
    call ncd_putatt(ncid, ncd_global, 'ilun_crop'                               , 2)
    call ncd_putatt(ncid, ncd_global, 'ilun_landice'                            , 3)
    call ncd_putatt(ncid, ncd_global, 'ilun_landice_multiple_elevation_classes' , 4)
    call ncd_putatt(ncid, ncd_global, 'ilun_deep_lake'                          , 5)
    call ncd_putatt(ncid, ncd_global, 'ilun_wetland'                            , 6)
    call ncd_putatt(ncid, ncd_global, 'ilun_urban_tbd'                          , 7)
    call ncd_putatt(ncid, ncd_global, 'ilun_urban_hd'                           , 8)
    call ncd_putatt(ncid, ncd_global, 'ilun_urban_md'                           , 9)

  end subroutine restFile_dimset

  !-----------------------------------------------------------------------
  subroutine restFile_dimcheck( ncid )
    !
    ! !DESCRIPTION:
    ! Check dimensions of restart file
    !
    ! !USES:
    use decompMod,  only : get_proc_global
    use clm_varpar, only : nlevsno, nlevlak, nlevgrnd, nlevurb
    use clm_varctl, only : single_column, nsrest, nsrStartup
    implicit none
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid
    !
    ! !LOCAL VARIABLES:
    integer :: numg     ! total number of gridcells across all processors
    integer :: numl     ! total number of landunits across all processors
    integer :: numc     ! total number of columns across all processors
    integer :: nump     ! total number of pfts across all processors
    character(len=32) :: subname='restFile_dimcheck' ! subroutine name
    !-----------------------------------------------------------------------

    ! Get relevant sizes

    if ( .not. single_column .or. nsrest /= nsrStartup )then
       call get_proc_global(numg, numl, numc, nump)
       call check_dim(ncid, 'gridcell', numg)
       call check_dim(ncid, 'landunit', numl)
       call check_dim(ncid, 'column'  , numc)
       call check_dim(ncid, 'pft'     , nump)
    end if
    call check_dim(ncid, 'levsno'  , nlevsno)
    call check_dim(ncid, 'levgrnd' , nlevgrnd)
    call check_dim(ncid, 'levurb'  , nlevurb)
    call check_dim(ncid, 'levlak'  , nlevlak) 

  end subroutine restFile_dimcheck

  !-----------------------------------------------------------------------
  subroutine restFile_enddef( ncid )
    !
    ! !DESCRIPTION:
    ! Read a CLM restart file.
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    !-----------------------------------------------------------------------

    call ncd_enddef(ncid)

  end subroutine restFile_enddef

  !-----------------------------------------------------------------------
  subroutine restFile_close( ncid )
    !
    ! !DESCRIPTION:
    ! Read a CLM restart file.
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t), intent(inout) :: ncid
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname='restFile_close' ! subroutine name
    !-----------------------------------------------------------------------

    call ncd_pio_closefile(ncid)

  end subroutine restFile_close

end module restFileMod



