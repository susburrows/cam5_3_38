module reweightMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handles modifications and error-checks related to changing subgrid weights
  !
  ! ----- Requirements for subgrid weights that are enforced here -----
  !
  ! (These requirements are checked in checkWeights/weightsOkay)
  !
  ! Note: in the following, 'active' refers to a pft, column, landunit or grid cell over
  ! which computations are performed, and 'inactive' refers to a pft, column or landunit
  ! where computations are NOT performed (grid cells are always active).
  ! 
  ! (1) For all columns, landunits and grid cells, the sum of all subgrid weights of its
  !     children (or grandchildren, etc.) is equal to 1. For example:
  !     - For all columns, the sum of all pft weights on the column equals 1
  !     - For all landunits, the sum of all col weights on the landunit equals 1
  !     - For all grid cells, the sum of all pft weights on the grid cell equals 1
  !     - etc.
  ! 
  ! (2) For all ACTIVE columns, landunits and grid cells, the sum of all subgrid weights of
  !     its ACTIVE children (or grandchildren, etc.) is equal to 1. For example:
  !     - For all active columns, the sum of all pft weights on the column equals 1 when
  !       just considering active pfts
  !     - For all active landunits, the sum of all col weights on the landunit equals 1 when
  !       just considering active cols
  !     - For ALL grid cells, the sum of all pft weights on the grid cell equals 1 when
  !       just considering active pfts -- note that all grid cells are considered active!
  !     - etc.
  !
  ! (3) For all INACTIVE columns, landunits and grid cells, the sum of all subgrid weights of
  !     its ACTIVE children, grandchildren, etc. are equal to either 0 or 1. For example:
  !     - For all inactive columns, the sum of all pft weights on the column equals either 0
  !       or 1 when just considering active pfts
  !     - For all inactive landunits, the sum of all col weights on the landunit equals
  !       either 0 or 1 when just considering active cols
  !     - etc.
  !
  ! Another way of stating (2) and (3) is that the sum of weights of all ACTIVE pfts, cols
  ! or landunits on their parent/grandparent/etc. is always equal to either 0 or 1 -- and
  ! must be equal to 1 if this parent/grandparent, etc. is itself active.
  !
  ! Note that, together, conditions (1) and (2) imply that any pft, col or landunit whose
  ! weight on the grid cell is non-zero must be active. In addition, these conditions imply
  ! that any pft whose weight on the column is non-zero must be active if the column is
  ! active (and similarly for any pft on an active landunit, and any col on an active
  ! landunit).
  !
  !
  ! ----- Implications of these requirements for computing subgrid averages -----
  !
  ! The preferred way to average from, say, pft to col is:
  !    colval(c) = 0
  !    do p = pfti(c), pftf(c)
  !       if (active(p)) colval(c) = colval(c) + pftval(p) * wtcol(p)
  ! (where wtcol(p) is the weight of the pft on the column)
  ! If column c is active, then the above conditions guarantee that the pwtcol values
  ! included in the above sum will sum to 1. If column c is inactive, then the above
  ! conditions guarantee that the pwtcol values included in the above sum will sum to
  ! either 1 or 0; if they sum to 0, then colval(c) will remain 0.
  !
  ! Another acceptable method is the following; this method accommodates some unknown
  ! fraction of pftval's being set to spval, and leaves colval set at spval if there are no
  ! valid pft values:
  !    colval(c) = spval
  !    sumwt(c) = 0
  !    do p = pfti(c), pftf(c)
  !       if (active(p) .and. wtcol(p) /= 0) then
  !          if (pftval(p) /= spval) then
  !             if (sumwt(c) == 0) colval(c) = 0
  !             colval(c) = colval(c) + pftval(p) * wtcol(p)
  !             sumwt(c) = sumwt(c) + wtcol(p)
  !          end if
  !       end if
  !    end do
  !    if (sumwt(c) /= 0) then
  !       colval(c) = colval(c) / sumwt(c)
  !    end if
  ! Note that here we check the condition (active(p) .and. wtcol(p) /= 0). We need to
  ! include a check for wtcol(p) /= 0 because we don't want to set colval(c) = 0 for zero-
  ! weight pfts in this line:
  !             if (sumwt(c) == 0) colval(c) = 0
  ! And we include a check for active(p) because we don't want to assume that pftval(p) has
  ! been set to spval for inactive pfts -- we want to allow for the possibility that
  ! pftval(p) will be NaN for inactive pfts.
  !
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog, all_active
  use decompMod      , only : bounds_type
  use shr_assert_mod , only : shr_assert
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  !
  ! PUBLIC TYPES:
  implicit none
  save

  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: compute_higher_order_weights  ! given p2c, c2l and l2g weights, compute other weights
  public :: reweightWrapup                ! do modifications and error-checks after modifying subgrid weights
  !
  ! !REVISION HISTORY:
  ! Created by Bill Sacks
  !
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: setActive                       ! set 'active' flags at pft, column & landunit level
  private :: is_active_l                     ! determine whether the given landunit is active
  private :: is_active_c                     ! determine whether the given column is active
  private :: is_active_p                     ! determine whether the given pft is active
  private :: is_gcell_of_landunit_all_istice ! determine whether a grid cell is 100% covered by the istice landunit
  private :: checkWeights                    ! check subgrid weights
  private :: weightsOkay                     ! determine if sum of weights satisfies requirements laid out above
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine compute_higher_order_weights(bounds)
    !
    ! !DESCRIPTION:
    ! Assuming pft%wtcol, col%wtlunit and lun%wtgcell have already been computed, compute
    ! the "higher-order" weights: pft%wtlunit, pft%wtgcell and col%wtgcell, for all p and c
    !
    ! !USES:
    use clmtype
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! clump bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p, c, l      ! indices for pft, col & landunit
    !------------------------------------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       col%wtgcell(c) = col%wtlunit(c) * lun%wtgcell(l)
    end do

    do p = bounds%begp, bounds%endp
       c = pft%column(p)
       pft%wtlunit(p) = pft%wtcol(p) * col%wtlunit(c)
       pft%wtgcell(p) = pft%wtcol(p) * col%wtgcell(c)
    end do
  end subroutine compute_higher_order_weights

  !-----------------------------------------------------------------------
  subroutine reweightWrapup(bounds)
    !
    ! !DESCRIPTION:
    ! Do additional modifications and error-checks that should be done after modifying subgrid
    ! weights
    !
    ! This should be called whenever any weights change (e.g., pft weights on the column,
    ! landunit weights on the grid cell, etc.).
    !
    ! !USES:
    use filterMod, only : setFilters
    use decompMod, only : BOUNDS_LEVEL_CLUMP
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! clump bounds
    !------------------------------------------------------------------------

    call shr_assert(bounds%level == BOUNDS_LEVEL_CLUMP, errMsg(__FILE__, __LINE__))

    call setActive(bounds)
    call checkWeights(bounds, active_only=.false.)
    call checkWeights(bounds, active_only=.true.)
    call setFilters(bounds)

  end subroutine reweightWrapup

  !-----------------------------------------------------------------------
  subroutine setActive(bounds)
    !
    ! !DESCRIPTION:
    ! Set 'active' flags at the pft, column and landunit level
    ! (note that grid cells are always active)
    !
    ! This should be called whenever any weights change (e.g., pft weights on the column,
    ! landunit weights on the grid cell, etc.).
    !
    ! Ensures that we don't have any active pft on an inactive column, or an active column on an
    ! inactive landunit (since these conditions could lead to garbage data)
    !
    ! !USES:
    use clmtype
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: l,c,p       ! loop counters
    logical :: error_found ! true if we find an error

    character(len=*), parameter :: subname = 'setActive'
    !------------------------------------------------------------------------

    error_found = .false.

    do l = bounds%begl,bounds%endl
       lun%active(l) = is_active_l(l)
    end do

    do c = bounds%begc,bounds%endc
       l = col%landunit(c)
       col%active(c) = is_active_c(c)
       if (col%active(c) .and. .not. lun%active(l)) then
          write(iulog,*) trim(subname),' ERROR: active column found on inactive landunit', &
                         'at c = ', c, ', l = ', l
          error_found = .true. 
       end if
       if (error_found) then
          call endrun(decomp_index=c, clmlevel=namec, msg=errMsg(__FILE__, __LINE__))
       end if
    end do

    do p = bounds%begp,bounds%endp
       c = pft%column(p)
       pft%active(p) = is_active_p(p)
       if (pft%active(p) .and. .not. col%active(c)) then
          write(iulog,*) trim(subname),' ERROR: active pft found on inactive column', &
                         'at p = ', p, ', c = ', c
          error_found = .true. 
       end if
       if (error_found) then
          call endrun(decomp_index=p, clmlevel=namep, msg=errMsg(__FILE__, __LINE__))
       end if
    end do

  end subroutine setActive

  !-----------------------------------------------------------------------
  logical function is_active_l(l)
    !
    ! !DESCRIPTION:
    ! Determine whether the given landunit is active
    !
    ! !USES:
    use clmtype
    use clm_varcon, only : istsoil, istice_mec
    use domainMod , only : ldomain
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: l   ! landunit index
    !
    ! !LOCAL VARIABLES:
    integer :: g  ! grid cell index
    !------------------------------------------------------------------------

    if (all_active) then
       is_active_l = .true.

    else
       g =lun%gridcell(l)

       is_active_l = .false.

       ! ------------------------------------------------------------------------
       ! General conditions under which is_active_l NEEDS to be true in order to satisfy
       ! the requirements laid out at the top of this module:
       ! ------------------------------------------------------------------------
       if (lun%wtgcell(l) > 0) is_active_l = .true.

       ! ------------------------------------------------------------------------
       ! Conditions under which is_active_p is set to true because we want extra virtual landunits:
       ! ------------------------------------------------------------------------

       ! Always run over ice_mec landunits within the glcmask, because this is where glc
       ! might need input from virtual (0-weight) landunits
       if (lun%itype(l) == istice_mec .and. ldomain%glcmask(g) == 1) is_active_l = .true.

       ! In general, include a virtual natural vegetation landunit. This aids
       ! initialization of a new landunit; and for runs that are coupled to CISM, this
       ! provides bare land SMB forcing even if there is no vegetated area.
       !
       ! However, we do NOT include a virtual vegetated column in grid cells that are 100%
       ! standard (non-mec) glacier. This is for performance reasons: for FV 0.9x1.25,
       ! excluding these virtual vegetated columns (mostly over Antarctica) leads to a ~
       ! 6% performance improvement (the performance improvement is much less for ne30,
       ! though). In such grid cells, we do not need the forcing to CISM (because if we
       ! needed forcing to CISM, we'd be using an istice_mec point rather than plain
       ! istice). Furthermore, standard glacier landunits cannot retreat (only istice_mec
       ! points can retreat, due to coupling with CISM), so we don't need to worry about
       ! the glacier retreating in this grid cell, exposing new natural veg area. The
       ! only thing that could happen is the growth of some special landunit - e.g., crop
       ! - in this grid cell, due to dynamic landunits. We'll live with the fact that
       ! initialization of the new crop landunit will be initialized in an un-ideal way
       ! in this rare situation.
       if (lun%itype(l) == istsoil .and. .not. is_gcell_of_landunit_all_istice(l)) then
          is_active_l = .true.
       end if

    end if

  end function is_active_l

  !-----------------------------------------------------------------------
  logical function is_active_c(c)
    !
    ! !DESCRIPTION:
    ! Determine whether the given column is active
    !
    ! !USES:
    use clmtype
    use clm_varcon, only : istice_mec, isturb_MIN, isturb_MAX
    use domainMod , only : ldomain
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: c   ! column index
    !
    ! !LOCAL VARIABLES:
    integer :: l  ! landunit index
    integer :: g  ! grid cell index
    !------------------------------------------------------------------------

    if (all_active) then
       is_active_c = .true.

    else
       l =col%landunit(c)
       g =col%gridcell(c)

       is_active_c = .false.

       ! ------------------------------------------------------------------------
       ! General conditions under which is_active_c NEEDS to be true in order to satisfy
       ! the requirements laid out at the top of this module:
       ! ------------------------------------------------------------------------
       if (lun%active(l) .and. col%wtlunit(c) > 0._r8) is_active_c = .true.

       ! ------------------------------------------------------------------------
       ! Conditions under which is_active_c is set to true because we want extra virtual columns:
       ! ------------------------------------------------------------------------

       ! always run over all ice_mec columns within the glcmask, because this is where glc
       ! might need input from virtual (0-weight) columns
       if (lun%itype(l) == istice_mec .and. ldomain%glcmask(g) == 1) is_active_c = .true.

       ! We don't really need to run over 0-weight urban columns. But because of some
       ! messiness in the urban code (many loops are over the landunit filter, then drill
       ! down to columns - so we would need to add 'col%active(c)' conditionals in many
       ! places) it keeps the code cleaner to run over 0-weight urban columns. This generally
       ! shouldn't add much computation time, since in most places, all urban columns are
       ! non-zero weight if the landunit is non-zero weight.
       if (lun%active(l) .and. (lun%itype(l) >= isturb_MIN .and. lun%itype(l) <= isturb_MAX)) then
          is_active_c = .true.
       end if
    end if

  end function is_active_c

  !-----------------------------------------------------------------------
  logical function is_active_p(p)
    !
    ! !DESCRIPTION:
    ! Determine whether the given pft is active
    !
    ! !USES:
    use clmtype
    use clm_varcon, only : istice_mec
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p   ! pft index
    !
    ! !LOCAL VARIABLES:
    integer :: c  ! column index
    !------------------------------------------------------------------------

    if (all_active) then
       is_active_p = .true.

    else
       c =pft%column(p)
    
       is_active_p = .false.

       ! ------------------------------------------------------------------------
       ! General conditions under which is_active_p NEEDS to be true in order to satisfy
       ! the requirements laid out at the top of this module:
       ! ------------------------------------------------------------------------
       if (col%active(c) .and. pft%wtcol(p) > 0._r8) is_active_p = .true.

    end if

  end function is_active_p

  !-----------------------------------------------------------------------
  function is_gcell_of_landunit_all_istice(l) result(all_istice)
    !
    ! !DESCRIPTION:
    ! For a given landunit, determine whether its grid cell is 100% covered by the istice
    ! landunit (not to be confused with istice_mec)
    !
    ! !USES:
    use clmtype    , only : lun, grc
    use clm_varcon , only : istice, ispval
    !
    ! !ARGUMENTS:
    implicit none
    logical :: all_istice        ! function result
    integer, intent(in) :: l     ! landunit index (i.e., index into landunit-level arrays)
    !
    ! !LOCAL VARIABLES:
    integer :: g        ! grid cell index
    integer :: l_istice ! index of the istice landunit in this grid cell

    real(r8), parameter :: tolerance = 1.e-13_r8  ! tolerance for checking whether landunit's weight is 1
    character(len=*), parameter :: subname = 'is_gcell_of_landunit_all_istice'
    !------------------------------------------------------------------------------

    g = lun%gridcell(l)
    l_istice = grc%landunit_indices(istice, g)
    if (l_istice == ispval) then
       ! There is no istice landunit on this grid cell
       all_istice = .false.
    else if (lun%wtgcell(l_istice) >= (1._r8 - tolerance)) then
       all_istice = .true.
    else
       all_istice = .false.
    end if

  end function is_gcell_of_landunit_all_istice

  !------------------------------------------------------------------------------
  subroutine checkWeights (bounds, active_only)
    !
    ! !DESCRIPTION:
    ! Check subgrid weights.
    !
    ! This routine operates in two different modes, depending on the value of active_only. If
    ! active_only is true, then we check the sum of weights of the ACTIVE children,
    ! grandchildren, etc. of a given point. If active_only is false, then we check the sum of
    ! weights of ALL children, grandchildren, etc. of a given point. 
    !
    ! Normally this routine will be called twice: once with active_only=false, and once with
    ! active_only=true.
    !
    ! !USES
    use clmtype
    !
    ! !ARGUMENTS
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    logical, intent(in) :: active_only ! true => check sum of weights just of ACTIVE children, grandchildren, etc.
    !
    ! !LOCAL VARIABLES:
    integer :: g,l,c,p     ! loop counters
    real(r8), allocatable :: sumwtcol(:), sumwtlunit(:), sumwtgcell(:)
    logical :: error_found                ! true if we find an error
    character(len=*), parameter :: subname = 'checkWeights'
    !------------------------------------------------------------------------------

    allocate(sumwtcol(bounds%begc:bounds%endc))
    allocate(sumwtlunit(bounds%begl:bounds%endl))
    allocate(sumwtgcell(bounds%begg:bounds%endg))

    error_found = .false.

    ! Check PFT-level weights
    sumwtcol(bounds%begc : bounds%endc) = 0._r8
    sumwtlunit(bounds%begl : bounds%endl) = 0._r8
    sumwtgcell(bounds%begg : bounds%endg) = 0._r8

    do p = bounds%begp,bounds%endp
       c = pft%column(p)
       l = pft%landunit(p)
       g = pft%gridcell(p)

       if ((active_only .and. pft%active(p)) .or. .not. active_only) then 
          sumwtcol(c) = sumwtcol(c) + pft%wtcol(p)
          sumwtlunit(l) = sumwtlunit(l) + pft%wtlunit(p)
          sumwtgcell(g) = sumwtgcell(g) + pft%wtgcell(p)
       end if
    end do

    do c = bounds%begc,bounds%endc
       if (.not. weightsOkay(sumwtcol(c), active_only, col%active(c))) then
          write(iulog,*) trim(subname),' ERROR: at c = ',c,'total PFT weight is ',sumwtcol(c), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do

    do l = bounds%begl,bounds%endl
       if (.not. weightsOkay(sumwtlunit(l), active_only, lun%active(l))) then
          write(iulog,*) trim(subname),' ERROR: at l = ',l,'total PFT weight is ',sumwtlunit(l), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do

    do g = bounds%begg,bounds%endg
       if (.not. weightsOkay(sumwtgcell(g), active_only, i_am_active=.true.)) then
          write(iulog,*) trim(subname),' ERROR: at g = ',g,'total PFT weight is ',sumwtgcell(g), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do

    ! Check col-level weights
    sumwtlunit(bounds%begl : bounds%endl) = 0._r8
    sumwtgcell(bounds%begg : bounds%endg) = 0._r8

    do c = bounds%begc,bounds%endc
       l = col%landunit(c)
       g = col%gridcell(c)

       if ((active_only .and. col%active(c)) .or. .not. active_only) then
          sumwtlunit(l) = sumwtlunit(l) + col%wtlunit(c)
          sumwtgcell(g) = sumwtgcell(g) + col%wtgcell(c)
       end if
    end do

    do l = bounds%begl,bounds%endl
       if (.not. weightsOkay(sumwtlunit(l), active_only, lun%active(l))) then
          write(iulog,*) trim(subname),' ERROR: at l = ',l,'total col weight is ',sumwtlunit(l), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do
    
    do g = bounds%begg,bounds%endg
       if (.not. weightsOkay(sumwtgcell(g), active_only, i_am_active=.true.)) then
          write(iulog,*) trim(subname),' ERROR: at g = ',g,'total col weight is ',sumwtgcell(g), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do

    ! Check landunit-level weights
    sumwtgcell(bounds%begg : bounds%endg) = 0._r8

    do l = bounds%begl,bounds%endl
       g = lun%gridcell(l)
       if ((active_only .and. lun%active(l)) .or. .not. active_only) then
          sumwtgcell(g) = sumwtgcell(g) + lun%wtgcell(l)
       end if
    end do

    do g = bounds%begg,bounds%endg
       if (.not. weightsOkay(sumwtgcell(g), active_only, i_am_active=.true.)) then
          write(iulog,*) trim(subname),' ERROR: at g = ',g,'total lunit weight is ',sumwtgcell(g), &
                         'active_only = ', active_only
          error_found = .true.
       end if
    end do

    deallocate(sumwtcol, sumwtlunit, sumwtgcell)

    if (error_found) then
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! Success

  end subroutine checkWeights

  !-----------------------------------------------------------------------
  logical function weightsOkay(sumwts, active_weights_only, i_am_active)
    !
    ! !DESCRIPTION:
    ! Determine if sumwts (the sum of weights of children, grandchildren or
    ! great-grandchildren of a column, landunit or grid cell) satisfies the requirements laid
    ! out above.
    !
    ! The way this is determined depends on the values of two other variables:
    ! - active_weights_only: does sumwts just include weights of active children,
    !   grandchildren or great-grandchilden? (alternative is that it includes weights of ALL
    !   children, grandchildren or great-grandchildren)
    ! - i_am_active: true if the column, landunit or grid cell of interest is active
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: sumwts              ! sum of weights of children, grandchildren or great-grandchildren
    logical , intent(in) :: active_weights_only ! true if sumwts just includes active children, etc.
    logical , intent(in) :: i_am_active         ! true if the current point is active
    !
    ! !LOCAL VARIABLES:
    logical :: weights_equal_1
    real(r8), parameter :: tolerance = 1.e-12_r8  ! tolerance for checking whether weights sum to 1
    !------------------------------------------------------------------------

    weights_equal_1 = (abs(sumwts - 1._r8) <= tolerance)

    if (active_weights_only) then
       if (i_am_active) then        ! condition (2) above
          weightsOkay = weights_equal_1
       else                         ! condition (3) above
          weightsOkay = (sumwts == 0._r8 .or. weights_equal_1)
       end if
    else                            ! condition (1) above
       ! (note that i_am_active is irrelevant in this case)
       weightsOkay = weights_equal_1
    end if

  end function weightsOkay

end module reweightMod
