!==========================================================================================!
!==========================================================================================!
!     Module with utilities to calculate the horizontal shading factor due to taller       !
! patches, and split patches based on the probability of being shaded by a taller          !
! neighbouring patch.                                                                      !
!------------------------------------------------------------------------------------------!
module hrzshade_utils
   implicit none

   !=======================================================================================!
   !=======================================================================================!



   contains


   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine calculates some derived quantities needed to run the horizontal   !
   ! shading routines.                                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine init_cci_variables()
      use canopy_radiation_coms, only : cci_gapmin     & ! intent(in)
                                      , cci_pixres     & ! intent(in)
                                      , cci_gapsize    & ! intent(in)
                                      , cci_nretn      & ! intent(in)
                                      , cci_gaparea    & ! intent(out)
                                      , rls_nxy        & ! intent(out)
                                      , rls_npixel     & ! intent(out)
                                      , rls_ngap       & ! intent(out)
                                      , rls_length     & ! intent(out)
                                      , rls_area       & ! intent(out)
                                      , rls_x          & ! intent(out)
                                      , rls_y          & ! intent(out)
                                      , rls_ztch       & ! intent(out)
                                      , rls_cci        & ! intent(out)
                                      , rls_fbeam      & ! intent(out)
                                      , rls_igp        & ! intent(out)
                                      , gap_npixel     & ! intent(out)
                                      , gap_x0         & ! intent(out)
                                      , gap_y0         & ! intent(out)
                                      , gap_ipa        & ! intent(out)
                                      , gap_fbeam      & ! intent(out)
                                      , gap_idx        & ! intent(out)
                                      , gap_mask       ! ! intent(out)
      use disturb_coms         , only : min_patch_area ! ! intent(in)
      !----- Local variables. -------------------------------------------------------------!
      integer  :: gap_nxy ! Total number of gaps in each direction.
      integer  :: i       ! Counter
      integer  :: j       ! Counter
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Make sure that cci_gapsize is a perfect multiple of cci_pixres.  In case it   !
      ! is not, adjust cci_pixres to so it becomes a perfect multiple.                     !
      !------------------------------------------------------------------------------------!
      if (cci_gapsize <= 0. .or. cci_pixres <= 0.) then
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
         write(unit=*,fmt='(a)') ' Invalid settings! Gapsize and pixres must be positive.'
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
         write(unit=*,fmt='(a,1x,f8.3)') ' GAPSIZE = ',cci_gapsize
         write(unit=*,fmt='(a,1x,f8.3)') ' PIXRES  = ',cci_pixres
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
         call fatal_error('Invalid gapsize/pixres settings','init_cci_variables'           &
                         ,'hrzshade_utils.f90')
      elseif (mod(cci_gapsize,cci_pixres) /= 0.) then
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
         write(unit=*,fmt='(a)') ' PIXRES is not a divisor of GAPSIZE.'
         write(unit=*,fmt='(a)') 'Adjusting PIXRES.'
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
         write(unit=*,fmt='(a,1x,f8.3)') ' GAPSIZE         = ',cci_gapsize
         write(unit=*,fmt='(a,1x,f8.3)') ' ORIGINAL PIXRES = ',cci_pixres
         cci_pixres = cci_gapsize / floor( cci_gapsize / cci_pixres)
         write(unit=*,fmt='(a,1x,f8.3)') ' UPDATED  PIXRES = ',cci_pixres
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
      end if
      !------------------------------------------------------------------------------------!



      !----- Initialise dimensions. -------------------------------------------------------!
      gap_nxy     = ceiling(sqrt(cci_gapmin/min_patch_area))
      rls_nxy     = gap_nxy    * nint(cci_gapsize/cci_pixres)
      rls_ngap    = gap_nxy    * gap_nxy
      rls_npixel  = rls_nxy    * rls_nxy
      rls_length  = rls_nxy    * cci_pixres
      cci_gaparea = cci_gapsize * cci_gapsize
      rls_area    = rls_length * rls_length
      gap_npixel  = nint(cci_gapsize/cci_pixres) * nint(cci_gapsize/cci_pixres)
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Allocate maps.                                                                 !
      !------------------------------------------------------------------------------------!
      allocate(rls_x     (rls_nxy,rls_nxy))
      allocate(rls_y     (rls_nxy,rls_nxy))
      allocate(rls_ztch  (rls_nxy,rls_nxy))
      allocate(rls_cci   (rls_nxy,rls_nxy))
      allocate(rls_igp   (rls_nxy,rls_nxy))
      allocate(rls_fbeam (rls_nxy,rls_nxy))
      allocate(gap_x0    (rls_ngap))
      allocate(gap_y0    (rls_ngap))
      allocate(gap_idx   (rls_ngap))
      allocate(gap_ipa   (rls_ngap))
      allocate(gap_mask  (rls_ngap))
      allocate(gap_fbeam (rls_ngap))
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !       Assign "raster" coordinates, and map each pixel to the corresponding gap.    !
      !------------------------------------------------------------------------------------!
      do j=1,rls_nxy
         do i=1,rls_nxy
            !----- Coordinates. -----------------------------------------------------------!
            rls_x(i,j) = (i - 1) * cci_pixres
            rls_y(i,j) = (j - 1) * cci_pixres
            !------------------------------------------------------------------------------!

            !----- Check to which gap this point corresponds to. --------------------------!
            rls_igp(i,j) = 1 + floor( rls_x(i,j) / cci_gapsize)                            &
                             + floor( rls_y(i,j) / cci_gapsize) * gap_nxy
            !------------------------------------------------------------------------------!
         end do
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Assign southwestern corner coordinates for all sites.                         !
      !------------------------------------------------------------------------------------!
      do i=1,rls_ngap
         gap_x0(i) = real(mod(i-1,gap_nxy)) * cci_gapsize
         gap_y0(i) = real((i-1)/gap_nxy)    * cci_gapsize
      end do
      !------------------------------------------------------------------------------------!


      return
   end subroutine init_cci_variables
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine resets any horizontal shading that might exist in the patch      !
   ! attribution.  Possibly unnecessary, definitely reassuring.                            !
   !---------------------------------------------------------------------------------------!
   subroutine reset_hrzshade(csite)
      use ed_state_vars, only : sitetype ! ! structure
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype), target :: csite ! Current site
      !----- Local variables. -------------------------------------------------------------!
      integer                :: ipa   ! Patch counter
      !------------------------------------------------------------------------------------!


      !----- Overwrite light information for all sites. -----------------------------------!
      do ipa=1,csite%npatches
         csite%fbeam     (ipa) = 1.0
         csite%light_type(ipa) = 1
      end do
      !------------------------------------------------------------------------------------!


      return
   end subroutine reset_hrzshade
   !---------------------------------------------------------------------------------------!


   !=======================================================================================!
   !=======================================================================================!
   !     This calculates the crown closure index and the correction factor to be applied   !
   ! to incoming direct radiation based on shading from neighbouring patches.  To find     !
   ! this, we develop a pseudo-landscape and determine the probability of a gap belonging  !
   ! to each patch to be shaded by neighbouring patches.                                   !
   !---------------------------------------------------------------------------------------!
   subroutine split_hrzshade(csite)
      use ed_state_vars         , only : sitetype            & ! structure
                                       , patchtype           & ! structure
                                       , allocate_sitetype   & ! subroutine
                                       , deallocate_sitetype & ! subroutine
                                       , copy_sitetype_mask  & ! subroutine
                                       , copy_sitetype       ! ! subroutine
      use canopy_radiation_coms , only : cci_radius          & ! intent(in)
                                       , cci_pixres          & ! intent(in)
                                       , cci_gapsize         & ! intent(in)
                                       , cci_gapmin          & ! intent(in)
                                       , cci_nretn           & ! intent(in)
                                       , rls_nxy             & ! intent(in)
                                       , rls_npixel          & ! intent(in)
                                       , rls_ngap            & ! intent(in)
                                       , rls_length          & ! intent(in)
                                       , rls_area            & ! intent(in)
                                       , rls_x               & ! intent(in)
                                       , rls_y               & ! intent(in)
                                       , rls_igp             & ! intent(in)
                                       , gap_x0              & ! intent(in)
                                       , gap_y0              & ! intent(in)
                                       , cci_gaparea         & ! intent(in)
                                       , gap_npixel          & ! intent(in)
                                       , at_bright           & ! intent(in)
                                       , at_dark             & ! intent(in)
                                       , rls_ztch            & ! intent(out)
                                       , rls_cci             & ! intent(out)
                                       , rls_fbeam           & ! intent(out)
                                       , gap_ipa             & ! intent(out)
                                       , gap_idx             & ! intent(out)
                                       , gap_mask            & ! intent(out)
                                       , gap_fbeam           ! ! intent(out)
      use allometry             , only : dbh2ca              & ! function
                                       , h2crownbh           ! ! function
      use random_utils          , only : fsample             & ! sub-routine
                                       , isample             & ! sub-routine
                                       , runif               & ! sub-routine
                                       , runif_sca           & ! function
                                       , ipickone            & ! function
                                       , fpickone            ! ! function
      use consts_coms           , only : pii                 & ! intent(in)
                                       , twopi               ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)              , target      :: csite       ! Current site
      !----- Local variables. -------------------------------------------------------------!
      type(sitetype)              , pointer     :: tsite       ! Scratch site
      type(patchtype)             , pointer     :: cpatch      ! Current patch
      integer     , dimension(:,:), allocatable :: ilight_ipa  ! Light type of subpatches
      integer     , dimension(:)  , allocatable :: gap_pool    ! Sequential
      integer     , dimension(:)  , allocatable :: ipa_ngaps   ! Sequential
      integer     , dimension(:)  , allocatable :: ipa_seq     ! Sequential
      integer     , dimension(:)  , allocatable :: ipa_shf     ! Shuffler
      integer                                   :: ico         ! Cohort counter
      integer                                   :: igp         ! Gap counter
      integer                                   :: iii         ! Aux counter
      integer                                   :: ipa         ! Patch counter
      integer                                   :: ipl         ! Individual plant counter
      integer                                   :: isf         ! Shuffle counter
      integer                                   :: ix          ! X index of point
      integer                                   :: iy          ! Y index of point
      integer                                   :: ngap_diff   ! Mismatch count
      integer                                   :: nindiv      ! Number of individuals
      integer                                   :: npa         ! Patch counter
      integer                                   :: npat_new    ! Number of new patches
      integer                                   :: npoints     ! Number of "returns"
      real(kind=4), dimension(:,:), allocatable :: cciarea_ipa ! Area by type of patch
      real(kind=4), dimension(:,:), allocatable :: fbeam_ipa   ! Absorption at the top
      real(kind=4)                              :: a_ptc       ! Angle of point
      real(kind=4)                              :: ca_ind      ! Crown area
      real(kind=4)                              :: gap_area    ! Total gap area
      real(kind=4)                              :: rh_ind      ! Crown horizontal radius
      real(kind=4)                              :: rh_ptc      ! Hor. radius of the point
      real(kind=4)                              :: rv_ind      ! Crown vertical radius
      real(kind=4)                              :: rv_ptc      ! Vert. radius of the point
      real(kind=4)                              :: wa_fbeam    ! Wgt. Avg. of "fbeam"
      real(kind=4)                              :: x_ind       ! X pos. of the individual
      real(kind=4)                              :: x_ptc       ! X position of point
      real(kind=4)                              :: y_ind       ! Y pos. of the individual
      real(kind=4)                              :: y_ptc       ! Y position of point
      real(kind=4)                              :: z_ptc       ! Z position of point
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !       Allocate temporary vectors.                                                  !
      !------------------------------------------------------------------------------------!
      allocate(ipa_shf    (  csite%npatches))
      allocate(ipa_seq    (  csite%npatches))
      allocate(ipa_ngaps  (  csite%npatches))
      allocate(fbeam_ipa  (3,csite%npatches))
      allocate(cciarea_ipa(3,csite%npatches))
      allocate(ilight_ipa (3,csite%npatches))
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Reset variables.                                                               !
      !------------------------------------------------------------------------------------!
      rls_ztch   (:,:) = 0.0
      rls_cci    (:,:) = 0.0
      rls_fbeam  (:,:) = 0.0
      gap_fbeam    (:) = 0.0
      fbeam_ipa  (:,:) = 0.0
      cciarea_ipa(:,:) = 0.0
      ilight_ipa (:,:) = 0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Create sequential vector, and shuffle it.  The latter will only be used in     !
      ! case rounding errors leave gaps unassigned or too many gaps assigned.  This should !
      ! be a minor adjustment only (~ 1-2 gaps).                                           !
      !------------------------------------------------------------------------------------!
      do ipa=1,csite%npatches
         ipa_seq(ipa) = ipa
      end do
      call isample(csite%npatches,ipa_seq,csite%npatches,ipa_shf,.false.)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Find the total number of gaps per patch (1st guess), then correct patches in  !
      ! case the total does not match the number of gaps.                                  !
      !------------------------------------------------------------------------------------!
      ipa_ngaps(:) = nint(csite%area(:) * rls_ngap)
      ngap_diff    = sum(ipa_ngaps) - rls_ngap
      if (ngap_diff > 0) then
         !---------------------------------------------------------------------------------!
         !     Loop through all patches in random order, adding one gap to the count until !
         ! all gaps have an associated patch.                                              !
         !---------------------------------------------------------------------------------!
         iii = 0
         isf = 0
         addgap_loop: do
            !------ Check whether we still have gaps to go. -------------------------------!
            if (iii >= ngap_diff) exit addgap_loop
            !------------------------------------------------------------------------------!


            !------ Update shuffling index and patch. -------------------------------------!
            isf = 1 + mod(isf,csite%npatches)
            ipa = ipa_shf(isf)
            !------------------------------------------------------------------------------!

            !------ Add one to the counter. -----------------------------------------------!
            ipa_ngaps(ipa) = ipa_ngaps(ipa) + 1
            iii            = iii + 1
            !------------------------------------------------------------------------------!
         end do addgap_loop
         !---------------------------------------------------------------------------------!
      elseif (ngap_diff < 0) then
         !---------------------------------------------------------------------------------!
         !     Loop through all patches in random order, deleting one gap to the count     !
         ! until we don't have more gap indices than gaps.                                 !
         !---------------------------------------------------------------------------------!
         iii = 0
         isf = 0
         delgap_loop: do
            !------ Check whether we still have gaps to go. -------------------------------!
            if (iii >= ngap_diff) exit delgap_loop
            !------------------------------------------------------------------------------!


            !------ Update shuffling index. -----------------------------------------------!
            isf = 1 + mod(isf,csite%npatches)
            ipa = ipa_shf(isf)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Subtract gap count only if it keep gap count above minimum.             !
            !------------------------------------------------------------------------------!
            if (ipa_ngaps(ipa) > nint(cci_gapmin)) then
               !------ Add one to the counter. --------------------------------------------!
               ipa_ngaps(ipa) = ipa_ngaps(ipa) - 1
               iii = iii + 1
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do delgap_loop
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !       Assign gap indices in order (we will shuffle afterwards).                    !
      !------------------------------------------------------------------------------------!
      ipa = 1
      iii = 0
      do igp=1,rls_ngap
         !----- Update indices. -----------------------------------------------------------!
         if (iii == ipa_ngaps(ipa)) then
            ipa = ipa + 1
            iii = 1
         else
            iii = iii + 1
         end if
         !---------------------------------------------------------------------------------!

         !----- Assign patch index. -------------------------------------------------------!
         gap_idx(igp) = ipa
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Shuffle gaps.                                                                  !
      !------------------------------------------------------------------------------------!
      call isample(rls_ngap,gap_idx,rls_ngap,gap_ipa,.false.)
      !----- Replace gap_idx by the gap indices. ------------------------------------------!
      do igp=1,rls_ngap
         gap_idx(igp) = igp
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Loop through all patches.                                                      !
      !------------------------------------------------------------------------------------!
      patchloop_1st: do ipa=1,csite%npatches
         !----- Get current patch. --------------------------------------------------------!
         cpatch => csite%patch(ipa)
         !---------------------------------------------------------------------------------!


         !---- Update mask and create subset of selectable gaps. --------------------------!
         gap_mask(:) = gap_ipa(:) == ipa
         allocate(gap_pool(ipa_ngaps(ipa)))
         gap_pool = pack(gap_idx,gap_mask)
         !---------------------------------------------------------------------------------!


         !----- Loop through all cohorts. -------------------------------------------------!
         cohloop: do ico=1,cpatch%ncohorts

            !----- Find horizontal and vertical radii. ------------------------------------!
            ca_ind = dbh2ca(cpatch%dbh(ico),cpatch%hite(ico),cpatch%sla(ico)               &
                           ,cpatch%pft(ico))
            rh_ind = sqrt(ca_ind * pii)
            rv_ind = 0.5 * ( cpatch%hite(ico) - h2crownbh(cpatch%hite(ico),cpatch%pft(ico)))
            !------------------------------------------------------------------------------!


            !----- Find out how many individuals to be placed in the landscape. -----------!
            nindiv = nint(cpatch%nplant(ico) * real(ipa_ngaps(ipa)) * cci_gaparea)
            !------------------------------------------------------------------------------!
            indloop: do ipl=1,nindiv
               !----- Find the position of this particular tree. --------------------------!
               igp    = ipickone(ipa_ngaps(ipa),gap_pool)
               x_ind  = gap_x0(igp) + runif_sca(0.,cci_gapsize)
               y_ind  = gap_y0(igp) + runif_sca(0.,cci_gapsize)
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !      Generate "point cloud returns" according to the crown area.          !
               !---------------------------------------------------------------------------!
               npoints = nint(cci_nretn * ca_ind)
               ptcloop: do iii=1,npoints
                  a_ptc           = runif_sca(0.,twopi)
                  rh_ptc          = rh_ind * (1.0 - runif_sca(0.,1.))
                  rv_ptc          = rv_ind * sqrt(1.0-rh_ptc*rh_ptc/(rh_ind*rh_ind))
                  x_ptc           = x_ind + rh_ptc * cos(a_ptc)
                  y_ptc           = y_ind + rh_ptc * sin(a_ptc)
                  z_ptc           = cpatch%hite(ico) - rv_ind + rv_ptc
                  ix              = 1 + mod(floor(x_ptc/cci_pixres),rls_nxy)
                  iy              = 1 + mod(floor(y_ptc/cci_pixres),rls_nxy)
                  !----- Update TCH height in case the point is higher. -------------------!
                  if (z_ptc > rls_ztch(ix,iy)) then
                     rls_ztch(ix,iy) = z_ptc
                  end if
                  !------------------------------------------------------------------------!
               end do ptcloop
               !---------------------------------------------------------------------------!
            end do indloop
            !------------------------------------------------------------------------------!
         end do cohloop
         !---------------------------------------------------------------------------------!

         !----- Free memory before moving to the next patch. ------------------------------!
         deallocate(gap_pool)
         !---------------------------------------------------------------------------------!
      end do patchloop_1st
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find crown closure index.                                                     !
      !------------------------------------------------------------------------------------!
      call cci_lieberman(rls_npixel,rls_x,rls_y,rls_ztch,rls_length,rls_cci,.true.)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find the uncorrected light illumination factor.                               !
      !------------------------------------------------------------------------------------!
      call  cci_abstop(rls_npixel,rls_cci,rls_fbeam)
      do iy=1,rls_nxy
         do ix=1,rls_nxy
            igp = rls_igp(ix,iy)
            gap_fbeam(igp) = gap_fbeam(igp) + rls_fbeam(ix,iy)
         end do
      end do
      gap_fbeam(:) = gap_fbeam(:) / gap_npixel
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Find the mean absorption at the top and area under each category.             !
      !------------------------------------------------------------------------------------!
      do igp=1,rls_ngap
         ipa = gap_ipa(igp)
         if (gap_fbeam(igp) > at_bright) then
            fbeam_ipa  (1,ipa) = fbeam_ipa  (1,ipa) + gap_fbeam(igp)
            cciarea_ipa(1,ipa) = cciarea_ipa(1,ipa) + 1.0
         else if (gap_fbeam(igp) > at_dark) then
            fbeam_ipa  (2,ipa) = fbeam_ipa  (2,ipa) + gap_fbeam(igp)
            cciarea_ipa(2,ipa) = cciarea_ipa(2,ipa) + 1.0
         else
            fbeam_ipa  (3,ipa) = fbeam_ipa  (3,ipa) + gap_fbeam(igp)
            cciarea_ipa(3,ipa) = cciarea_ipa(3,ipa) + 1.0
         end if
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Eliminate under-represented light environments.                                !
      !------------------------------------------------------------------------------------!
      do ipa=1,csite%npatches
         do iii=1,3
            if (cciarea_ipa(iii,ipa) > 0. .and. cciarea_ipa(iii,ipa) < cci_gapmin) then
               select case (iii)
               case (1)
                  !----- Find a patch to send data. ---------------------------------------!
                  if (cciarea_ipa(2,ipa) == 0.) then
                     !----- Merge bright and dark. ----------------------------------------!
                     fbeam_ipa  (3,ipa) = ( fbeam_ipa  (1,ipa) * cciarea_ipa(1,ipa)        &
                                          + fbeam_ipa  (3,ipa) * cciarea_ipa(3,ipa) )      &
                                        / ( cciarea_ipa(1,ipa) + cciarea_ipa(3,ipa) )
                     cciarea_ipa(3,ipa) = cciarea_ipa(1,ipa) + cciarea_ipa(3,ipa)
                     !---------------------------------------------------------------------!
                  else
                     !----- Merge bright and mid. -----------------------------------------!
                     fbeam_ipa  (2,ipa) = ( fbeam_ipa  (1,ipa) * cciarea_ipa(1,ipa)        &
                                          + fbeam_ipa  (2,ipa) * cciarea_ipa(2,ipa) )      &
                                        / ( cciarea_ipa(1,ipa) + cciarea_ipa(2,ipa) )
                     cciarea_ipa(2,ipa) = cciarea_ipa(1,ipa) + cciarea_ipa(2,ipa)
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!



                  !----- Remove information from bright patch. ----------------------------!
                  fbeam_ipa  (1,ipa) = 0.
                  cciarea_ipa(1,ipa) = 0.
                  !------------------------------------------------------------------------!
               case (2)
                  !----- Find a patch to send data. ---------------------------------------!
                  if (cciarea_ipa(1,ipa) == 0.) then
                     !----- Merge mid and dark. -------------------------------------------!
                     fbeam_ipa  (3,ipa) = ( fbeam_ipa  (2,ipa) * cciarea_ipa(2,ipa)        &
                                          + fbeam_ipa  (3,ipa) * cciarea_ipa(3,ipa) )      &
                                        / ( cciarea_ipa(2,ipa) + cciarea_ipa(3,ipa) )
                     cciarea_ipa(3,ipa) = cciarea_ipa(2,ipa) + cciarea_ipa(3,ipa)
                     !---------------------------------------------------------------------!
                  else
                     !----- Merge mid and bright. -----------------------------------------!
                     fbeam_ipa  (1,ipa) = ( fbeam_ipa  (2,ipa) * cciarea_ipa(2,ipa)        &
                                          + fbeam_ipa  (1,ipa) * cciarea_ipa(1,ipa) )      &
                                        / ( cciarea_ipa(2,ipa) + cciarea_ipa(1,ipa) )
                     cciarea_ipa(1,ipa) = cciarea_ipa(2,ipa) + cciarea_ipa(1,ipa)
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!



                  !----- Remove information from bright patch. ----------------------------!
                  fbeam_ipa  (2,ipa) = 0.
                  cciarea_ipa(2,ipa) = 0.
                  !------------------------------------------------------------------------!
               case (3)
                  !----- Find a patch to send data. ---------------------------------------!
                  if (cciarea_ipa(1,ipa) == 0.) then
                     !----- Merge mid and dark. -------------------------------------------!
                     fbeam_ipa  (2,ipa) = ( fbeam_ipa  (3,ipa) * cciarea_ipa(3,ipa)        &
                                          + fbeam_ipa  (2,ipa) * cciarea_ipa(2,ipa) )      &
                                        / ( cciarea_ipa(3,ipa) + cciarea_ipa(2,ipa) )
                     cciarea_ipa(2,ipa) = cciarea_ipa(3,ipa) + cciarea_ipa(2,ipa)
                     !---------------------------------------------------------------------!
                  else
                     !----- Merge bright and dark. ----------------------------------------!
                     fbeam_ipa  (1,ipa) = ( fbeam_ipa  (3,ipa) * cciarea_ipa(3,ipa)        &
                                          + fbeam_ipa  (1,ipa) * cciarea_ipa(1,ipa) )      &
                                        / ( cciarea_ipa(3,ipa) + cciarea_ipa(1,ipa) )
                     cciarea_ipa(1,ipa) = cciarea_ipa(3,ipa) + cciarea_ipa(1,ipa)
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!



                  !----- Remove information from bright patch. ----------------------------!
                  fbeam_ipa  (3,ipa) = 0.
                  cciarea_ipa(3,ipa) = 0.
                  !------------------------------------------------------------------------!
               end select
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Convert sum to mean.                                                           !
      !------------------------------------------------------------------------------------!
      where (cciarea_ipa > 0.0)
         fbeam_ipa   = fbeam_ipa / cciarea_ipa
      elsewhere
         fbeam_ipa   = 0.
      end where
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the light type index.                                                     !
      !------------------------------------------------------------------------------------!
      where (fbeam_ipa > at_bright)
         ilight_ipa   = 1
      elsewhere (fbeam_ipa > at_dark)
         ilight_ipa   = 2
      elsewhere (fbeam_ipa > 0.)
         ilight_ipa   = 3
      elsewhere
         ilight_ipa   = 0
      end where
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the area of all new patches.                                              !
      !------------------------------------------------------------------------------------!
      do ipa=1,csite%npatches
         cciarea_ipa(:,ipa) = csite%area(ipa) * cciarea_ipa(:,ipa)/sum(cciarea_ipa(:,ipa))
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Find the weighted mean of fbeam_ipa, and re-normalise so the weighted average !
      ! becomes always one.                                                                !
      !------------------------------------------------------------------------------------!
      wa_fbeam = 0.0
      do ipa=1,csite%npatches
         do iii=1,3
            wa_fbeam = wa_fbeam + fbeam_ipa(iii,ipa) * cciarea_ipa(iii,ipa)
         end do
      end do
      wa_fbeam       = wa_fbeam       / sum(cciarea_ipa)
      fbeam_ipa(:,:) = fbeam_ipa(:,:) / wa_fbeam
      !------------------------------------------------------------------------------------!



      !----- Allocate the temporary site that will host the original patches. -------------!
      nullify (tsite)
      allocate(tsite)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Allocate the new patch.  We will retain the order (tallest to shortest), and  !
      ! for the same original patch, we will organise it from brightest to darkest.        !
      !------------------------------------------------------------------------------------!
      npat_new = count(cciarea_ipa > 0.0)
      call allocate_sitetype(tsite,npat_new)
      npa = 0
      do ipa=1,csite%npatches
         do iii=1,3
            !----- Check that this candidate patch has non-zero area. ---------------------!
            if (cciarea_ipa(iii,ipa) > 0.) then
               npa = npa + 1
               !----- Copy all properties to the sub-patch. -------------------------------!
               call copy_sitetype(csite,tsite,ipa,ipa,npa,npa)
               !---------------------------------------------------------------------------!


               !----- Adjust patch area and assign illumination index. --------------------!
               tsite%area      (npa) = cciarea_ipa(iii,ipa)
               tsite%fbeam     (npa) = fbeam_ipa  (iii,ipa)
               tsite%light_type(npa) = ilight_ipa (iii,ipa)
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!


      !----- Copy temporary site back to the ED structure. --------------------------------!
      call deallocate_sitetype(csite)
      call allocate_sitetype(csite,npat_new)
      call copy_sitetype(tsite,csite,1,npat_new,1,npat_new)
      call deallocate_sitetype(tsite)
      !------------------------------------------------------------------------------------!


      !------ Free memory before leaving the sub-routine. ---------------------------------!
      deallocate(tsite      )
      deallocate(ipa_shf    )
      deallocate(ipa_seq    )
      deallocate(ipa_ngaps  )
      deallocate(fbeam_ipa  )
      deallocate(cciarea_ipa)
      deallocate(ilight_ipa )
      !------------------------------------------------------------------------------------!

      return
   end subroutine split_hrzshade
   !=======================================================================================!
   !=======================================================================================!



   !=======================================================================================!
   !=======================================================================================!
   !    Function cci.lieberman.                                                            !
   !    This function computes the crown closure index for each individual, following the  !
   ! method presented by                                                                   !
   !                                                                                       !
   !     Lieberman, M., D. Lieberman, R. Peralta, G. S. Hartshorn, 1995: Canopy closure    !
   !        and the distribution of tropical forest tree species at La Selva, Costa Rica.  !
   !        J. Trop. Ecol., 11 (2), 161--177.                                              !
   !---------------------------------------------------------------------------------------!
   subroutine cci_lieberman(nxyz,x,y,z,xymax,cci,circular)
      use canopy_radiation_coms, only : cci_radius  ! ! intent(in)
      use rk4_coms             , only : tiny_offset ! ! intent(in)
      implicit none
      !----- Variable declaration. --------------------------------------------------------!
      integer                      , intent(in)            :: nxyz
      real(kind=4), dimension(nxyz), intent(in)            :: x
      real(kind=4), dimension(nxyz), intent(in)            :: y
      real(kind=4), dimension(nxyz), intent(in)            :: z
      real(kind=4)                 , intent(in)            :: xymax
      real(kind=4), dimension(nxyz), intent(out)           :: cci
      logical                      , intent(in) , optional :: circular
      !----- Local variables. -------------------------------------------------------------!
      integer                                              :: m
      integer                                              :: n
      real(kind=8), dimension(nxyz)                        :: x8
      real(kind=8), dimension(nxyz)                        :: y8
      real(kind=8), dimension(nxyz)                        :: z8
      real(kind=8), dimension(nxyz)                        :: cci8
      real(kind=8)                                         :: xymax8
      real(kind=8)                                         :: radius8
      real(kind=8)                                         :: dx8
      real(kind=8)                                         :: dy8
      real(kind=8)                                         :: dz8
      real(kind=8)                                         :: dr8
      logical                                              :: circ_use
      !----- External function. -----------------------------------------------------------!
      real(kind=4)                 , external              :: sngloff
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     In case circular is not provided, assume boundary conditions are not circular. !
      !------------------------------------------------------------------------------------!
      if (present(circular)) then
         circ_use = circular
      else
         circ_use = .false.
      end if
      !------------------------------------------------------------------------------------!



      !----- Initialise double precision placeholders. ------------------------------------!
      x8(:)   = dble(x(:))
      y8(:)   = dble(y(:))
      z8(:)   = dble(z(:))
      cci8(:) = 0.d0
      radius8 = dble(cci_radius)
      xymax8  = dble(xymax)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Loop through each element, find CCI.                                           !
      !------------------------------------------------------------------------------------!
      oloop: do m=1,nxyz
         !---------------------------------------------------------------------------------!
         !     Check every grid element.                                                   !
         !---------------------------------------------------------------------------------!
         iloop: do n=1,nxyz
            !------------------------------------------------------------------------------!
            !     Find distances.  In case the landscape is assumed circular, we must      !
            ! check distances with an offset of the loop.                                  !
            !------------------------------------------------------------------------------!
            if (circ_use) then
               dx8 = min( abs(x8(n) - x8(m)          )                                     &
                        , abs(x8(n) - x8(m) + xymax8 )                                     &
                        , abs(x8(n) - x8(m) - xymax8 ) )
               dy8 = min( abs(y8(n) - y8(m)          )                                     &
                        , abs(y8(n) - y8(m) + xymax8 )                                     &
                        , abs(y8(n) - y8(m) - xymax8 ) )
            else
               dx8 = x8(n) - x8(m)
               dy8 = y8(n) - y8(m)
            end if
            dz8 = x8(n) - x8(m)
            dr8 = sqrt(dx8*dx8 + dy8*dy8)
            !------------------------------------------------------------------------------!


            !----- Check whether point is within radius and taller than current point. ----!
            if (dr8 > 0.d0 .and. dr8 <= radius8 .and. dz8 > 0.d0) then
               cci8(m) = cci8(m) + sin(dz8 / sqrt(dz8*dz8 + dr8*dr8))
            end if
            !------------------------------------------------------------------------------!
         end do iloop
         !---------------------------------------------------------------------------------!
      end do oloop
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Loop through each element, find CCI.                                           !
      !------------------------------------------------------------------------------------!
      do m=1,nxyz
         cci(m) = sngloff(cci8(m),tiny_offset)
      end do
      !------------------------------------------------------------------------------------!

      return
   end subroutine cci_lieberman
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      Find the corresponding light absorption factor.  The function and coefficients   !
   ! are deliberately hard-coded because this will change in the near future.              !
   !---------------------------------------------------------------------------------------!
   subroutine cci_abstop(nxyz,cci,abstop)
      use canopy_radiation_coms, only : at08        & ! intent(in)
                                      , at18        ! ! intent(in)
      use rk4_coms             , only : tiny_offset ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      integer                      , intent(in)  :: nxyz
      real(kind=4), dimension(nxyz), intent(in)  :: cci
      real(kind=4), dimension(nxyz), intent(out) :: abstop
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)                               :: cci8
      real(kind=8)                               :: abstop8
      integer                                    :: n
      !----- External function. -----------------------------------------------------------!
      real(kind=4)                 , external    :: sngloff
      !------------------------------------------------------------------------------------!


      !----- Translate CCI into absorption at the top of the canopy. ----------------------!
      do n=1,nxyz
         cci8       = dble(cci(n))
         abstop8    = exp(at08 + at18 * cci8)
         abstop (n) = sngloff(abstop8,tiny_offset)
      end do
      !------------------------------------------------------------------------------------!
   end subroutine cci_abstop
   !---------------------------------------------------------------------------------------!
end module hrzshade_utils
!==========================================================================================!
!==========================================================================================!
