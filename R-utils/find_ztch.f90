!==========================================================================================!
!==========================================================================================!
!    Function find.ztch.                                                                   !
!                                                                                          !
!    This function estimates the canopy height model given the height, crown area, and     !
! crown depth, and geolocation of each individual.  It returns both the canopy height      !
! model and the list of trees that are partially or entirely at the canopy.                !
!------------------------------------------------------------------------------------------!
subroutine find_ztch(nind,nxtch,nytch,nseed,nretn,xyres,seed,height,rh,rv,x0,y0,ztch,canopy)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                             , intent(in)   :: nind
   integer                             , intent(in)   :: nxtch
   integer                             , intent(in)   :: nytch
   integer                             , intent(in)   :: nseed
   integer                             , intent(in)   :: nretn
   real(kind=8)                        , intent(in)   :: xyres
   integer     , dimension(nseed)      , intent(in)   :: seed
   real(kind=8), dimension(nind)       , intent(in)   :: height
   real(kind=8), dimension(nind)       , intent(in)   :: rh
   real(kind=8), dimension(nind)       , intent(in)   :: rv
   real(kind=8), dimension(nind)       , intent(in)   :: x0
   real(kind=8), dimension(nind)       , intent(in)   :: y0
   real(kind=8), dimension(nxtch,nytch), intent(out)  :: ztch
   integer     , dimension(nind)       , intent(out)  :: canopy
   !----- Local variables. ----------------------------------------------------------------!
   integer                                            :: nseed_def
   integer                                            :: nseed_use
   integer                                            :: ind
   integer                                            :: pts
   integer                                            :: npts
   integer                                            :: ix
   integer                                            :: iy
   real(kind=8)                                       :: azim
   real(kind=8)                                       :: a_ptc
   real(kind=8)                                       :: rh_ptc
   real(kind=8)                                       :: rv_ptc
   real(kind=8)                                       :: x_ptc
   real(kind=8)                                       :: y_ptc
   real(kind=8)                                       :: z_ptc
   real(kind=8)                                       :: x_max
   real(kind=8)                                       :: y_max
   real(kind=4)                                       :: rnow
   integer     , dimension(nxtch,nytch)               :: itch
   !----- Local constants. ----------------------------------------------------------------!
   real(kind=8)                        , parameter    :: pi = 3.141592653589793d0
   !---------------------------------------------------------------------------------------!


   !----- Set the random seed. ------------------------------------------------------------!
   call random_seed(size=nseed_def)
   nseed_use = min(nseed,nseed_def)
   call random_seed(put=seed(1:nseed_use))
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Initialise output variables.                                                      !
   !---------------------------------------------------------------------------------------!
   ztch  (:,:) = 0.d0
   canopy  (:) = 0
   !---------------------------------------------------------------------------------------!


   !----- Define maximum dimensions of the TCH model. -------------------------------------!
   x_max = dble(nxtch) * xyres
   y_max = dble(nytch) * xyres
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Loop through individuals, and randomly assign points to represent the canopy.     !
   !---------------------------------------------------------------------------------------!
   ind_loop: do ind=1,nind
      npts = nint(nretn * pi * rh(ind) * rh(ind))
   
      ptc_loop: do pts=1,npts
         call random_number(rnow)
         a_ptc  = 2.d0 * pi * dble(rnow)
         call random_number(rnow)
         rh_ptc = rh(ind) * (1.d0 - dble(rnow))
         rv_ptc = rv(ind) * sqrt(1.d0-rh_ptc*rh_ptc/(rh(ind)*rh(ind)))
         x_ptc  = modulo(x0(ind) + rh_ptc * cos(a_ptc),x_max)
         y_ptc  = modulo(y0(ind) + rh_ptc * sin(a_ptc),y_max)
         z_ptc  = height(ind) - rv(ind) + rv_ptc
         ix     = 1 + modulo(floor(x_ptc/xyres),nxtch)
         iy     = 1 + modulo(floor(y_ptc/xyres),nytch)
         !---------------------------------------------------------------------------------!
         !    Update TCH height in case the point is the highest recorded.                 !
         !---------------------------------------------------------------------------------!
         if (z_ptc > ztch(ix,iy)) then
            ztch(ix,iy) = z_ptc
            itch(ix,iy) = ind
         end if
         !---------------------------------------------------------------------------------!
      end do ptc_loop
   end do ind_loop
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Flag all the trees that are partially or totally in the canopy.                  !
   !---------------------------------------------------------------------------------------!
   can_loop: do ind=1,nind
      if (any(itch == ind)) canopy(ind) = 1
   end do can_loop


   return
end subroutine find_ztch
!==========================================================================================!
!==========================================================================================!
