!==========================================================================================!
!==========================================================================================!
!     This subroutine takes space-time averaged precipitation, and performs downscaling    !
! techniques to resolve a realization of point precipitation at a smaller time step.  The  !
! resulting precipitation is only a realization, because it is achieved through a          !
! stochastic process.  It also relies on parameters which describe the probability density !
! functions.  These parameters may vary by time and space with different dependencies.     !
! For practical application, we will use region averaged parameters which vary by season   !
! and hour of day.  This follows Lammering (1999) routines that use multidimensional       !
! bounded cascades following Menabde, Seed and others may be added in the future.          !
!                                                                                          !
!     The subroutine accepts 6 arguments:                                                  !
!  - nprecip_in  : Number of points of the input precipitation.                            !
!  - nprecip_out : Number of points of the output precipitation.                           !
!  - imo         : Current month (1=jan, 2=feb...)                                         !
!  - ireali      : Current realisation (just for the banner printing...                    !
!  - precip_in   : Vector containing the input precipitation data.                         !
!  - precip_out  : Vector containing the output precipitation data.                        !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine downscale_precip(nprecip_in,nprecip_out,imo,ireali,precip_in,precip_out)
   use rconstants, only : day_sec          ! ! intent(in)
   use mod_ioopts, only : ninputs_day      & ! intent(in)
                        , inpfrq           & ! intent(in)
                        , radfrq           & ! intent(in)
                        , radratio         ! ! intent(in)
   use mod_interp, only : interp_buffer    & ! intent(inout)
                        , npdf             & ! intent(in)
                        , nlocpcp          & ! intent(in)
                        , nlocpcpi         & ! intent(in)
                        , max_local_precip & ! intent(in)
                        , local_precip     & ! intent(in)
                        , frac_u           ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)  :: nprecip_in
   integer                     , intent(in)  :: nprecip_out
   integer                     , intent(in)  :: imo
   integer                     , intent(in)  :: ireali
   real, dimension(nprecip_in) , intent(in)  :: precip_in
   real, dimension(nprecip_out), intent(out) :: precip_out
   !----- Local variables. ----------------------------------------------------------------!
   integer                  :: imet
   integer                  :: ic
   integer                  :: ip
   integer                  :: isub
   real, dimension(nlocpcp) :: lcdf       ! The local cum. dens. fctn of precip
   real                     :: scdf       ! Sum of cdf, used for normalization
   real                     :: cdfy       ! Random number on the y-axis of the cdf
   real                     :: correction ! This is applied after each month resampling, 
                                          !     to conserve volume
   real                     :: frac_ug    ! The native grid resolution rainfall fraction
   !---------------------------------------------------------------------------------------!
   !     Hour Of Day Index.  This is the index of the input frequency block.  For example, !
   ! if the input data is every six hours, then the index should be 1-4, for hours 1-6,    !
   ! 7-12, 13-18, and 19-24.                                                               !
   !---------------------------------------------------------------------------------------!
   integer                  :: hodi
   !----- Locally saved variables. --------------------------------------------------------!
   logical, save            :: first_time = .true.
   !---------------------------------------------------------------------------------------!


   !----- This part needs to be called only once. -----------------------------------------!
   if (first_time) then

      !----- Filling the local precipitation vector. --------------------------------------!
      do ic=1,nlocpcp
         local_precip(ic) = max_local_precip*real(ic)/real(nlocpcp)
      end do

      !----- Copying the namelist variable into the array-based variable. -----------------!
      interp_buffer%frac_u = reshape(frac_u(1:npdf), (/ ninputs_day, 12 /) )

      first_time = .false.
   end if
   !---------------------------------------------------------------------------------------!


   write (unit=*,fmt='(a,1x,i5,a)')                                                        &
                              '     - Running the precipitation downscaling, ',ireali,'...'

   !----- Initialize output. --------------------------------------------------------------!
   precip_out = 0.0

   if (sum(precip_in) > 0.) then
 
      !----- Calculate the grid rainfall precipitation fraction. --------------------------!
      frac_ug=0.0
      do imet=1,nprecip_in
         frac_ug  = frac_ug + ceiling(precip_in(imet)/100.0)
      end do
      frac_ug=frac_ug / real(nprecip_in)

      !----- Keep looping until downscaled is more than zero. -----------------------------!
      checksum: do

         !----- This is the loop where the reassignment is conducted. ---------------------!
         metloop: do imet = 1,nprecip_in

            if (precip_in(imet) > 0.0) then

               hodi = 1+floor(mod(real((imet-1)*inpfrq),day_sec)/inpfrq)

               !----- Construct the cumulative density function. --------------------------!
               scdf=0.
               do ic=1,nlocpcp
                  lcdf(ic) = interp_buffer%frac_u(hodi,imo) / precip_in(imet)              &
                           * exp(max( -interp_buffer%frac_u(hodi,imo)*local_precip(ic)     &
                                    / precip_in(imet), -80.0 ))

                  scdf = scdf + lcdf(ic)
               end do

               lcdf(1) = lcdf(1) / scdf
               do ic=2,nlocpcp
                  lcdf(ic) = lcdf(ic-1) + lcdf(ic)/scdf
               end do
               lcdf(nlocpcp)=1.0

               !----- Randomly select from this distribution for each sub-step. -----------!
               subloop1: do isub=1,radratio
                  ip = (imet-1)*radratio+isub
                  call random_number(cdfy)
                  
                  if(cdfy > (1.0-min(1.0,interp_buffer%frac_u(hodi,imo)/frac_ug))) then
                     call random_number(cdfy)
                     ic = 0
                     selloop: do 
                        ic = ic + 1
                        if (lcdf(ic) > cdfy) exit selloop
                     end do selloop

                     precip_out(ip) = local_precip(ic)

                  else
                     precip_out(ip) = 0.0

                  end if
               end do subloop1
            else
               subloop2: do isub=1,radratio
                  ip = (imet-1) * radratio + isub
                  precip_out(ip) = 0.0
               end do subloop2
            end if
         end do metloop
         if ( sum(precip_out) > 0. ) exit checksum
      end do checksum

      !----- Apply correction to conserve monthly rainfall volume. ------------------------!
      correction = radratio*sum(precip_in)/sum(precip_out)
      do imet = 1,nprecip_in
         do isub=1,radratio
            ip             = (imet-1) * radratio + isub
            precip_out(ip) = precip_out(ip) * correction
         end do
      end do
   end if  ! Check on positive precip

   return

end subroutine downscale_precip
!==========================================================================================!
!==========================================================================================!
