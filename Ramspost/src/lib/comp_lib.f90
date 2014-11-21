!==========================================================================================!
!==========================================================================================!
!     The following subroutines used to be "entries", which is deprecated.  They have been !
! superseded  by subroutines.                                                              !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_zero(n1,n2,n3,a)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=0.
         end do
      end do
   end do
   return
end subroutine RAMS_comp_zero
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_1minus(n1,n2,n3,a)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!

   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=1.-a(i,j,k)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_1minus
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_mults(n1,n2,n3,a,s)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   real(kind=4)                     , intent(in)    :: s
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k) * s
         end do
      end do
   end do
   return
   return
end subroutine RAMS_comp_mults
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_adds(n1,n2,n3,a,s)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   real(kind=4)                     , intent(in)    :: s
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!

   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k) + s
         end do
      end do
   end do
   return
end subroutine RAMS_comp_adds
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_accum(n1,n2,n3,a,b)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   real(kind=4), dimension(n1,n2,n3), intent(in)    :: b
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!

   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)+b(i,j,k)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_accum
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_noneg(n1,n2,n3,a)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!

   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=max(a(i,j,k),0.)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_noneg
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_nopos(n1,n2,n3,a)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!
   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=min(a(i,j,k),0.)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_nopos
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_subt(n1,n2,n3,a,b)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   real(kind=4), dimension(n1,n2,n3), intent(in)    :: b
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!

   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)-b(i,j,k)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_subt
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_mult(n1,n2,n3,a,b)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   real(kind=4), dimension(n1,n2,n3), intent(in)    :: b
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!

   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)*b(i,j,k)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_mult
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_z(n1,n2,n3,geo,topt,ngrd)
   use somevars, only : myztn  & ! intent(in)
                      , myzmn  & ! intent(in)
                      , mynnzp ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: geo
   real(kind=4), dimension(n1,n2)   , intent(in)    :: topt
   integer                          , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!

   do k=1,n3
      do j=1,n2
         do i=1,n1
            geo(i,j,k) = topt(i,j) + myztn(k,ngrd)*(1.-topt(i,j)/myzmn(mynnzp(1)-1,1))
         end do
      end do
   end do
   return
end subroutine RAMS_comp_z
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_tempk(n1,n2,n3,inth_outt,exner)
   use therm_lib, only : extheta2temp ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: inth_outt
   real(kind=4), dimension(n1,n2,n3), intent(in)    :: exner
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!

   do k=1,n3
      do j=1,n2
         do i=1,n1
            inth_outt(i,j,k)= extheta2temp(exner(i,j,k),inth_outt(i,j,k))
         end do
      end do
   end do
   return
end subroutine RAMS_comp_tempk
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_press(n1,n2,n3,inex_outp)
   use therm_lib, only : exner2press ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: inex_outp
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!

   do k=1,n3
      do j=1,n2
         do i=1,n1
            inex_outp(i,j,k) = exner2press(inex_outp(i,j,k)) * .01
         end do
      end do
   end do
   return
end subroutine RAMS_comp_press
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_avgw(n1,n2,n3,a)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!
   do k=n3,2,-1
      do j=1,n2
         do i=1,n1
            a(i,j,k)=0.5*(a(i,j,k)+a(i,j,k-1))
         end do
      end do
   end do
   return
end subroutine RAMS_comp_avgw
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_avgu(n1,n2,n3,a)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!
   do k=1,n3
      do j=1,n2
         do i=n1,2,-1
            a(i,j,k)=0.5*(a(i,j,k)+a(i-1,j,k))
         end do
      end do
   end do
   return
end subroutine RAMS_comp_avgu
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_avgv(n1,n2,n3,a)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!

   do k=1,n3
      do j=n2,2,-1
         do i=1,n1
            a(i,j,k)=0.5*(a(i,j,k)+a(i,j-1,k))
         end do
      end do
   end do
   return
end subroutine RAMS_comp_avgv
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_speed(n1,n2,n3,inu_outmag,v)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: inu_outmag
   real(kind=4), dimension(n1,n2,n3), intent(in)    :: v
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!
   do k=1,n3
      do j=1,n2
         do i=1,n1
            inu_outmag(i,j,k)=sqrt(inu_outmag(i,j,k)**2+v(i,j,k)**2)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_speed
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_dir(n1,n2,n3,inu_outdir,v,ngrd)
   use somevars, only : myplatn & ! intent(in)
                      , myplonn & ! intent(in)
                      , myxtn   & ! intent(in)
                      , myytn   ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: inu_outdir
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: v
   integer                          , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   real                                             :: qlat
   real                                             :: qlon
   real                                             :: ff
   real                                             :: unow
   real                                             :: vnow
   !---------------------------------------------------------------------------------------!

   do k=1,n3
      do j=1,n2
         do i=1,n1
            call xy_ll(qlat,qlon,myplatn(ngrd),myplonn(ngrd),myxtn(i,ngrd),myytn(j,ngrd))
            unow = inu_outdir(i,j,k)
            vnow = v(i,j,k)
            call uvtoueve(unow,vnow,inu_outdir(i,j,k),v(i,j,k),qlat,qlon                   &
                         ,myplatn(ngrd),myplonn(ngrd))
            call winddf(inu_outdir(i,j,k),ff,inu_outdir(i,j,k),v(i,j,k))
         end do
      end do
   end do
   return
end subroutine RAMS_comp_dir
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_dewk(n1,n2,n3,a,b,c)
   use therm_lib, only : exner2press   & ! function
                       , extheta2temp  & ! function
                       , rslif         & ! function
                       , dewfrostpoint ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: b
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: c
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   real                                             :: xpress
   real                                             :: xtemp
   real                                             :: xwatsat
   !---------------------------------------------------------------------------------------!

   do k=1,n3
      do j=1,n2
         do i=1,n1
            xpress=exner2press(b(i,j,k))
            xtemp=extheta2temp(b(i,j,k),c(i,j,k))
            xwatsat=rslif(xpress,xtemp)
            a(i,j,k)=dewfrostpoint(xpress,min(a(i,j,k),xwatsat) )
         end do
      end do
   end do
   return
end subroutine RAMS_comp_dewk
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_thetv(n1,n2,n3,a,b,c)
   use therm_lib, only : virtt   ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: b
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: c
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!

   do k=1,n3
      do j=1,n2
         do i=1,n1
            a(i,j,k)=virtt(a(i,j,k),b(i,j,k),c(i,j,k))
         end do
      end do
   end do
   return
end subroutine RAMS_comp_thetv
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_rh(n1,n2,n3,a,b,c)
   use therm_lib, only : extheta2temp   & ! function
                       , exner2press    & ! function
                       , rehuil         ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: b
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: c
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   real                                             :: xpress
   real                                             :: xtemp
   !---------------------------------------------------------------------------------------!
   do k=1,n3
      do j=1,n2
         do i=1,n1
            xtemp    = extheta2temp(b(i,j,k),c(i,j,k))
            xpress   = exner2press(b(i,j,k))
            a(i,j,k) = 100. * min(1.,rehuil(xpress,xtemp,a(i,j,k),.false.))
         end do
      end do
   end do
   return
end subroutine RAMS_comp_rh
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_vegclass(nmax,a)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                      , intent(in)    :: nmax
   real(kind=4), dimension(nmax), intent(inout) :: a
   !----- Local variables. ----------------------------------------------------------------!
   integer                                      :: n
   !---------------------------------------------------------------------------------------!
   do n=1,nmax
      a(n) = real(nint(a(n)))
   end do
   return
end subroutine RAMS_comp_vegclass
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_raintemp(n1,n2,n3,a)
   use rconstants, only : tsupercool_liq   & ! intent(in)
                        , cliqi            & ! intent(in)
                        , t00              ! ! intent(in)
   use rout_coms , only : undefflg         ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   real                                             :: xpress
   real                                             :: xtemp
   !---------------------------------------------------------------------------------------!
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if (a(i,j,k) > 0.) then
               a(i,j,k) = tsupercool_liq + a(i,j,k) * cliqi - t00
            else
               a(i,j,k) = undefflg
            end if
         end do
      end do
   end do
   return
end subroutine RAMS_comp_raintemp
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_qtcpcp(n1,n2,n3,a)
   use therm_lib , only : uint2tl  ! ! function
   use rconstants, only : t00      ! ! intent(in)
   use rout_coms , only : undefflg ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   real                                             :: temptemp
   real                                             :: fracliq
   !---------------------------------------------------------------------------------------!
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if (a(i,j,k) > 0.) then
               call uint2tl(a(i,j,k),temptemp,fracliq)
               a(i,j,k) = temptemp - t00
            else
               a(i,j,k) = undefflg
            end if
         end do
      end do
   end do
   return
end subroutine RAMS_comp_qtcpcp
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_fracliq(n1,n2,n3,a)
   use therm_lib , only : uint2tl  ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   real                                             :: temptemp
   real                                             :: fracliq
   !---------------------------------------------------------------------------------------!
   do k=1,n3
      do j=1,n2
         do i=1,n1
            call uint2tl(a(i,j,k),temptemp,fracliq)
            a(i,j,k) = fracliq
         end do
      end do
   end do
return
end subroutine RAMS_comp_fracliq
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_fracice(n1,n2,n3,a)
   use therm_lib , only : uint2tl  ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   real                                             :: temptemp
   real                                             :: fracliq
   !---------------------------------------------------------------------------------------!
   do k=1,n3
      do j=1,n2
         do i=1,n1
            call uint2tl(a(i,j,k),temptemp,fracliq)
            a(i,j,k) = 1.0 - fracliq
         end do
      end do
   end do
   return
end subroutine RAMS_comp_fracice
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_iceliq(n1,n2,n3,energy,mass,m_ice,m_liq)
   use therm_lib , only : uint2tl  ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(in)    :: energy
   real(kind=4), dimension(n1,n2,n3), intent(in)    :: mass
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: m_ice
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: m_liq
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   real                                             :: temptemp
   real                                             :: fracliq
   !---------------------------------------------------------------------------------------!
   do k=1,n3
      do j=1,n2
         do i=1,n1
            call uint2tl(energy(i,j,k),temptemp,fracliq)
            m_ice(i,j,k) = ( 1.0 - fracliq ) * mass(i,j,k)
            m_liq(i,j,k) =         fracliq   * mass(i,j,k)
         end do
      end do
   end do
return
end subroutine RAMS_comp_iceliq
!..........................................................................................!


!..........................................................................................!
subroutine RAMS_comp_hydrodiam(n1,n2,n3,a,c,ccfmas,ppwmas)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: c
   real                             , intent(in)    :: ccfmas
   real                             , intent(in)    :: ppwmas
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   real                                             :: rpwmas
   !----- Local constants. ----------------------------------------------------------------!
   real                             , parameter     :: negligible = 1.e-10
   !---------------------------------------------------------------------------------------!
   rpwmas = 1. / ppwmas
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if(a(i,j,k) > negligible .and. c(i,j,k) > negligible)then
               a(i,j,k) = (a(i,j,k) / (c(i,j,k) * ccfmas))**rpwmas
            else
               a(i,j,k) = 0.
            endif
         end do
      end do
   end do
   return
end subroutine RAMS_comp_hydrodiam
!..........................................................................................!


!..........................................................................................!
subroutine rams_sum_snowlayers(nx,ny,nz,np,a)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                             , intent(in)    :: nx
   integer                             , intent(in)    :: ny
   integer                             , intent(in)    :: nz
   integer                             , intent(in)    :: np
   real(kind=4), dimension(nx,ny,nz,np), intent(inout) :: a
   !----- Local variables. ----------------------------------------------------------------!
   integer                                             :: i
   integer                                             :: j
   integer                                             :: k
   integer                                             :: p
   !---------------------------------------------------------------------------------------!
   do i=1,nx
      do j=1,ny
         do k=2,nz
            do p=1,np
               a(i,j,1,p) = a(i,j,1,p) + a(i,j,k,p)
            end do
         end do
      end do
   end do
   return
end subroutine rams_sum_snowlayers
!..........................................................................................!


!..........................................................................................!
subroutine rams_fill_sst(n1,n2,n3,kp,a,c)
   use therm_lib , only : uint2tl ! ! function
   use rconstants, only : t00     ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2)   , intent(inout) :: a
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: c
   integer                          , intent(in)    :: kp
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   real                                             :: temptemp
   real                                             :: fracliq
   !---------------------------------------------------------------------------------------!
   do j=1,n2
      do i = 1,n1
         call uint2tl(c(i,j,kp),temptemp,fracliq)
         a(i,j) = temptemp-t00
      end do
   end do
   return
end subroutine rams_fill_sst
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     Patch multiplication.                                                                !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_multap(n1,n2,n3,n4,a,b)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   integer                          , intent(in)    :: n4
   real(kind=4), dimension(n1,n2,n4), intent(inout) :: a
   real(kind=4), dimension(n1,n2,n3), intent(in)    :: b
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   !---------------------------------------------------------------------------------------!

   do k=1,n4
      do j=1,n2
         do i=1,n1
            a(i,j,k)=a(i,j,k)*b(i,j,1)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_multap
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine get_leaf_soil(n1,n2,n4,n5,a,a2)
   implicit none
   integer, intent(in) :: n1,n2,n4,n5
   real, dimension(n1,n2,n4,n5), intent(out) :: a2
   real, dimension(n1,n2,n4*n5), intent(in)  :: a
   integer :: kip, k,i,j,ip
   kip=0
   do ip=1,n5
      do k=1,n4
         kip=kip+1
         do j=1,n2
            do i=1,n1
               a2(i,j,k,ip)=a(i,j,kip)
            end do
         end do
      end do
   end do
   return
end subroutine get_leaf_soil
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine converts a 3-D array into a cloud-dependent, 4-D array.             !
!------------------------------------------------------------------------------------------!
subroutine get_cumulus(n1,n2,n3,n6,a,a6)
   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   integer                     , intent(in)  :: n1
   integer                     , intent(in)  :: n2
   integer                     , intent(in)  :: n3
   integer                     , intent(in)  :: n6
   real, dimension(n1,n2,n3,n6), intent(out) :: a6
   real, dimension(n1,n2,n3*n6), intent(in)  :: a
   !----- Local variables. ----------------------------------------------------------------!
   integer                                   :: kic
   integer                                   :: k
   integer                                   :: i
   integer                                   :: j
   integer                                   :: ic
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   kic=0
   do ic=1,n6
      do k=1,n3
         kic=kic+1
         do j=1,n2
            do i=1,n1
               a6(i,j,k,ic)=a(i,j,kic)
            end do
         end do
      end do
   end do
   return
end subroutine get_cumulus
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_richardson(n1,n2,n3,np,rib,z0,speed,thetav_atm,thetav_can,topt,ngrd)

   use somevars
   use rconstants
   use therm_lib, only : uint2tl, uextcm2tl, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3,np
   real   , dimension(n1,n2,np), intent(inout) :: rib
   real   , dimension(n1,n2,np), intent(in)    :: z0,thetav_can
   real   , dimension(n1,n2,n3), intent(in)    :: speed,thetav_atm
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   integer                                     :: i,j,p
   real                                        :: zedtop,zagl,spd

   zedtop = myzmn(mynnzp(1)-1,1)
   do j=1,n2
      do i=1,n1
         zagl=myztn(2,ngrd)*(1.-topt(i,j)/zedtop)
         do p=1,np
            spd = max(speed(i,j,2),0.65)
            rib(i,j,p) = grav * (zagl-z0(i,j,p)) * (thetav_atm(i,j,2)-thetav_can(i,j,p))   &
                       / (0.5 * (thetav_can(i,j,2)+thetav_can(i,j,p)) * spd**2)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_richardson
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_dn0(n1,n2,n3,a,b,c,topt,ngrd)

   use somevars
   use rconstants
   use therm_lib, only : uint2tl, uextcm2tl, dewfrostpoint, rslif, virtt, thetaeiv
   implicit none
   integer                     , intent(in)    :: n1,n2,n3
   real   , dimension(n1,n2,n3), intent(inout) :: a,b,c
   real   , dimension(n1,n2)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   integer                                     :: i,j,k
   real   , dimension(n3)                      :: scratch2,scratch11,scratch12
   real                                        :: zedtop,c1,c2,c3

   zedtop = myzmn(mynnzp(1)-1,1)
   do j=1,n2
      do i=1,n1
         do k=1,n3
           scratch2(k)=myztn(k,ngrd)*(1.-topt(i,j)/zedtop)+topt(i,j)
         enddo
         call htint(n3,mypi01dn(1,ngrd),myztn(1,ngrd),n3,scratch11,scratch2)
         call htint(n3,myth01dn(1,ngrd),myztn(1,ngrd),n3,scratch12,scratch2)       
         do k=1,n3
            b(i,j,k)=scratch12(k)
         enddo
         a(i,j,n3) = scratch11(n3)

         c1=grav*2.*(1.-topt(i,j)/zedtop)
         c2=(1-cpor)
         c3=cpdry**c2
         do k=n3-1,1,-1
            a(i,j,k)=a(i,j,k+1) +c1/((b(i,j,k)+b(i,j,k+1))*mydzmn(k,ngrd))
         enddo
         do k=1,n3
            c(i,j,k)=(c3*p00)/(rdry*b(i,j,k)*a(i,j,k)**c2)
         enddo

      enddo
   enddo
   return
end subroutine RAMS_comp_dn0
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine calculates the zonal component of the relative vorticity.           !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_relvortx(nx,ny,nz,curlx,vwnd,wwnd,dwdy,dvdz,aux,topt,ngrd)

   use somevars  , only : myxmn       & ! intent(in)
                        , myxtn       & ! intent(in)
                        , myymn       & ! intent(in)
                        , myytn       & ! intent(in)
                        , myzmn       & ! intent(in)
                        , myztn       & ! intent(in)
                        , mydeltayn   & ! intent(in)
                        , mydzmn      & ! intent(in)
                        , mydztn      & ! intent(in)
                        , mynnzp      & ! intent(in)
                        , myjdim      & ! intent(in)
                        , myihtran    ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(inout) :: curlx
   real   , dimension(nx,ny,nz), intent(inout) :: vwnd
   real   , dimension(nx,ny,nz), intent(inout) :: wwnd
   real   , dimension(nx,ny,nz), intent(inout) :: dvdz
   real   , dimension(nx,ny,nz), intent(inout) :: dwdy
   real   , dimension(nx,ny,nz), intent(inout) :: aux
   real   , dimension(nx,ny)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   integer                                     :: y1
   integer                                     :: y2
   integer                                     :: z1
   integer                                     :: z2
   real                                        :: factor
   real   , dimension(nx+ny+nz)                :: dum1
   real   , dimension(nx+ny+nz)                :: dum2
   !---------------------------------------------------------------------------------------!



   !----- Bottom boundary condition for meridional wind. ----------------------------------!
   factor = myztn(1,ngrd) / myztn(2,ngrd)
   do y=1,ny
      do x=1,nx
         vwnd(x,y,1) = vwnd(x,y,2) * factor
      end do
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Initialise the auxiliary and the output variables.                                 !
   !---------------------------------------------------------------------------------------!
   do y=1,ny
      do x=1,nx
         do z=1,nz
            dwdy (x,y,z) = 0.
            dvdz (x,y,z) = 0.
            aux  (x,y,z) = 0.
            curlx(x,y,z) = 0.
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the derivatives of wind components.                                          !
   !---------------------------------------------------------------------------------------!
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,wwnd,dwdy,'ydir','wpnt',topt,myxmn(:,ngrd)            &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,vwnd,dvdz,'zdir','vpnt',topt,myxmn(:,ngrd)            &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the curl in the staggered coordinates.                                       !
   !---------------------------------------------------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
            aux(x,y,z) = dwdy(x,y,z) - dvdz(x,y,z)
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Update the relative vorticity at the thermodynamic points.                         !
   !---------------------------------------------------------------------------------------!
   do y = 2,ny-1
      y1 = max(y-1,2)
      y2 = min(y,ny-1)
      do x = 2,nx-1
         do z = 1,nz
            z1           = max(z-1,1)
            z2           = min(z,nz-1)
            curlx(x,y,z) = 0.25 * ( aux(x,y1,z1) + aux(x,y1,z2)                            &
                                  + aux(x,y2,z1) + aux(x,y2,z2) )
         end do
      end do
   end do

  return
end subroutine RAMS_comp_relvortx
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine calculates the meridional component of the relative vorticity.      !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_relvorty(nx,ny,nz,curly,uwnd,wwnd,dudz,dwdx,aux,topt,ngrd)
   use somevars  , only : myxmn       & ! intent(in)
                        , myxtn       & ! intent(in)
                        , myymn       & ! intent(in)
                        , myytn       & ! intent(in)
                        , myzmn       & ! intent(in)
                        , myztn       & ! intent(in)
                        , mydeltayn   & ! intent(in)
                        , mydzmn      & ! intent(in)
                        , mydztn      & ! intent(in)
                        , mynnzp      & ! intent(in)
                        , myjdim      & ! intent(in)
                        , myihtran    ! ! intent(in)
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(inout) :: curly
   real   , dimension(nx,ny,nz), intent(inout) :: uwnd
   real   , dimension(nx,ny,nz), intent(inout) :: wwnd
   real   , dimension(nx,ny,nz), intent(inout) :: dudz
   real   , dimension(nx,ny,nz), intent(inout) :: dwdx
   real   , dimension(nx,ny,nz), intent(inout) :: aux
   real   , dimension(nx,ny)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   integer                                     :: x1
   integer                                     :: x2
   integer                                     :: z1
   integer                                     :: z2
   real                                        :: factor
   real   , dimension(nx+ny+nz)                :: dum1
   real   , dimension(nx+ny+nz)                :: dum2
   !---------------------------------------------------------------------------------------!



   !----- Bottom boundary condition for meridional wind. ----------------------------------!
   factor = myztn(1,ngrd) / myztn(2,ngrd)
   do y=1,ny
      do x=1,nx
         uwnd(x,y,1) = uwnd(x,y,2) * factor
      end do
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Initialise the auxiliary and the output variables.                                 !
   !---------------------------------------------------------------------------------------!
   do y=1,ny
      do x=1,nx
         do z=1,nz
            dwdx (x,y,z) = 0.
            dudz (x,y,z) = 0.
            aux  (x,y,z) = 0.
            curly(x,y,z) = 0.
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Find the derivatives of wind components.                                          !
   !---------------------------------------------------------------------------------------!
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,wwnd,dwdx,'xdir','wpnt',topt,myxmn(:,ngrd)            &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,uwnd,dudz,'zdir','upnt',topt,myxmn(:,ngrd)            &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !     Find the curl in the staggered coordinates.                                       !
   !---------------------------------------------------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
            aux(x,y,z) = dudz(x,y,z) - dwdx(x,y,z)
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Update the relative vorticity at the thermodynamic points.                         !
   !---------------------------------------------------------------------------------------!
   do y = 2,ny-1
      do x = 2,nx-1
         x1 = max(x-1,2)
         x2 = min(x,nx-1)
         do z = 1,nz
            z1 = max(z-1,1)
            z2 = min(z,nz-1)
            curly(x,y,z) = 0.25 * ( aux(x1,y,z1) + aux(x1,y,z2)                            &
                                  + aux(x2,y,z1) + aux(x2,y,z2) )
         end do
      end do
   end do

  return
end subroutine RAMS_comp_relvorty
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine calculates the vertical component of the relative vorticity.        !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_relvortz(nx,ny,nz,curlz,uwnd,vwnd,dudy,dvdx,aux,topt,ngrd)
   use somevars  , only : myxmn       & ! intent(in)
                        , myxtn       & ! intent(in)
                        , myymn       & ! intent(in)
                        , myytn       & ! intent(in)
                        , myzmn       & ! intent(in)
                        , myztn       & ! intent(in)
                        , mydeltayn   & ! intent(in)
                        , mydzmn      & ! intent(in)
                        , mydztn      & ! intent(in)
                        , mynnzp      & ! intent(in)
                        , myjdim      & ! intent(in)
                        , myihtran    ! ! intent(in)
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(inout) :: curlz
   real   , dimension(nx,ny,nz), intent(inout) :: uwnd
   real   , dimension(nx,ny,nz), intent(inout) :: vwnd
   real   , dimension(nx,ny,nz), intent(inout) :: dudy
   real   , dimension(nx,ny,nz), intent(inout) :: dvdx
   real   , dimension(nx,ny,nz), intent(inout) :: aux
   real   , dimension(nx,ny)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   integer                                     :: x1
   integer                                     :: x2
   integer                                     :: y1
   integer                                     :: y2
   real                                        :: factor
   real   , dimension(nx+ny+nz)                :: dum1
   real   , dimension(nx+ny+nz)                :: dum2
   !---------------------------------------------------------------------------------------!



   !----- Bottom boundary condition for meridional wind. ----------------------------------!
   factor = myztn(1,ngrd) / myztn(2,ngrd)
   do y=1,ny
      do x=1,nx
         uwnd(x,y,1) = uwnd(x,y,2) * factor
         vwnd(x,y,1) = vwnd(x,y,2) * factor
      end do
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    Initialise the auxiliary and the output variables.                                 !
   !---------------------------------------------------------------------------------------!
   do y=1,ny
      do x=1,nx
         do z=1,nz
            dvdx (x,y,z) = 0.
            dudy (x,y,z) = 0.
            aux  (x,y,z) = 0.
            curlz(x,y,z) = 0.
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Find the derivatives of wind components.                                          !
   !---------------------------------------------------------------------------------------!
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,vwnd,dvdx,'xdir','vpnt',topt,myxmn(:,ngrd)            &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,uwnd,dudy,'ydir','upnt',topt,myxmn(:,ngrd)            &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !     Find the curl in the staggered coordinates.                                       !
   !---------------------------------------------------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
            aux(x,y,z) = dvdx(x,y,z) - dudy(x,y,z)
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!







   !---------------------------------------------------------------------------------------!
   !    Update the relative vorticity at the thermodynamic points.                         !
   !---------------------------------------------------------------------------------------!
   do y = 2,ny-1
      y1 = max(y-1,2)
      y2 = min(y,ny-1)
      do x = 2,nx-1
         x1 = max(x-1,2)
         x2 = min(x,nx-1)
         do z = 1,nz
            curlz(x,y,z) = 0.25 * ( aux(x1,y1,z) + aux(x1,y2,z)                            &
                                  + aux(x2,y1,z) + aux(x2,y2,z) )
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine RAMS_comp_relvortz
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     Sub-routine for absolute vorticity.                                                  !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_totvortz(nx,ny,nz,curlz,glat)
   use rconstants, only : omega   & ! intent(in)
                        , pio180  ! ! intent(in)
   implicit none 
   !------ Arguments. ---------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(inout) :: curlz
   real   , dimension(nx,ny)   , intent(in)    :: glat
   !------ Local variables. ---------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   real                                        :: fcor
   !---------------------------------------------------------------------------------------!



   yloop: do y = 1,ny
      xloop: do x = 1,nx
         fcor = 2. * omega * sin(glat(x,y) * pio180)
         zloop: do z = 1,nz
            curlz(x,y,z) = curlz(x,y,z) + fcor
         end do zloop
      end do xloop
   end do yloop


   return
end subroutine RAMS_comp_totvortz
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_potvort(nx,ny,nz,pvort,curlx,curly,curlz,theta,dens,dthdx,dthdy,dthdz &
                            ,topt,ngrd)
   use somevars
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(inout) :: pvort
   real   , dimension(nx,ny,nz), intent(inout) :: curlx
   real   , dimension(nx,ny,nz), intent(inout) :: curly
   real   , dimension(nx,ny,nz), intent(inout) :: curlz
   real   , dimension(nx,ny,nz), intent(inout) :: theta
   real   , dimension(nx,ny,nz), intent(in)    :: dens
   real   , dimension(nx,ny,nz), intent(inout) :: dthdx
   real   , dimension(nx,ny,nz), intent(inout) :: dthdy
   real   , dimension(nx,ny,nz), intent(inout) :: dthdz
   real   , dimension(nx,ny)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   real   , dimension(nx+ny+nz)                :: dum1
   real   , dimension(nx+ny+nz)                :: dum2
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Initialise the gradients.                                                         !
   !---------------------------------------------------------------------------------------!
   do x=1,nx
      do y=1,ny
         do z=1,nz
            dthdx(x,y,z) = 0.
            dthdy(x,y,z) = 0.
            dthdz(x,y,z) = 0.
            pvort(x,y,z) = 0.
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Find the gradient components.                                                     !
   !---------------------------------------------------------------------------------------!
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,theta,dthdx,'xdir','tpnt',topt,myxmn(:,ngrd)          &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)

   call gradr(nx,ny,nz,2,nx-1,2,ny-1,theta,dthdy,'ydir','tpnt',topt,myxmn(:,ngrd)          &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)

   call gradr(nx,ny,nz,2,nx-1,2,ny-1,theta,dthdz,'zdir','tpnt',topt,myxmn(:,ngrd)          &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Find the potential vorticity.  Currently all components are considered.          !
   !---------------------------------------------------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
           pvort(x,y,z) = ( curlx(x,y,z) * dthdx(x,y,z)                                    &
                          + curly(x,y,z) * dthdy(x,y,z)                                    &
                          + curlz(x,y,z) * dthdz(x,y,z) ) / dens(x,y,z)
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine RAMS_comp_potvort
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_xydiv(nx,ny,nz,dive,uwnd,vwnd,dudx,dvdy,aux,topt,ngrd)

   use somevars
   use rconstants
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(inout) :: dive
   real   , dimension(nx,ny,nz), intent(inout) :: uwnd
   real   , dimension(nx,ny,nz), intent(inout) :: vwnd
   real   , dimension(nx,ny,nz), intent(inout) :: dudx
   real   , dimension(nx,ny,nz), intent(inout) :: dvdy
   real   , dimension(nx,ny,nz), intent(inout) :: aux
   real   , dimension(nx,ny)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   integer                                     :: x1
   integer                                     :: x2
   integer                                     :: y1
   integer                                     :: y2
   real                                        :: factor
   real   , dimension(nx+ny+nz)                :: dum1
   real   , dimension(nx+ny+nz)                :: dum2
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Scale down the wind direction to make boundary conditions.                        !
   !---------------------------------------------------------------------------------------!
   factor = myztn(1,ngrd) / myztn(2,ngrd)
   do y=1,ny
      do x=1,nx
         uwnd(x,y,1) = uwnd(x,y,2) * factor
         vwnd(x,y,1) = vwnd(x,y,2) * factor
      enddo
   enddo
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Reset divergence and wind gradients.                                               !
   !---------------------------------------------------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
            dudx(x,y,z) = 0.
            dvdy(x,y,z) = 0.
            dive(x,y,z) = 0.
            aux (x,y,z) = 0.
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the wind gradients.                                                          !
   !---------------------------------------------------------------------------------------!
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,uwnd,dudx,'xdir','upnt',topt,myxmn(:,ngrd)            &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,vwnd,dvdy,'ydir','vpnt',topt,myxmn(:,ngrd)            &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the divergence in the XY plane.                                              !
   !---------------------------------------------------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
            aux(x,y,z) = dudx(x,y,z) + dvdy(x,y,z)
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Update the divergence at the XY plane at the thermodynamic points.                 !
   !---------------------------------------------------------------------------------------!
   do y = 2,ny-1
      y1 = max(y-1,2)
      y2 = min(y,ny-1)
      do x = 2,nx-1
         x1 = max(x-1,2)
         x2 = min(x,nx-1)
         do z = 1,nz
            dive(x,y,z) = 0.25 * ( aux(x1,y1,z) + aux(x1,y2,z)                             &
                                 + aux(x2,y1,z) + aux(x2,y2,z) )
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine RAMS_comp_xydiv
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_solenoidx(nx,ny,nz,alpha,press,solex,topt,ngrd)
   use somevars
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(inout) :: alpha ! Specific Volume
   real   , dimension(nx,ny,nz), intent(inout) :: press
   real   , dimension(nx,ny,nz), intent(inout) :: solex
   real   , dimension(nx,ny)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   real   , dimension(nx,ny,nz)                :: dady
   real   , dimension(nx,ny,nz)                :: dadz
   real   , dimension(nx,ny,nz)                :: dpdy
   real   , dimension(nx,ny,nz)                :: dpdz
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   real   , dimension(nx+ny+nz)                :: dum1
   real   , dimension(nx+ny+nz)                :: dum2
   !---------------------------------------------------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
            solex(x,y,z) = 0.
         end do
      end do
   end do

   call gradr(nx,ny,nz,2,nx-1,2,ny-1,alpha,dadz,'zdir','tpnt',topt,myxmn(:,ngrd)           &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,press,dpdy,'ydir','tpnt',topt,myxmn(:,ngrd)           &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)

   do z=1,nz-1
      do y=2,ny-1
         do x=2,nx-1
            solex(x,y,z) = dadz(x,y,z) * dpdy(x,y,z)
         end do
      end do
   end do


   call gradr(nx,ny,nz,2,nx-1,2,ny-1,press,dpdz,'zdir','tpnt',topt,myxmn(:,ngrd)           &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,alpha,dady,'ydir','tpnt',topt,myxmn(:,ngrd)           &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)


   do z=1,nz-1
      do y=2,ny-1
         do x=2,nx-1
            solex(x,y,z) = solex(x,y,z) - dpdz(x,y,z) * dady(x,y,z)
         end do
      end do
   end do

  return
end subroutine RAMS_comp_solenoidx
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_solenoidy(nx,ny,nz,alpha,press,soley,topt,ngrd)
   use somevars
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(inout) :: alpha ! Specific Volume
   real   , dimension(nx,ny,nz), intent(inout) :: press
   real   , dimension(nx,ny,nz), intent(inout) :: soley
   real   , dimension(nx,ny)   , intent(in)    :: topt
   integer                     , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   real   , dimension(nx,ny,nz)                :: dadx
   real   , dimension(nx,ny,nz)                :: dadz
   real   , dimension(nx,ny,nz)                :: dpdx
   real   , dimension(nx,ny,nz)                :: dpdz
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   real   , dimension(nx+ny+nz)                :: dum1
   real   , dimension(nx+ny+nz)                :: dum2
   !----- Local constants. ----------------------------------------------------------------!
   logical                     , parameter     :: printdbg = .false.
   !---------------------------------------------------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
            soley(x,y,z) = 0.
         end do
      end do
   end do

   call gradr(nx,ny,nz,2,nx-1,2,ny-1,alpha,dadx,'xdir','tpnt',topt,myxmn(:,ngrd)           &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,press,dpdz,'zdir','tpnt',topt,myxmn(:,ngrd)           &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1)  &
              ,myjdim,myihtran)

   do z=1,nz-1
      do y=2,ny-1
         do x=2,nx-1
            soley(x,y,z) = dadx(x,y,z) * dpdz(x,y,z)
         end do
      end do
   end do

   call gradr(nx,ny,nz,2,nx-1,2,ny-1,press,dpdx,'xdir','tpnt',topt,myxmn(:,ngrd)           &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)
   call gradr(nx,ny,nz,2,nx-1,2,ny-1,alpha,dadz,'zdir','tpnt',topt,myxmn(:,ngrd)           &
             ,myxtn(:,ngrd),myymn(:,ngrd),myytn(:,ngrd),myzmn(:,ngrd),myztn(:,ngrd)        &
             ,mydeltayn(ngrd),mydzmn(:,ngrd),mydztn(:,ngrd),dum1,dum2,myzmn(mynnzp(1)-1,1) &
             ,myjdim,myihtran)

   do z=1,nz-1
      do y=2,ny-1
         do x=2,nx-1
            soley(x,y,z) = soley(x,y,z) - dadz(x,y,z) * dpdx(x,y,z)
         end do
      end do
   end do
   
   if (printdbg) then
      write (unit=61,fmt='(3(a5,1x),7(a12,1x))') '    X','    Y','    Z','       ALPHA'    &
                                                ,'       PRESS','        DADX'             &
                                                ,'        DADZ','        DPDX'             &
                                                ,'        DPDZ','       SOLEY'

      do z=1,nz-1
         do y=2,ny-1
            do x=2,nx-1
               write(unit=61,fmt='(3(i5,1x),7(es12.5,1x))') x,y,z,alpha(x,y,z)             &
                                                                 ,press(x,y,z)             &
                                                                 , dadx(x,y,z)             &
                                                                 , dadz(x,y,z)             &
                                                                 , dpdx(x,y,z)             &
                                                                 , dpdz(x,y,z)             &
                                                                 ,soley(x,y,z)
            end do
         end do
      end do
   end if

   return
end subroutine RAMS_comp_solenoidy
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_sfcwmeantemp(n1,n2,ns,np,a,b,c,d,e)
   use rconstants
   use rout_coms, only  : undefflg
   use therm_lib, only: uint2tl
   implicit none 
   integer :: n1,n2,ns,np,nlev,i,j,ip,k
   real, dimension(n1,n2,np)    :: a,d,e
   real, dimension(n1,n2,ns,np) :: b,c
   real :: temptemp,fracliq,snowarea,xmasstot
   !a  area
   !b energy
   !c  mass
   !d  nlev
   !e integrated value
   do j=1,n2
     do i=1,n1
        e(i,j,1) = 0.
        if (a(i,j,1) <= 0.99) then
           do ip=2,np
              xmasstot = 0.
              nlev=nint(d(i,j,ip))
              e(i,j,ip) = 0.
              if (nlev > 0 ) then
                 do k=1,nlev
                    if (c(i,j,k,ip) > 1.e-6) then
                       xmasstot  = xmasstot + c(i,j,k,ip)
                       call uint2tl(b(i,j,k,ip),temptemp,fracliq)
                       e(i,j,ip) = e(i,j,ip) + temptemp*c(i,j,k,ip)
                    end if
                 end do
                 if (xmasstot > 1.e-6) e(i,j,ip) = e(i,j,ip) / xmasstot - t00
              else
                 e(i,j,ip) = 0.
              end if
           end do

           snowarea = 0.
           do ip=2,np
              if( nint(d(i,j,ip)) > 0) then
                 snowarea = snowarea + a(i,j,ip)
                 e(i,j,1) = e(i,j,1) + e(i,j,ip)*a(i,j,ip)
              end if
           end do
           if (snowarea > 0.) then
              e(i,j,1) = e(i,j,1) / snowarea
           else
              e(i,j,1) = undefflg
           end if
        else
           e(i,j,1) = undefflg
        end if
     end do
   end do
   return
end subroutine RAMS_comp_sfcwmeantemp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_thetaeiv(n1,n2,n3,xxx,temp,pres,rv,rtp)
   use therm_lib, only: thetaeiv ! ! function
   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   integer                     , intent(in)    :: n1
   integer                     , intent(in)    :: n2
   integer                     , intent(in)    :: n3
   real   , dimension(n1,n2,n3), intent(inout) :: xxx
   real   , dimension(n1,n2,n3), intent(inout) :: temp
   real   , dimension(n1,n2,n3), intent(inout) :: pres
   real   , dimension(n1,n2,n3), intent(inout) :: rv
   real   , dimension(n1,n2,n3), intent(inout) :: rtp
   !------ Local variables. ---------------------------------------------------------------!
   integer                                     :: i
   integer                                     :: j
   integer                                     :: k
   !---------------------------------------------------------------------------------------!


   !----- xxx comes as thil, gets out as theta_e_iv. --------------------------------------!
   do k=1,n3
      do j=1,n2
         do i=1,n1
            if (rtp(i,j,k) < rv(i,j,k)) rtp(i,j,k) = rv(i,j,k)
            xxx(i,j,k)=thetaeiv(xxx(i,j,k),pres(i,j,k),temp(i,j,k),rv(i,j,k),rtp(i,j,k)    &
                               ,.true.)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_thetaeiv
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_thetaeivs(n1,n2,n3,xxx,temp,pres,rliq,rice)
   use therm_lib, only : thetaeivs  & ! function
                       , rslif      ! ! function
   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   integer                     , intent(in)    :: n1
   integer                     , intent(in)    :: n2
   integer                     , intent(in)    :: n3
   real   , dimension(n1,n2,n3), intent(inout) :: xxx
   real   , dimension(n1,n2,n3), intent(inout) :: temp
   real   , dimension(n1,n2,n3), intent(inout) :: pres
   real   , dimension(n1,n2,n3), intent(inout) :: rliq
   real   , dimension(n1,n2,n3), intent(inout) :: rice
   !------ Local variables. ---------------------------------------------------------------!
   integer                                     :: i
   integer                                     :: j
   integer                                     :: k
   real                                        :: rsat
   !---------------------------------------------------------------------------------------!


   !----- xxx comes as thil, gets out as theta_e_ivs. -------------------------------------!
   do k=1,n3
      do j=1,n2
         do i=1,n1
            rsat       = rslif(pres(i,j,k),temp(i,j,k))
            xxx(i,j,k) = thetaeivs(xxx(i,j,k),temp(i,j,k),rsat,rliq(i,j,k),rice(i,j,k))
         end do
      end do
   end do
   return
end subroutine RAMS_comp_thetaeivs
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_sfcwinteg(n1,n2,ns,np,a,c,d,e)
   use rout_coms, only : undefflg
   implicit none 
   integer :: n1,n2,ns,np,nlev,i,j,ip,k
   real, dimension(n1,n2,np)    :: a,d,e
   real, dimension(n1,n2,ns,np) :: c
   !a  area                                                   
   !c mass/depth                                             
   !d  nlev                                                   
   !e  integrated value                                       
      do j=1,n2                                               
        do i=1,n1                                             
           if (a(i,j,1) > 0.99) then                          
              e(i,j,1) = undefflg                           
           else
              e(i,j,1) = 0.                                  
              do ip=2,np                                     
                 nlev=nint(d(i,j,ip))
                 e(i,j,ip) = 0.                               
                 do k=1,nlev                                  
                    e(i,j,ip) = e(i,j,ip) + c(i,j,k,ip)      
                 end do                                       
                 e(i,j,1) = e(i,j,1) + e(i,j,ip)*a(i,j,ip)    
              end do                                          
              e(i,j,1) = e(i,j,1) / (1.-a(i,j,1))            
           end if                                             
        end do                                                
      end do                                                  
   return        
end subroutine RAMS_comp_sfcwinteg
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This routine is for quantities defined for all patches.                               !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_patchsum(nx,ny,nz,np,iovar)
   use leaf_coms, only : min_patch_area
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                         , intent(in)    :: nx
   integer                         , intent(in)    :: ny
   integer                         , intent(in)    :: nz
   integer                         , intent(in)    :: np
   real   , dimension(*)           , intent(inout) :: iovar
   !----- Local variables. ----------------------------------------------------------------!
   integer                                         :: x
   integer                                         :: y
   integer                                         :: z
   integer                                         :: p
   integer                                         :: iareaa
   integer                                         :: iareaz
   integer                                         :: ivala
   integer                                         :: ivalz
   integer                                         :: iouta
   integer                                         :: ioutz
   integer                                         :: narea
   integer                                         :: nvals
   integer                                         :: nout
   real                                            :: totarea
   real   , dimension(nx,ny,nz,np)                 :: patval
   real   , dimension(nx,ny,nz)                    :: psum
   real   , dimension(nx,ny,np)                    :: pfarea
   !---------------------------------------------------------------------------------------!

   !----- Define the indices for copying in and out. --------------------------------------!
   narea  = nx * ny * np
   nvals  = nx * ny * nz * np
   nout   = nx * ny * nz
   iareaa = 1
   iareaz = narea
   ivala  = iareaz + 1
   ivalz  = iareaz + nvals
   iouta  = 1
   ioutz  = nout

   !----- Copy the patch fraction area to a scratch array. --------------------------------!
   call atob(narea,iovar(iareaa:iareaz),pfarea)
   
   !----- Copy the patch-structured data to a scratch array. ------------------------------!
   call atob(nvals,iovar(ivala:ivalz),patval)

   !----- Compute the patch weigthed average, including water. ----------------------------!
   do z = 1,nz
      do y = 1,ny
         do x = 1,nx
            psum(x,y,z) = 0.
            totarea     = 0.
            ploop: do p = 1,np
               if (pfarea(x,y,p) <= min_patch_area) cycle ploop
               
               psum(x,y,z) = psum(x,y,z) + pfarea(x,y,p) * patval(x,y,z,p)
               totarea     = totarea + pfarea(x,y,p)
            end do ploop
            psum(x,y,z) = psum(x,y,z) / totarea
         end do
      end do
   end do

   !----- Copy psum into iovar, which will be used for output. ----------------------------!
   call atob(nout,psum,iovar(iouta:ioutz))

   return
end subroutine RAMS_comp_patchsum
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This routine is for quantities that are not defined for water patches.              !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_patchsum_l(nx,ny,nz,np,iovar)
   use leaf_coms, only : min_patch_area
   use rout_coms, only : undefflg
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                         , intent(in)    :: nx
   integer                         , intent(in)    :: ny
   integer                         , intent(in)    :: nz
   integer                         , intent(in)    :: np
   real   , dimension(*)           , intent(inout) :: iovar
   !----- Local variables. ----------------------------------------------------------------!
   integer                                         :: x
   integer                                         :: y
   integer                                         :: z
   integer                                         :: p
   integer                                         :: iareaa
   integer                                         :: iareaz
   integer                                         :: ivala
   integer                                         :: ivalz
   integer                                         :: iouta
   integer                                         :: ioutz
   integer                                         :: narea
   integer                                         :: nvals
   integer                                         :: nout
   real                                            :: landarea
   real   , dimension(nx,ny,nz,np)                 :: patval
   real   , dimension(nx,ny,nz)                    :: psum
   real   , dimension(nx,ny,np)                    :: pfarea
   !---------------------------------------------------------------------------------------!


   !----- Define the indices for copying in and out. --------------------------------------!
   narea  = nx * ny * np
   nvals  = nx * ny * nz * np
   nout   = nx * ny * nz
   iareaa = 1
   iareaz = narea
   ivala  = iareaz + 1
   ivalz  = iareaz + nvals
   iouta  = 1
   ioutz  = nout

   !----- Copy the patch fraction area to a scratch array. --------------------------------!
   call atob(narea,iovar(iareaa:iareaz),pfarea)
   
   !----- Copy the patch-structured data to a scratch array. ------------------------------!
   call atob(nvals,iovar(ivala:ivalz),patval)

   !----- Compute the patch weigthed average, excluding water. ----------------------------!
   do z = 1,nz
      do y = 1,ny
         do x = 1,nx
            if (pfarea(x,y,1) < 1.-min_patch_area) then
               psum(x,y,z) = 0.
               landarea    = 0.
               do p = 2,np
                  psum(x,y,z) = psum(x,y,z) + pfarea(x,y,p) * patval(x,y,z,p)
                  landarea    = landarea    + pfarea(x,y,p)
               end do
               psum(x,y,z)    = psum(x,y,z) / landarea
            else
               psum(x,y,z)    = undefflg
            end if
         end do
      end do
   end do

   !----- Copy psum into iovar, which will be used for output. ----------------------------!
   call atob(nout,psum,iovar(iouta:ioutz))

   return
end subroutine RAMS_comp_patchsum_l
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      Extract value from largest patch.                                                   !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_bigpatch(n1,n2,n3,n4,a,f,bpat)
   use somevars, only : mynbig
   implicit none
   integer, intent(in) :: n1,n2,n3,n4
   real   , dimension(n1,n2,n3,n4) , intent(in)    :: a
   real   , dimension(n1,n2,mynbig), intent(inout) :: f
   real   , dimension(n1,n2,n3)    , intent(out)   :: bpat
   integer                                         :: i,j,k,ip
   do k = 1,n3
      do j = 1,n2
         do i = 1,n1
            if (f(i,j,2) >= f(i,j,1)) then
               bpat(i,j,k) = a(i,j,k,2)
            else
               bpat(i,j,k) = a(i,j,k,1)
            end if
         end do
      end do
   end do

   !---------------------------------------------------------------------------------------!
   !      Copy bpat into f, which was passed in as a(1).                                   !
   !---------------------------------------------------------------------------------------!

   do k = 1,n3
      do j = 1,n2
         do i = 1,n1
            f(i,j,k) = bpat(i,j,k)
         end do
      end do
   end do
   return
end subroutine RAMS_comp_bigpatch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routines finds the vegetation temperature.                                  !
!                                                                                          !
! Input:  a = Veg. Energy        in     J/m2                                               !
!         b = Veg. Water         in    kg/m2                                               !
!         c = Veg. Heat capacity in   J/m2/K                                               !
!         e = Canopy Theta       in  Celsius                                               !
! Output: a = Temperature in Celsius                                                       !
!         b = Fraction in liquid phase                                                     !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_tvegc(n1,n2,n3,a,b,c,e)
   use rconstants, only : t00       ! ! intent(in)
   use therm_lib , only : uextcm2tl ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: n1
   integer                          , intent(in)    :: n2
   integer                          , intent(in)    :: n3
   real(kind=4), dimension(n1,n2,n3), intent(in)    :: c
   real(kind=4), dimension(n1,n2,n3), intent(in)    :: e
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: a
   real(kind=4), dimension(n1,n2,n3), intent(inout) :: b
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: i
   integer                                          :: j
   integer                                          :: k
   real(kind=4)                                     :: temptemp
   real(kind=4)                                     :: fracliq
   !---------------------------------------------------------------------------------------!



   do k=1,n3
      do j=1,n2
         do i=1,n1
            !------------------------------------------------------------------------------!
            !   Compute tveg only if there is enough heat capacity, otherwise assign       !
            ! canopy temperature.                                                          !
            !------------------------------------------------------------------------------!
            if (c(i,j,k) > 10.) then
               call uextcm2tl(a(i,j,k),b(i,j,k),c(i,j,k),temptemp,fracliq)
               a(i,j,k) = temptemp-t00
               b(i,j,k) = fracliq
            else
               a(i,j,k) = e(i,j,k)
               b(i,j,k) = 0.5
            end if
            !------------------------------------------------------------------------------!
         end do
      end do
   end do
   return
end subroutine RAMS_comp_tvegc
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will compute the mixing ratio (water vapour or CO2 or any other      !
! tracer), using the similarity theory.  It currently uses either the Louis (1979), the    !
! Oncley and Dudhia (1995) or the Beljaars-Holtslag (1991) methods, and the choice is done !
! based on the ISTAR used in the run.                                                      !
!                                                                                          !
! 1. Based on L79;                                                                         !
! 2. Based on: OD95, but with some terms computed as in L79 and B71 to avoid singular-     !
!    ities (now using the iterative method to find zeta).                                  !
! 3. Based on BH91, using an iterative method to find zeta, and using the modified         !
!    equation for stable layers.                                                           !
! 4. Based on CLM04, with special functions for very stable and very stable case, even     !
!    though we use a different functional form for very unstable case for momentum.        !
!    This is ensure that phi_m decreases monotonically as zeta becomes more negative.      !
!    We use a power law of order of -1/6 instead.                                          !
!                                                                                          !
! References:                                                                              !
! B71.  BUSINGER, J.A, et. al; Flux-Profile relationships in the atmospheric surface       !
!           layer. J. Atmos. Sci., 28, 181-189, 1971.                                      !
! L79.  LOUIS, J.F.; Parametric Model of vertical eddy fluxes in the atmosphere.           !
!           Boundary-Layer Meteor., 17, 187-202, 1979.                                     !
! BH91. BELJAARS, A.C.M.; HOLTSLAG, A.A.M.; Flux parameterization over land surfaces for   !
!           atmospheric models. J. Appl. Meteor., 30, 327-341, 1991.                       !
! OD95. ONCLEY, S.P.; DUDHIA, J.; Evaluation of surface fluxes from MM5 using observa-     !
!           tions.  Mon. Wea. Rev., 123, 3344-3357, 1995.                                  !
! CLM04. OLESON, K. W., et al.; Technical description of the community land model (CLM)    !
!           NCAR Technical Note NCAR/TN-461+STR, Boulder, CO, May 2004.                    !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine RAMS_reduced_prop(nx,ny,nz,np,ng,which,topt,theta_atm,rvap_atm,co2_atm,uspd_atm &
                            ,theta_can,rvap_can,co2_can,prss_can,zout,rough,rib,zeta,parea &
                            ,ustar,tstar,rstar,cstar,varred)
   use rpost_coms, only : isfcl          ! ! intent(in)
   use somevars  , only : myztn          & ! intent(in)
                        , myzmn          & ! intent(in)
                        , mynnzp         & ! intent(in)
                        , myistar        ! ! intent(in)
   use rconstants, only : grav           & ! intent(in)
                        , p00i           & ! intent(in)
                        , rocp           & ! intent(in)
                        , vonk           & ! intent(in)
                        , ep             & ! intent(in)
                        , toodry         ! ! intent(in)
   use therm_lib , only : virtt          & ! function
                        , eslif          & ! function
                        , rslif          & ! function
                        , tslif          & ! function
                        , thetaeiv       ! ! function
   use leaf_coms , only : ustmin         & ! intent(in)
                        , ubmin          & ! intent(in)
                        , bl79           & ! intent(in)
                        , csm            & ! intent(in)
                        , csh            & ! intent(in)
                        , dl79           & ! intent(in)
                        , ribmax         & ! intent(in)
                        , tprandtl       & ! intent(in)
                        , z0moz0h        & ! intent(in)
                        , z0hoz0m        & ! intent(in)
                        , min_patch_area & ! intent(in)
                        , psim           & ! function
                        , psih           & ! function
                        , zoobukhov      ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in)    :: nx
   integer                  , intent(in)    :: ny
   integer                  , intent(in)    :: nz
   integer                  , intent(in)    :: np
   integer                  , intent(in)    :: ng
   character(len=4)         , intent(in)    :: which
   real, dimension(nx,ny,nz), intent(in)    :: theta_atm
   real, dimension(nx,ny,nz), intent(in)    :: rvap_atm
   real, dimension(nx,ny,nz), intent(in)    :: co2_atm
   real, dimension(nx,ny,nz), intent(in)    :: uspd_atm
   real, dimension(nx,ny,np), intent(in)    :: theta_can
   real, dimension(nx,ny,np), intent(in)    :: rvap_can
   real, dimension(nx,ny,np), intent(in)    :: co2_can
   real, dimension(nx,ny,np), intent(in)    :: prss_can
   real, dimension(nx,ny)   , intent(in)    :: topt
   real                     , intent(in)    :: zout
   real, dimension(nx,ny,np), intent(in)    :: rough
   real, dimension(nx,ny,np), intent(inout) :: rib
   real, dimension(nx,ny,np), intent(inout) :: zeta
   real, dimension(nx,ny,np), intent(in)    :: parea
   real, dimension(nx,ny,np), intent(inout) :: ustar
   real, dimension(nx,ny,np), intent(inout) :: tstar
   real, dimension(nx,ny,np), intent(inout) :: rstar
   real, dimension(nx,ny,np), intent(inout) :: cstar
   real, dimension(nx,ny)   , intent(inout) :: varred
   !----- Local variables. ----------------------------------------------------------------!
   integer           :: x            ! Longitude counter
   integer           :: y            ! Latitude counter
   integer           :: p            ! Patch counter
   real              :: zgrd         ! Grid bottom
   real              :: ztop         ! Grid top
   real              :: zref         ! Reference height
   real              :: rtgt         ! Terrain-following coordinate correction factor. 
   real              :: thetav_atm   ! Atmospheric virtual potential temperature
   real              :: thetav_can   ! Canopy air space virtual potential temperature
   logical           :: stable       ! Stable state
   logical           :: is_ed2       ! This is an ED-2 run
   real              :: zroz0m       ! zref/rough(momentum)
   real              :: lnzroz0m     ! ln[zref/rough(momentum)]
   real              :: zroz0h       ! zref/rough(heat)
   real              :: lnzroz0h     ! ln[zref/rough(heat)]
   real              :: zooz0m       ! zout/rough(momentum)
   real              :: lnzooz0m     ! ln[zout/rough(momentum)]
   real              :: zooz0h       ! zout/rough(heat)
   real              :: lnzooz0h     ! ln[zout/rough(heat)]
   real              :: uref         ! Reference wind speed.
   real              :: ured         ! Wind reduced to the level of interest.
   real              :: rvapp        ! Output variable for this patch.
   real              :: thetap        ! Output variable for this patch.
   real              :: redp         ! Output variable for this patch.
   real              :: validarea    ! Total area where we have results.
   !----- Local variables, used by L79. ---------------------------------------------------!
   real              :: a2r          ! Drag coefficient in neutral conditions
   real              :: a2o          ! Drag coefficient in neutral conditions
   real              :: fhr          ! Stability parameter for heat at z = zref
   real              :: fmr          ! Stability parameter for momentum at z = zref
   real              :: fho          ! Stability parameter for heat at z = zout
   real              :: fmo          ! Stability parameter for momentum at z = zout
   real              :: c2           ! Part of the c coefficient common to momentum & heat.
   real              :: c3           ! Another auxiliary variable.
   real              :: multh        ! Factor to be multiplied to get the heat/water.
   real              :: cm           ! c coefficient times |Rib|^1/2 for momentum.
   real              :: ch           ! c coefficient times |Rib|^1/2 for heat.
   real              :: ee           ! (z/z0)^1/3 -1. for eqn. 20 w/o assuming z/z0 >> 1.
   !----- Local variables, used by OD95 and/or BH91. --------------------------------------!
   real              :: zetaom       ! (zout + roughness(momentum))/(Obukhov length).
   real              :: zetaoh       ! (zout + roughness(heat)    )/(Obukhov length).
   real              :: zeta0m       ! roughness(momentum)/(Obukhov length).
   real              :: zeta0h       ! roughness(heat)/(Obukhov length).
   real              :: ribold       ! Bulk richardson number.
   real              :: tdmax        ! Maximum dew point temperature.
   real              :: rvmax        ! Maximum mixing ratio.
   !----- External functions. -------------------------------------------------------------!
   real, external    :: cbrt         ! Cubic root
   !---------------------------------------------------------------------------------------!

   !----- Decide whether this is an ED-2 run or not. --------------------------------------!
   is_ed2 = isfcl == 5


   !----- Define grid bottom and top. -----------------------------------------------------!
   zgrd = myztn(2,ng)
   ztop = myzmn(mynnzp(1)-1,1)

   yloop: do y = 1,ny
      xloop: do x = 1,nx
         rtgt = 1. - topt(x,y) / ztop
         zref = zgrd * rtgt
         
         !----- Compute the virtual potential temperature at the model first level. -------!
         thetav_atm = virtt(theta_atm(x,y,2),rvap_atm(x,y,2),rvap_atm(x,y,2))

         !----- Initialise the output variable. -------------------------------------------!
         varred(x,y) = 0.
         validarea   = 0.

         ploop: do p = 1,np
            !----- Skip patch if the area is tiny. ----------------------------------------!
            if (parea(x,y,p) <= min_patch_area) cycle ploop

            !----- Compute the virtual pot. temperature at the canopy air space (CAS). ----!
            thetav_can = virtt(theta_can(x,y,p),rvap_can(x,y,p),rvap_can(x,y,p))

            !----- Compute the reference wind speed. --------------------------------------!
            uref = max(uspd_atm(x,y,2),ubmin)

            !----- Re-compute the Richardson number if this is an ED-2 run. ---------------!
            if (is_ed2) then
               rib(x,y,p) = 2.0 * grav * (zref-rough(x,y,p)) * (thetav_atm-thetav_can)     &
                          / ( (thetav_atm+thetav_can) * uref * uref)
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Find the bulk Richardson number and determine whether the layer is       !
            ! stable or not.                                                               !
            !------------------------------------------------------------------------------!
            stable     = rib(x,y,p) > 0.0
            !------------------------------------------------------------------------------!


            !------ Check whether we must correct the bulk Richardson number and stars. ---!
            if (rib(x,y,p) > ribmax .and. myistar /= 1 .and. is_ed2) then
               uref = sqrt(rib(x,y,p) / ribmax) * uref
               rib(x,y,p) = ribmax
            end if


            !------------------------------------------------------------------------------!
            !     Find some variables common to all methods.  Notice that, unlike          !
            ! leaf_stars, we here use the output height, not the reference height.         !
            !------------------------------------------------------------------------------!
            zroz0m      = (zref)/rough(x,y,p)
            lnzroz0m    = log(zroz0m)
            zroz0h      = z0moz0h * zroz0m
            lnzroz0h    = log(zroz0h)
            zooz0m      = (zout+rough(x,y,p))/rough(x,y,p)
            lnzooz0m    = log(zooz0m)
            zooz0h      = z0moz0h * zooz0m
            lnzooz0h    = log(zooz0h)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Here we find the standard functions of heat and momentum to integrate    !
            ! the result.                                                                  !
            !------------------------------------------------------------------------------!
            select case (myistar)
            case (1)
               !---------------------------------------------------------------------------!
               !     Here we will use L79 model, the BRAMS default.                        !
               !---------------------------------------------------------------------------!

               !----- Compute the a-square factor and the coefficient to find theta*. -----!
               a2r  = (vonk / lnzroz0m) ** 2.
               a2o  = (vonk / lnzooz0m) ** 2.
               !---------------------------------------------------------------------------!

               if (stable) then
                  !------------------------------------------------------------------------!
                  !     Stable case.                                                       !
                  !------------------------------------------------------------------------!
                  fmr = 1.0 / (1.0 + (2.0*bl79 * rib(x,y,p) / sqrt(1.0 + dl79*rib(x,y,p))))
                  fhr = 1.0 / (1.0 + (3.0*bl79 * rib(x,y,p) * sqrt(1.0 + dl79*rib(x,y,p))))
                  fmo = fmr
                  fho = fhr
                  !------------------------------------------------------------------------!

               else
                  !------------------------------------------------------------------------!
                  !     Unstable case.  The only difference from the original method is    !
                  ! that we no longer assume z >> z0, so the "c" coefficient uses the full !
                  ! z/z0 term.                                                             !
                  !------------------------------------------------------------------------!
                  ee  = cbrt(zroz0m) - 1.
                  c2  = bl79 * a2r * ee * sqrt(ee * abs(rib(x,y,p)))
                  cm  = csm * c2
                  ch  = csh * c2
                  fmr = (1.0 - 2.0 * bl79 * rib(x,y,p) / (1.0 + 2.0 * cm))
                  fhr = (1.0 - 3.0 * bl79 * rib(x,y,p) / (1.0 + 3.0 * ch))
                  ee  = cbrt(zooz0m) - 1.
                  c2  = bl79 * a2o * ee * sqrt(ee * abs(rib(x,y,p)))
                  cm  = csm * c2
                  ch  = csh * c2
                  fmo = (1.0 - 2.0 * bl79 * rib(x,y,p) / (1.0 + 2.0 * cm))
                  fho = (1.0 - 3.0 * bl79 * rib(x,y,p) / (1.0 + 3.0 * ch))
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!

               if (is_ed2) then
                  !----- Re-compute the stars if this is an ED-2 run. ---------------------!
                  ustar(x,y,p) = max(ustmin,uref * sqrt(a2r * fmr))
                  !----- Finding the coefficient to scale the other stars. ----------------!
                  c3 = a2r * uref * fhr / ustar(x,y,p)
                  !----- Computing the other scales. --------------------------------------!
                  rstar(x,y,p) = c3 * (rvap_atm (x,y,2) - rvap_can (x,y,p)   )
                  tstar(x,y,p) = c3 * (theta_atm(x,y,2) - theta_can(x,y,p)   )
                  cstar(x,y,p) = c3 * (co2_atm  (x,y,2) - co2_can  (x,y,p)   )

                  !----- Compute zeta from u* and T* --------------------------------------!
                  zeta(x,y,p)  = grav * vonk * c3 * (thetav_atm - thetav_can)              &
                               / (thetav_atm * ustar(x,y,p) * ustar(x,y,p))
               end if

               ured  = max(ubmin, ustar(x,y,p) * lnzooz0m / (vonk * sqrt(fmo)))
               multh = tprandtl * ustar(x,y,p) * lnzooz0m / (vonk * ured * fho)

            case default
               !---------------------------------------------------------------------------!
               ! 2. Here we use the model proposed by OD95, the standard for MM5, but with !
               !    some terms that were computed in B71 (namely, the "0" terms), which    !
               !    prevent singularities.                                                 !
               !    However we know zeta, so zeta0 can be written as z0/z * zeta.          !
               ! 3. Here we use the model proposed by BH91, which is almost the same as    !
               !    the OD95 method, except that the stable functions are computed in a    !
               !    more generic way.  BH91 claim that the oft-used approximation          !
               !    (-beta*zeta) can cause poor ventilation of the stable layer, leading   !
               !    to decoupling between the atmosphere and the canopy air space and      !
               !    excessive cooling.                                                     !
               ! 4. Here we use a similar approach as in CLM04, excepth that the momentum  !
               !    flux gradient function for the unstable case for momentum is switched  !
               !    by a power of -1/6 (kind of the square of the heat one).  This is to   !
               !    guarantee that the psi function doesn't have local maxima/minima.      !
               !---------------------------------------------------------------------------!
               if (is_ed2) then 
                  !----- Make sure that the bulk Richardson number is not above ribmax. ---!
                  zeta(x,y,p) = zoobukhov(rib(x,y,p),zref,rough(x,y,p),zroz0m,lnzroz0m     &
                                         ,zroz0h,lnzroz0h,stable,myistar)
               end if

               zetaom = (zout + rough(x,y,p)) * zeta(x,y,p) / zref
               zetaoh = zetaom
               zeta0m = rough(x,y,p) * zeta(x,y,p) / zref
               zeta0h = zeta0m
               !---------------------------------------------------------------------------!


               !----- Re-compute the stars if this is an ED-2 run. ------------------------!
               if (is_ed2) then
                  ustar(x,y,p) = max (ustmin, vonk * uref                                  &
                                            / (lnzroz0m - psim(zeta(x,y,p),stable,myistar) &
                                                        + psim(zeta0m,stable,myistar)    ))

                  !----- Finding the coefficient to scale the other stars. ----------------!
                  c3    = vonk / (tprandtl * (lnzroz0m - psih(zeta(x,y,p),stable,myistar)  &
                                                       + psih(zeta0m,stable,myistar)     ))
                  !----- Computing the other scales. --------------------------------------!
                  rstar(x,y,p) = c3 * (rvap_atm (x,y,2) - rvap_can (x,y,p)   )
                  tstar(x,y,p) = c3 * (theta_atm(x,y,2) - theta_can(x,y,p)   )
                  cstar(x,y,p) = c3 * (co2_atm  (x,y,2) - co2_can  (x,y,p)   )
               end if
               !---------------------------------------------------------------------------!


               ured  = ustar(x,y,p) * ( lnzooz0m - psim(zetaom,stable,myistar)             &
                                      + psim(zeta0m,stable,myistar) ) / vonk
               multh = tprandtl     * ( lnzooz0m - psih(zetaoh,stable,myistar)             &
                                      + psih(zeta0h,stable,myistar) ) / vonk

            end select
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !      Find the temperature and saturation mixing ratio.  Both will be used to !
            ! constrain humidity under 100%.                                               !
            !------------------------------------------------------------------------------!
            tdmax = (theta_can(x,y,p) + tstar(x,y,p) * multh)                              &
                  * (p00i * prss_can(x,y,p)) ** rocp
            rvmax = rslif(prss_can(x,y,p),tdmax)
            !------------------------------------------------------------------------------!


            !----- We now compute the reference variable for this patch. ------------------!
            select case (which)
            case('WIND')
               redp = ured
            case('THET')
               redp = theta_can(x,y,p) + tstar(x,y,p) * multh
            case('TEMP')
               !----- Tdmax is the reduced temperature above. -----------------------------!
               redp = tdmax
            case('RVAP')
               redp = max(toodry,min(rvmax,rvap_can(x,y,p)  + rstar(x,y,p) * multh))
            case('THEE')
               !----- Find the potential temperature. -------------------------------------!
               thetap = theta_can(x,y,p) + tstar(x,y,p) * multh
               rvapp  = max(toodry,min(rvmax,rvap_can(x,y,p)  + rstar(x,y,p) * multh))
               !----- Convert it to temperature. ------------------------------------------!
               redp   = thetap * (p00i * prss_can(x,y,p)) ** rocp
               !----- Find thetae_iv. -----------------------------------------------------!
               redp   = thetaeiv(thetap,prss_can(x,y,p),redp,rvapp,rvapp)
            case('CO_2')
               redp = max(toodry,co2_can(x,y,p)   + cstar(x,y,p) * multh)
            case('TDEW')
               !----- Find the mixing ratio. ----------------------------------------------!
               redp = max(toodry,min(rvmax,rvap_can(x,y,p)  + rstar(x,y,p) * multh))
               !----- Convert it to vapour partial pressure. ------------------------------!
               redp = prss_can(x,y,p) * redp / (ep + redp)
               !----- Find the dew/frost point. -------------------------------------------!
               redp = tslif(redp)
            case('ZETA')
               redp = zeta(x,y,p)
            case('RICH')
               redp = rib(x,y,p)
            case('USTR')
               redp = ustar(x,y,p)
            case('TSTR')
               redp = tstar(x,y,p)
            case('RSTR')
               redp = rstar(x,y,p)
            case('CSTR')
               redp = cstar(x,y,p)
            end select

            validarea   = validarea   + parea(x,y,p)
            varred(x,y) = varred(x,y) + parea(x,y,p) * redp
         end do ploop

         varred(x,y) = varred(x,y) / validarea
      end do xloop
   end do yloop

   return
end subroutine RAMS_reduced_prop
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the relative humidity.                                      !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_relhum(nx,ny,no,tdinrhout,temp)
   use therm_lib, only : eslif
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: no
   real, dimension(nx,ny,no)   , intent(inout) :: tdinrhout
   real, dimension(nx,ny,no)   , intent(in)    :: temp
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: o
   real                                        :: e
   real                                        :: es
   !---------------------------------------------------------------------------------------!

   oloop: do o=1,no
      yloop: do y=1,ny
         xloop: do x=1,nx
            e  = eslif(tdinrhout(x,y,o))
            es = eslif(temp(x,y,o))
            tdinrhout(x,y,o) = max(0.,min(1.,e / es))
         end do xloop
      end do yloop
   end do oloop
   return
end subroutine RAMS_comp_relhum
!==========================================================================================!
!==========================================================================================!




!==========================================================================================!
!==========================================================================================!
!     This subroutine avoids using tiny numbers to some variables, by flushing them to     !
! zero if they are too small.  This avoids FPE, especially when patch integration is about !
! to happen.                                                                               !
!------------------------------------------------------------------------------------------!
subroutine RAMS_flush_to_zero(nx,ny,nz,np,myvar,threshold)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: np
   real, dimension(nx,ny,nz,np), intent(inout) :: myvar
   real                        , intent(in)    :: threshold
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   integer                                     :: p
   !---------------------------------------------------------------------------------------!
   ploop: do p=1,np
      zloop: do z=1,nz
         yloop: do y=1,ny
            xloop: do x=1,nx
               if (abs(myvar(x,y,z,p)) < threshold) myvar(x,y,z,p) = 0.
            end do xloop
         end do yloop
      end do zloop
   end do ploop
   return
end subroutine RAMS_flush_to_zero
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!       Will Cheng's code for calculating slp with mm5's GRAPH method.  Added for          !
!  calculating SLP from MM5 algorithm.                                                     !
!    The subroutine calculates SLP from an algorithm taken from  GRAPH, a post-processing  !
! package of MM5 V3.3                                                                      !
!                                                                                          !
!    Input: theta - potential temperature (K)         3D                                   !
!           pp    - Exner function        (J/kg K)    3D                                   !
!           z     - terrain               (m)         2D                                   !
!                                                                                          !
!    Ouput: SLP   - sea-level pressure    (hPa)       2D                                   !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_slpmm5(n1,n2,n3,theta,pp,z,slp)
   use rconstants, only : cpdryi      & ! intent(in)
                        , rdry        & ! intent(in)
                        , cpor        & ! intent(in)
                        , p00         ! ! intent(in)
   use therm_lib , only : exner2press ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                    , intent(in)    :: n1
   integer                    , intent(in)    :: n2
   integer                    , intent(in)    :: n3
   real, dimension(n1,n2,n3)  , intent(in)    :: theta
   real, dimension(n1,n2,n3)  , intent(in)    :: pp
   real, dimension(n1,n2)     , intent(in)    :: z
   real, dimension(n1,n2)     , intent(inout) :: slp
   !----- Local variables. ----------------------------------------------------------------!
   integer                                    :: i
   integer                                    :: j
   integer                                    :: k
   integer                                    :: kk
   real, dimension(n1,n2)                     :: sfp
   real, dimension(n1,n2)                     :: ts
   real, dimension(n1,n2,n3-1)                :: t_mm5
   real, dimension(n1,n2,n3-1)                :: p_mm5
   !---------------------------------------------------------------------------------------!
   
   do j = 1,n2
      do i = 1,n1
         !----- Calculate surface pressure. -----------------------------------------------!
         sfp(i,j) = exner2press(0.5*(pp(i,j,1)+pp(i,j,2))) *.01
         !----- Calculate surface temp. ---------------------------------------------------!
         ts(i,j)  = 0.5 * cpdryi * (theta(i,j,1)*pp(i,j,1) + theta(i,j,2)*pp(j,j,2))
      end do
   end do

   do k = 2,n3
      kk = n3-k+1
      do j = 1,n2
         do i = 1,n1
            !----- Flip arrays upside down for input to GRAPH subroutine. -----------------!
            t_mm5(i,j,kk) = theta(i,j,k) * pp(i,j,k) * cpdryi
            p_mm5(i,j,kk) = exner2press(pp(i,j,k)) * .01
         end do
      end do
   end do

   call seaprs_0(t_mm5,p_mm5,z,sfp,ts,n1,n2,n3-1,slp)

   return
end subroutine RAMS_comp_slpmm5
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes sea level pressure from the rule                            !
!              t1/t2=(p1/p2)**(gamma*r/g).                                                 !
!                                                                                          !
!     *** Levels go from top-down ***                                                      !
!                                                                                          !
!     Input       t        temperature (Kelvin)                3D                          !
!                 ter      terrain     (m)                     2D                          !
!                 sfp      surface pressure (hPa)              2D                          !
!                 imx      dot point dimension n-s                                         !
!                 jmx      dot point dimension e-w                                         !
!                 kx       number of vertical levels                                       !
!                                                                                          !
!     Output      slp      sea level pressure (hPa)            2D                          !
!                                                                                          !
!------------------------------------------------------------------------------------------!
subroutine seaprs_0(t,pp,ter,sfp,ts,imx,jmx,kx,slp)
   use rconstants, only : rdry & ! intent(in)
                        , grav & ! intent(in)
                        , t00  ! ! intent(in)
   implicit none
   !----- Local constants. ----------------------------------------------------------------!
   real                          , parameter     :: gamma  = 6.5e-3
   real                          , parameter     :: tcrit  = t00+17.5
   real                          , parameter     :: pconst = 100.
   real                          , parameter     :: xterm  = gamma * rdry / grav
   !----- Arguments. ----------------------------------------------------------------------!
   integer                       , intent(in)    :: imx
   integer                       , intent(in)    :: jmx
   integer                       , intent(in)    :: kx
   real   , dimension(imx,jmx,kx), intent(in)    :: t
   real   , dimension(imx,jmx,kx), intent(in)    :: pp
   real   , dimension(imx,jmx)   , intent(in)    :: ter
   real   , dimension(imx,jmx)   , intent(in)    :: sfp
   real   , dimension(imx,jmx)   , intent(inout) :: ts
   real   , dimension(imx,jmx)   , intent(inout) :: slp
   !----- Local variables. ----------------------------------------------------------------!
   integer                                       :: i
   integer                                       :: j
   integer                                       :: k
   integer                                       :: kupto
   integer                                       :: klo
   integer                                       :: khi
   logical                                       :: l1
   logical                                       :: l2
   logical                                       :: l3
   real   , dimension(imx,jmx)                   :: ps
   real   , dimension(imx,jmx)                   :: pl
   real   , dimension(imx,jmx)                   :: t0
   real   , dimension(imx,jmx)                   :: xklev
   real                                          :: xk
   real                                          :: xkhold
   real                                          :: plo
   real                                          :: phi
   real                                          :: tlo
   real                                          :: thi
   real                                          :: tl
   real                                          :: tbar
   real                                          :: t0hold
   real                                          :: hl
   !---------------------------------------------------------------------------------------!

   !------ Compute pressure at pconst mb above surface (pl). ------------------------------!
   kupto=kx/2


   mainloop: do
      do j=1, jmx
         do i=1,imx
            pl(i,j)=sfp(i,j)-pconst
            xklev(i,j)=0.
         end do
      end do

      !----- Find 2 levels on sigma surfaces surrounding pl at each i,j. ------------------!
      jloop: do j=1,jmx
         iloop: do i=1,imx
            kloop: do k=kx-1,kupto,-1
               xk     = real(k)
               xkhold = xklev(i,j)
               xklev(i,j) = merge(xk,xkhold                                                &
                                 ,((pp(i,j,k)   <  pl(i,j)) .and.                          &
                                   (pp(i,j,k+1) >= pl(i,j))       ))
            end do kloop

            if (xklev(i,j) < 1.) then
               write (unit=*,fmt='(a,1x,es12.5,1x,a)')                                     &
                  ' Error finding pressure level ',pconst,' mb above the surface!'
               write (unit=*,fmt='(a,1x,i5,a)') ' Last k level =',kupto,'...'

               if (kupto /= 1) then
                 write (unit=*,fmt='(a)') ' Trying again with kupto=1...'
                 kupto = 1
                 cycle mainloop
               else
                  write(unit=*,fmt='(a,1x,i5,1x)')     ' - I    =',i
                  write(unit=*,fmt='(a,1x,i5,1x)')     ' - J    =',i
                  write(unit=*,fmt='(a,1x,es12.5,1x)') ' - PL   =',pl(i,j)
                  write(unit=*,fmt='(a,1x,es12.5,1x)') ' - PSFC =',sfp(i,j)
                  stop
               end if
            end if
         end do iloop
      end do jloop
      !---- The default is to leave the loop... -------------------------------------------!
      exit mainloop
   end do mainloop

   !---------------------------------------------------------------------------------------!
   !      Get temperature at pl (tl), extrapolate t at surface (ts) and T at sea level     !
   ! (t0) with 6.5 k/km lapse rate.                                                        !
   !---------------------------------------------------------------------------------------!
   jloop2: do j=1,jmx
      iloop2: do i=1,imx
         klo     = nint(xklev(i,j))+1
         khi     = nint(xklev(i,j))
         plo     = pp(i,j,klo)
         phi     = pp(i,j,khi)
         tlo     = t(i,j,klo)
         thi     = t(i,j,khi)
         tl      = thi-(thi-tlo)*alog(pl(i,j)/phi)/alog(plo/phi)
         ts(i,j) = tl*(sfp(i,j)/pl(i,j))**xterm
         tbar    = (ts(i,j)+tl)*0.5
         hl      = ter(i,j)-rdry/grav*alog(pl(i,j)/sfp(i,j))*tbar
         t0(i,j) = tl+gamma*hl
      end do iloop2
   end do jloop2
   !---------------------------------------------------------------------------------------!



   !----- Correct sea level temperature if too hot. ---------------------------------------!
   jloop3: do j=1,jmx
      iloop3: do i=1,imx
         l1      = t0(i,j) <  tcrit
         l2      = ts(i,j) >= tcrit
         l3      = .not. l1
         t0hold  = t0(i,j)

         t0(i,j) = merge(t0hold,merge(tcrit,tcrit-0.005*(ts(i,j)-tcrit)**2,l2.and.l3)      &
                        ,l1.and.l2)
      end do iloop3
   end do jloop3
   !---------------------------------------------------------------------------------------!



   !----- Compute sea level pressure. -----------------------------------------------------!
   jloop4: do j=1,jmx
      iloop4: do i=1,imx
         slp(i,j)=sfp(i,j)*exp(2.*grav*ter(i,j)/(rdry*(ts(i,j)+t0(i,j))))
      end do iloop4
   end do jloop4
   !---------------------------------------------------------------------------------------!

   return
end subroutine seaprs_0
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the relative soil moisture.                                 !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_slmstf(nx,ny,nz,np,inwoutf,soil_text)
   use soil_coms, only : soil ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: np
   real, dimension(nx,ny,nz,np), intent(inout) :: inwoutf
   real, dimension(nx,ny,nz,np), intent(in)    :: soil_text
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   integer                                     :: p
   integer                                     :: nsoil
   real                                        :: soil_water
   real                                        :: soil_rmois
   !---------------------------------------------------------------------------------------!
   
   do p=2,np
      do z=1,nz
         do y=1,ny
            do x=1,nx
               !----- Copy point to temporary variable. -----------------------------------!
               soil_water = inwoutf(x,y,z,p)

               !----- Find the soil class. ------------------------------------------------!
               nsoil = nint(soil_text(x,y,z,p))

               select case (nsoil)
               case (0) ! Water
                  soil_rmois = 1.
               case default ! Soil
                  soil_rmois = (soil_water         - soil(nsoil)%soilcp)                   &
                             / (soil(nsoil)%slmsts - soil(nsoil)%soilcp)
               end select

               !----- Copy the result back to the array. ----------------------------------!
               inwoutf(x,y,z,p) = soil_rmois
            end do
         end do
      end do
   end do

   return
end subroutine RAMS_comp_slmstf
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the soil matric potential.                                  !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_smpot(nx,ny,nz,np,inwoutf,soil_text)
   use soil_coms , only : soil  ! ! intent(in)
   use rconstants, only : wdnsi & ! intent(in)
                        , grav  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: np
   real, dimension(nx,ny,nz,np), intent(inout) :: inwoutf
   real, dimension(nx,ny,nz,np), intent(in)    :: soil_text
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   integer                                     :: p
   integer                                     :: nsoil
   real                                        :: soil_water
   real                                        :: smpot
   !---------------------------------------------------------------------------------------!
   
   do p=2,np
      do z=1,nz
         do y=1,ny
            do x=1,nx
               !----- Copy point to temporary variable. -----------------------------------!
               soil_water = inwoutf(x,y,z,p)

               !----- Find the soil class. ------------------------------------------------!
               nsoil = nint(soil_text(x,y,z,p))

               select case (nsoil)
               case (0) ! Water
                  smpot = 0.
               case default ! Soil
                  smpot = grav * wdnsi * soil(nsoil)%slpots                                &
                          * (soil(nsoil)%slmsts / soil_water) ** soil(nsoil)%slbs
               end select

               !----- Copy the result back to the array. ----------------------------------!
               inwoutf(x,y,z,p) = smpot
            end do
         end do
      end do
   end do

   return
end subroutine RAMS_comp_smpot
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the temperature and liquid fraction given the internal      !
! energy and water content.                                                                !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_uextcm2tl(nx,ny,nz,np,inqoutt,inwoutl,soil_text)
   use rconstants, only : wdns ! ! intent(in)
   use soil_coms , only : soil ! ! intent(in)
   use therm_lib , only : uextcm2tl ! ! subroutine
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: np
   real, dimension(nx,ny,nz,np), intent(inout) :: inqoutt
   real, dimension(nx,ny,nz,np), intent(inout) :: inwoutl
   real, dimension(nx,ny,nz,np), intent(in)    :: soil_text
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   integer                                     :: p
   integer                                     :: nsoil
   real                                        :: energy
   real                                        :: water
   real                                        :: dryhcap
   real                                        :: temperature
   real                                        :: fracliq
   !---------------------------------------------------------------------------------------!

   !----- Loop through all grid points, skipping the water patch. -------------------------!
   do p=2,np
      do z=1,nz
         do y=1,ny
            do x=1,nx
               !----- Save variables into temporary places. -------------------------------!
               energy  = inqoutt(x,y,z,p)
               water   = inwoutl(x,y,z,p) * wdns
               nsoil   = nint(soil_text(x,y,z,p))
               dryhcap = soil(nsoil)%slcpd
               !----- Compute temperature and liquid water fraction. ----------------------!
               call uextcm2tl(energy,water,dryhcap,temperature,fracliq)
               !----- Save in the variables that will be returned. ------------------------!
               inqoutt(x,y,z,p) = temperature
               inwoutl(x,y,z,p) = fracliq
            end do
         end do
      end do
   end do

   return
end subroutine RAMS_comp_uextcm2tl
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_copysst(nx,ny,nz,inqoutt)
   use therm_lib , only : uint2tl  ! ! subroutine
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                    , intent(in)    :: nx
   integer                    , intent(in)    :: ny
   integer                    , intent(in)    :: nz
   real, dimension(nx,ny,nz)  , intent(inout) :: inqoutt
   !----- Local variables. ----------------------------------------------------------------!
   integer                                    :: x
   integer                                    :: y
   integer                                    :: z
   real                                       :: energy
   real, dimension(nx,ny)                     :: temperature
   real                                       :: fracliq
   !---------------------------------------------------------------------------------------!

   !----- Copy energy and compute the temperature, saving it into a 2-D array. ------------!
   do y=1,ny
      do x=1,nx
         energy = inqoutt(x,y,nz)
         call uint2tl(energy,temperature(x,y),fracliq)
      end do
   end do

   !----- Copy temperature to the output array. -------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
            inqoutt(x,y,z) = temperature(x,y)
         end do
      end do
   end do

   return
end subroutine RAMS_comp_copysst
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_tempC(nx,ny,nz,np,temp)
   use rconstants, only : t00 ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: np
   real, dimension(nx,ny,nz,np), intent(inout) :: temp
   !----- Local variables. ----------------------------------------------------------------!
   integer                                    :: x
   integer                                    :: y
   integer                                    :: z
   integer                                    :: p
   !---------------------------------------------------------------------------------------!

   do p=1,np
      do z=1,nz
         do y=1,ny
            do x=1,nx
               temp(x,y,z,p) = temp(x,y,z,p) - t00
            end do
         end do
      end do
   end do
   return
end subroutine RAMS_comp_tempC
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_pvap(nx,ny,nz,pres,inroute)
   use rconstants, only : ep ! ! intent(in) 
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(in)    :: pres
   real   , dimension(nx,ny,nz), intent(inout) :: inroute
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   real                                        :: rvap
   real                                        :: pvap
   !---------------------------------------------------------------------------------------!

   do x=1,nx
      do y=1,ny
         do z=1,nz
            rvap           = inroute(x,y,z)
            pvap           = pres(x,y,z) * rvap / (ep + rvap)
            inroute(x,y,z) = pvap
         end do
      end do
   end do

   return
end subroutine RAMS_comp_pvap
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_spvol(nx,ny,nz,intouta,pvap,pres)
   use rconstants, only : ep   & ! intent(in) 
                        , rdry ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(inout) :: intouta
   real   , dimension(nx,ny,nz), intent(in)    :: pvap
   real   , dimension(nx,ny,nz), intent(in)    :: pres
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   real                                        :: temp
   real                                        :: alpha
   !---------------------------------------------------------------------------------------!

   do z=1,nz
      do y=1,ny
         do x=1,nx
            temp           = intouta(x,y,z)
            alpha          = rdry * temp /(pres(x,y,z) - (1.-ep) * pvap(x,y,z))
            intouta(x,y,z) = alpha
         end do
      end do
   end do

   return
end subroutine RAMS_comp_spvol
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_theta2temp(nx,ny,nz,inthoutt,press)
   use rconstants, only : p00i & ! intent(in)
                        , rocp ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(inout) :: inthoutt
   real   , dimension(nx,ny,nz), intent(in)    :: press
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   real                                        :: theta
   real                                        :: temp
   !---------------------------------------------------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
            theta           = inthoutt(x,y,z)
            temp            = theta * (p00i * press(x,y,z)) ** rocp
            inthoutt(x,y,z) = temp
         end do
      end do
   end do
   return
end subroutine RAMS_comp_theta2temp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_zenith(nx,ny,cosz,zenith)
   use rconstants, only : onerad ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in)    :: nx
   integer                  , intent(in)    :: ny
   real   , dimension(nx,ny), intent(in)    :: cosz
   real   , dimension(nx,ny), intent(inout) :: zenith
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   real                                        :: theta
   real                                        :: temp
   !---------------------------------------------------------------------------------------!

   do y=1,ny
      do x=1,nx
         zenith(x,y) = acos(cosz(x,y)) * onerad
         if (zenith(x,y) > 90. ) zenith(x,y) = 90.
      end do
   end do

   return
end subroutine RAMS_comp_zenith
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_wflx2latent(nx,ny,np,wflx,temp)
   use rconstants, only : t3ple & ! intent(in)
                        , alvl3 & ! intent(in)
                        , alvi3 ! ! intent(in)
   use therm_lib , only : alvl  & ! function
                        , alvi  ! ! function
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: np
   real   , dimension(nx,ny,np), intent(inout) :: wflx
   real   , dimension(nx,ny,np), intent(in)    :: temp
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: p
   integer                                     :: x
   integer                                     :: y
   real                                        :: latent
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   do p=1,np
      do y=1,ny
         do x=1,nx
            if (temp(x,y,p) == t3ple) then
               latent = 0.5 * (alvl3 + alvi3)
            elseif (temp(x,y,p) > t3ple) then
               latent = alvl(temp(x,y,p))
            else
               latent = alvi(temp(x,y,p))
            end if
            wflx(x,y,p) = wflx(x,y,p) * latent
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine RAMS_comp_wflx2latent
!==========================================================================================!
!==========================================================================================!







!==========================================================================================!
!==========================================================================================!
!     This subroutine clones an array.                                                     !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_clone(nx,ny,nz,dest,src)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   real   , dimension(nx,ny,nz), intent(inout) :: dest
   real   , dimension(nx,ny,nz), intent(in)    :: src
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
            dest(x,y,z) = src(x,y,z)
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine RAMS_comp_clone
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine computes the rotated versions of zonal and meridional components.   !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_rotate(nx,ny,nz,uwnd,vwnd,ngrd)
   use somevars, only : myplatn & ! intent(in)
                      , myplonn & ! intent(in)
                      , myxtn   & ! intent(in)
                      , myytn   ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: nx
   integer                          , intent(in)    :: ny
   integer                          , intent(in)    :: nz
   real(kind=4), dimension(nx,ny,nz), intent(inout) :: uwnd
   real(kind=4), dimension(nx,ny,nz), intent(inout) :: vwnd
   integer                          , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: x
   integer                                          :: y
   integer                                          :: z
   real                                             :: qlat
   real                                             :: qlon
   real                                             :: utmp
   real                                             :: vtmp
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Main loop.                                                                       !
   !---------------------------------------------------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
            call xy_ll(qlat,qlon,myplatn(ngrd),myplonn(ngrd),myxtn(x,ngrd),myytn(y,ngrd))
            utmp = uwnd(x,y,z)
            vtmp = vwnd(x,y,z)
            call uvtoueve(utmp,vtmp,uwnd(x,y,z),vwnd(x,y,z),qlat,qlon,myplatn(ngrd)        &
                         ,myplonn(ngrd))
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!
   return
end subroutine RAMS_comp_rotate
!==========================================================================================!
!==========================================================================================!








!==========================================================================================!
!==========================================================================================!
!     This subroutine clones an array.                                                     !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_uvsfc(nx,ny,uorv,uvsfc,u,v)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                           , intent(in)    :: nx
   integer                           , intent(in)    :: ny
   character(len=1)                  , intent(in)    :: uorv
   real            , dimension(nx,ny), intent(inout) :: uvsfc
   real            , dimension(nx,ny), intent(in)    :: u
   real            , dimension(nx,ny), intent(in)    :: v
   !----- Local variables. ----------------------------------------------------------------!
   integer                                  :: x
   integer                                  :: y
   real                                     :: angle
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Select which component is sought.                                                  !
   !---------------------------------------------------------------------------------------!
   select case(uorv)
   case ('U','u','X','x')
      do y=1,ny
         do x=1,nx
            uvsfc(x,y) = uvsfc(x,y) * cos(atan2(v(x,y),u(x,y)))
         end do
      end do

   case('V','v','Y','y')
      do y=1,ny
         do x=1,nx
            uvsfc(x,y) = uvsfc(x,y) * sin(atan2(v(x,y),u(x,y)))
         end do
      end do

   case default
      write(unit=*,fmt="(a)")      " Invalid UORV entry!"
      write(unit=*,fmt="(a,1x,a)") "UORV =",uorv
      stop 

   end select
   !---------------------------------------------------------------------------------------!

   return
end subroutine RAMS_comp_uvsfc
!==========================================================================================!
!==========================================================================================!








!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the precipitable water.                                     !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_prwtr(nx,ny,nz,a,rtp,dn0,topt,ngrd)
   use somevars, only : myzmn  & ! intent(in)
                      , mynnzp ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                              , intent(in)    :: nx
   integer                              , intent(in)    :: ny
   integer                              , intent(in)    :: nz
   real            , dimension(nx,ny)   , intent(inout) :: a
   real            , dimension(nx,ny,nz), intent(in)    :: rtp
   real            , dimension(nx,ny,nz), intent(in)    :: dn0
   real            , dimension(nx,ny)   , intent(in)    :: topt
   integer                              , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   integer                                              :: x
   integer                                              :: y
   integer                                              :: z
   real                                                 :: ztop
   real                                                 :: rtgt
   real                                                 :: qtp
   !---------------------------------------------------------------------------------------!



   !----- Define grid bottom and top. -----------------------------------------------------!
   ztop = myzmn(mynnzp(1)-1,1)

   yloop: do y = 1,ny
      xloop: do x = 1,nx
         rtgt = 1. - topt(x,y) / ztop
         a(x,y) = 0.
         zloop: do z = 2,nz
             qtp    = rtp(x,y,z) / ( rtp(x,y,z) + 1.0 )
             a(x,y) = a(x,y) + qtp * dn0(x,y,z) * (myzmn(z,ngrd)-myzmn(z-1,ngrd)) * rtgt
         end do zloop
      end do xloop
   end do yloop

   return
end subroutine RAMS_comp_prwtr
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_cine_cape(nx,ny,nz,zref,which,energy,press,exner,theiv,rtp,rv,tempk   &
                              ,topt,ngrd)
   use therm_lib , only : thetaeiv2thil & ! function
                        , thil2tqall    & ! subroutine
                        , idealdens     ! ! function
   use somevars  , only : myzmn         & ! intent(in)
                        , mynnzp        ! ! intent(in)
   use rconstants, only : grav          & ! intent(in)
                        , cpdryi        ! ! intent(in)
   use rout_coms , only : undefflg      ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                              , intent(in)    :: nx
   integer                              , intent(in)    :: ny
   integer                              , intent(in)    :: nz
   integer                              , intent(in)    :: zref
   character(len=4)                     , intent(in)    :: which
   real            , dimension(nx,ny)   , intent(inout) :: energy
   real            , dimension(nx,ny,nz), intent(in)    :: press
   real            , dimension(nx,ny,nz), intent(in)    :: exner
   real            , dimension(nx,ny,nz), intent(in)    :: theiv
   real            , dimension(nx,ny,nz), intent(in)    :: rtp
   real            , dimension(nx,ny,nz), intent(in)    :: rv
   real            , dimension(nx,ny,nz), intent(in)    :: tempk
   real            , dimension(nx,ny)   , intent(in)    :: topt
   integer                              , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   integer                                              :: x
   integer                                              :: y
   integer                                              :: z
   real                                                 :: dz
   real                                                 :: ztop
   real                                                 :: rtgt
   real                                                 :: cine
   real                                                 :: cape
   real                                                 :: theiv_p
   real                                                 :: rtp_p
   real                                                 :: thil_p
   real                                                 :: tempk_p
   real                                                 :: rvap_p
   real                                                 :: rliq_p
   real                                                 :: rice_p
   real                                                 :: rsat_p
   real                                                 :: rho_p
   real                                                 :: rho_e
   logical                                              :: lfc
   !---------------------------------------------------------------------------------------!




   !----- Define top of the grid. ---------------------------------------------------------!
   ztop = myzmn(mynnzp(1)-1,1)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Loop over domain.                                                                 !
   !---------------------------------------------------------------------------------------!
   yloop: do y=1,ny
      xloop: do x=1,nx
         !------ Initiliase the properties for this layer. --------------------------------!
         rtgt = 1. - topt(x,y) / ztop
         lfc  = .false.
         cine = 0.
         cape = 0.
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     A parcel retains its theta_Eiv and rtp.                                     !
         !---------------------------------------------------------------------------------!
         theiv_p = theiv(x,y,zref)
         rtp_p   = rtp  (x,y,zref)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Loop over this layer.                                                       !
         !---------------------------------------------------------------------------------!
         zloop: do z = zref+1,nz
            dz = ( myzmn(z,ngrd) - myzmn(z-1,ngrd) ) * rtgt

            !----- Find the parcel temperature and vapour mixing ratio. -------------------!
            thil_p  = thetaeiv2thil(theiv_p,press(x,y,z),rtp_p)
            tempk_p = exner(x,y,z) * thil_p * cpdryi
            call thil2tqall(thil_p,exner(x,y,z),press(x,y,z),rtp_p,rliq_p,rice_p,tempk_p   &
                           ,rvap_p,rsat_p)
            rho_p = idealdens(press(x,y,z),tempk_p,rvap_p,rtp_p)
            rho_e = idealdens(press(x,y,z),tempk(x,y,z),rv(x,y,z),rtp(x,y,z))
            !------------------------------------------------------------------------------!


            !----- Find out whether to add to CINE or CAPE. -------------------------------!
            if (rho_p < rho_e) then
               !----- Buoyant.  Add energy to CAPE. ---------------------------------------!
               lfc  = .true.
               cape = cape + grav * (1. - rho_p / rho_e) * dz
               !---------------------------------------------------------------------------!
            else if (lfc) then
               !----- No longer buoyant, reached LNB, move on. ----------------------------!
               exit zloop
               !---------------------------------------------------------------------------!
            else
               !----- Beneath LFC, add to inhibition. -------------------------------------!
               cine = cine + grav * (rho_p / rho_e - 1.) * dz
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do zloop
         !---------------------------------------------------------------------------------!



         !----- Make CAPE and CINE undefined in case LFC was never found. -----------------!
         if (.not. lfc) then
            cine = undefflg
            cape = undefflg
         end if
         !---------------------------------------------------------------------------------!


         !----- Find out which energy goes to output. -------------------------------------!
         select case(which)
         case ('cape')
            energy(x,y) = cape
         case ('cine')
            energy(x,y) = cine
         end select
         !---------------------------------------------------------------------------------!
      end do xloop
   end do yloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine  RAMS_comp_cine_cape
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine calculates the SHOWALTER index.                                     !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_showalter(nx,ny,nz,showalter,press,exner,theiv,rtp,rv,tempk,topt,ngrd)
   use therm_lib , only : thetaeiv2thil & ! function
                        , thil2tqall    & ! subroutine
                        , idealdens     & ! function
                        , press2exner   ! ! function
   use rconstants, only : grav          & ! intent(in)
                        , cpdryi        ! ! intent(in)
   use rout_coms , only : undefflg      ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                              , intent(in)    :: nx
   integer                              , intent(in)    :: ny
   integer                              , intent(in)    :: nz
   real            , dimension(nx,ny)   , intent(inout) :: showalter
   real            , dimension(nx,ny,nz), intent(in)    :: press
   real            , dimension(nx,ny,nz), intent(in)    :: exner
   real            , dimension(nx,ny,nz), intent(in)    :: theiv
   real            , dimension(nx,ny,nz), intent(in)    :: rtp
   real            , dimension(nx,ny,nz), intent(in)    :: rv
   real            , dimension(nx,ny,nz), intent(in)    :: tempk
   real            , dimension(nx,ny)   , intent(in)    :: topt
   integer                              , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   integer                                              :: x
   integer                                              :: y
   integer                                              :: z
   real                                                 :: rtgt
   real                                                 :: theiv_p
   real                                                 :: rtp_p
   real                                                 :: thil_p
   real                                                 :: tempk_p
   real                                                 :: rvap_p
   real                                                 :: rliq_p
   real                                                 :: rice_p
   real                                                 :: rsat_p
   real                                                 :: tamb_500
   real                                                 :: epress
   real                                                 :: exn500
   !----- Constants. ----------------------------------------------------------------------!
   real                                 , parameter     :: p850 = 85000.
   real                                 , parameter     :: p500 = 50000.
   !---------------------------------------------------------------------------------------!



   !----- Exner function at 500hPa. -------------------------------------------------------!
   exn500 = press2exner(p500)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Loop over domain.                                                                 !
   !---------------------------------------------------------------------------------------!
   yloop: do y=1,ny
      xloop: do x=1,nx

         showalter(x,y) = undefflg

         !----- The Showalter index can't be determined at places above 850hPa. -----------!
         if (press(x,y,2) >= p850) then

            !------------------------------------------------------------------------------!
            !     Loop over height coordinates.                                            !
            !------------------------------------------------------------------------------!
            zloop: do z=3,nz
               if (press(x,y,z-1) >= p850 .and. press(x,y,z) < p850) then
                  if (press(x,y,z-1) == p850) then
                     !------ Grab data from 850hPa. ---------------------------------------!
                     theiv_p = theiv(x,y,z-1)
                     rtp_p   = rtp  (x,y,z-1)
                     !---------------------------------------------------------------------!
                  else
                     !------ Interpolate thetae and rtp to 850hPa. ------------------------!
                     epress  = log(p850/press(x,y,z-1))/log(press(x,y,z)/press(x,y,z-1))
                     theiv_p = theiv(x,y,z-1) * (theiv(x,y,z)/theiv(x,y,z-1))**epress
                     rtp_p   = rtp  (x,y,z-1) * (rtp  (x,y,z)/rtp  (x,y,z-1))**epress
                     !---------------------------------------------------------------------!
                  end if



                  !------ Find the values at 500 hPa. -------------------------------------!
                  thil_p  = thetaeiv2thil(theiv_p,p500,rtp_p)
                  tempk_p = exn500 * thil_p * cpdryi
                  call thil2tqall(thil_p,exn500,p500,rtp_p,rliq_p,rice_p,tempk_p,rvap_p    &
                                 ,rsat_p)
                  !------------------------------------------------------------------------!
               elseif (press(x,y,z-1) >= p500 .and. press(x,y,z) < p500) then
                  !----- Interpolate temperature to 500 hPa. ------------------------------!
                  if (press(x,y,z-1) == p500) then
                     tamb_500 = tempk(x,y,z-1)
                  else
                     !------ Interpolate temperature to 850hPa. ---------------------------!
                     epress   = log(p500/press(x,y,z-1))/log(press(x,y,z)/press(x,y,z-1))
                     tamb_500 = tempk(x,y,z-1) * (tempk(x,y,z)/tempk(x,y,z-1))**epress
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!

                  !----- Find the Showalter index. ----------------------------------------!
                  showalter(x,y) = tamb_500 - tempk_p
                  exit zloop
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end do zloop
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do xloop
   end do yloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine  RAMS_comp_showalter
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine computes the SWEAT index.                                           !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_sweat(nx,ny,nz,sweat,press,tempk,tdewk,uwnd,vwnd,uvel)
   use rout_coms , only : undefflg      ! ! intent(in)
   use rconstants, only : t00           ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: nx
   integer                          , intent(in)    :: ny
   integer                          , intent(in)    :: nz
   real(kind=4), dimension(nx,ny)   , intent(inout) :: sweat
   real(kind=4), dimension(nx,ny,nz), intent(inout) :: press
   real(kind=4), dimension(nx,ny,nz), intent(inout) :: tempk
   real(kind=4), dimension(nx,ny,nz), intent(inout) :: tdewk
   real(kind=4), dimension(nx,ny,nz), intent(inout) :: uwnd
   real(kind=4), dimension(nx,ny,nz), intent(inout) :: vwnd
   real(kind=4), dimension(nx,ny,nz), intent(inout) :: uvel
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: x
   integer                                          :: y
   integer                                          :: z
   real(kind=4)                                     :: temp_850
   real(kind=4)                                     :: temp_500
   real(kind=4)                                     :: tdew_850
   real(kind=4)                                     :: tdew_500
   real(kind=4)                                     :: uvel_850
   real(kind=4)                                     :: uvel_500
   real(kind=4)                                     :: uwnd_850
   real(kind=4)                                     :: uwnd_500
   real(kind=4)                                     :: vwnd_850
   real(kind=4)                                     :: vwnd_500
   real(kind=4)                                     :: udir_850
   real(kind=4)                                     :: udir_500
   real(kind=4)                                     :: ttotals
   real(kind=4)                                     :: wshear
   real(kind=4)                                     :: epress
   !----- Constants. ----------------------------------------------------------------------!
   real                             , parameter     :: p850 = 85000.
   real                             , parameter     :: p500 = 50000.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     First loop, wind magnitude.                                                       !
   !---------------------------------------------------------------------------------------!
   do z=1,nz
      do y=1,ny
         do x=1,nx
            uvel(x,y,z) = sqrt(uwnd(x,y,z)*uwnd(x,y,z)+vwnd(x,y,z)*vwnd(x,y,z))
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Second loop, find interpolated data.                                              !
   !---------------------------------------------------------------------------------------!
   do y=1,ny
      do x=1,nx
         !----- Initialise SWEAT with missing value at higher locations. ------------------!
         sweat(x,y) = undefflg
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Skip points when the surface is above 850hPa.                               !
         !---------------------------------------------------------------------------------!
         if (press(x,y,2) >= p850) then
            zloop: do z=1,nz
               if (press(x,y,z-1) >= p850 .and. press(x,y,z) < p850) then
                  if (press(x,y,z-1) == p850) then
                     !------ Grab data from 850hPa. ---------------------------------------!
                     temp_850 = tempk(x,y,z-1)
                     tdew_850 = tdewk(x,y,z-1)
                     uwnd_850 = uwnd (x,y,z-1)
                     vwnd_850 = vwnd (x,y,z-1)
                     uvel_850 = uvel (x,y,z-1)
                     !---------------------------------------------------------------------!
                  else
                     !------ Interpolate thetae and rtp to 850hPa. ------------------------!
                     epress  = log(p850/press(x,y,z-1))/log(press(x,y,z)/press(x,y,z-1))
                     temp_850 = tempk(x,y,z-1) * ( tempk(x,y,z) / tempk(x,y,z-1) )**epress
                     tdew_850 = tdewk(x,y,z-1) * ( tdewk(x,y,z) / tdewk(x,y,z-1) )**epress
                     uvel_850 = uvel (x,y,z-1) * ( uvel (x,y,z) / uvel (x,y,z-1) )**epress
                     uwnd_850 = uwnd (x,y,z-1) * ( uwnd (x,y,z) - uwnd (x,y,z-1) ) *epress
                     vwnd_850 = vwnd (x,y,z-1) * ( vwnd (x,y,z) - vwnd (x,y,z-1) ) *epress
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!

                  !----- Wind direction at 850 hPa. ---------------------------------------!
                  udir_850 = atan2(vwnd_850,uwnd_850)
                  !------------------------------------------------------------------------!
               elseif (press(x,y,z-1) >= p500 .and. press(x,y,z) < p500) then
                  if (press(x,y,z-1) == p500) then
                     !------ Grab data from 850hPa. ---------------------------------------!
                     temp_500 = tempk(x,y,z-1)
                     tdew_500 = tdewk(x,y,z-1)
                     uwnd_500 = uwnd (x,y,z-1)
                     vwnd_500 = vwnd (x,y,z-1)
                     uvel_500 = uvel (x,y,z-1)
                     !---------------------------------------------------------------------!
                  else
                     !------ Interpolate thetae and rtp to 850hPa. ------------------------!
                     epress   = log(p500/press(x,y,z-1))/log(press(x,y,z)/press(x,y,z-1))
                     temp_500 = tempk(x,y,z-1) * ( tempk(x,y,z) / tempk(x,y,z-1) )**epress
                     tdew_500 = tdewk(x,y,z-1) * ( tdewk(x,y,z) / tdewk(x,y,z-1) )**epress
                     uvel_500 = uvel (x,y,z-1) * ( uvel (x,y,z) / uvel (x,y,z-1) )**epress
                     uwnd_500 = uwnd (x,y,z-1) * ( uwnd (x,y,z) - uwnd (x,y,z-1) ) *epress
                     vwnd_500 = vwnd (x,y,z-1) * ( vwnd (x,y,z) - vwnd (x,y,z-1) ) *epress
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!

                  !----- Wind direction at 850 hPa. ---------------------------------------!
                  udir_500 = atan2(vwnd_500,uwnd_500)
                  !------------------------------------------------------------------------!


                  !----- Everything has been found, leave the zloop. ----------------------!
                  exit zloop
                  !------------------------------------------------------------------------!

               end if
               !---------------------------------------------------------------------------!
            end do zloop
            !------------------------------------------------------------------------------!


            !----- Find the Total totals. -------------------------------------------------!
            ttotals = max(49., temp_850 + tdew_850 - 2. * temp_500)
            !------------------------------------------------------------------------------!


            !----- Find the shear term. ---------------------------------------------------!
            wshear  = max(0., sin((udir_500-udir_850) + 0.2))
            !------------------------------------------------------------------------------!


            !----- Find the SWEAT index. --------------------------------------------------!
            sweat(x,y) =  12. * max(0., (tdew_850 - t00)) +  20. * (ttotals - 49.)         &
                       +   2. * uvel_850 * 1.94384 + uvel_500 * 1.94384 + 125. * wshear
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!


   return
end subroutine RAMS_comp_sweat
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine computes the SWEAT index.                                           !
!------------------------------------------------------------------------------------------!
subroutine RAMS_comp_kindex(nx,ny,nz,kindex,press,tempk,tdewk)
   use rout_coms , only : undefflg      ! ! intent(in)
   use rconstants, only : t00           ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                          , intent(in)    :: nx
   integer                          , intent(in)    :: ny
   integer                          , intent(in)    :: nz
   real(kind=4), dimension(nx,ny)   , intent(inout) :: kindex
   real(kind=4), dimension(nx,ny,nz), intent(inout) :: press
   real(kind=4), dimension(nx,ny,nz), intent(inout) :: tempk
   real(kind=4), dimension(nx,ny,nz), intent(inout) :: tdewk
   !----- Local variables. ----------------------------------------------------------------!
   integer                                          :: x
   integer                                          :: y
   integer                                          :: z
   real(kind=4)                                     :: temp_850
   real(kind=4)                                     :: temp_700
   real(kind=4)                                     :: temp_500
   real(kind=4)                                     :: tdew_850
   real(kind=4)                                     :: tdew_700
   real(kind=4)                                     :: tdew_500
   real(kind=4)                                     :: epress
   !----- Constants. ----------------------------------------------------------------------!
   real                             , parameter     :: p850 = 85000.
   real                             , parameter     :: p700 = 70000.
   real                             , parameter     :: p500 = 50000.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Second loop, find interpolated data.                                              !
   !---------------------------------------------------------------------------------------!
   do y=1,ny
      do x=1,nx
         !----- Initialise K-index with missing value at higher locations. ----------------!
         kindex(x,y) = undefflg
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Skip points when the surface is above 850hPa.                               !
         !---------------------------------------------------------------------------------!
         if (press(x,y,2) >= p850) then
            zloop: do z=1,nz
               if (press(x,y,z-1) >= p850 .and. press(x,y,z) < p850) then
                  if (press(x,y,z-1) == p850) then
                     !------ Grab data from 850hPa. ---------------------------------------!
                     temp_850 = tempk(x,y,z-1)
                     tdew_850 = tdewk(x,y,z-1)
                     !---------------------------------------------------------------------!
                  else
                     !------ Interpolate thetae and rtp to 850hPa. ------------------------!
                     epress  = log(p850/press(x,y,z-1))/log(press(x,y,z)/press(x,y,z-1))
                     temp_850 = tempk(x,y,z-1) * ( tempk(x,y,z) / tempk(x,y,z-1) )**epress
                     tdew_850 = tdewk(x,y,z-1) * ( tdewk(x,y,z) / tdewk(x,y,z-1) )**epress
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!

               elseif (press(x,y,z-1) >= p700 .and. press(x,y,z) < p700) then
                  if (press(x,y,z-1) == p700) then
                     !------ Grab data from 700hPa. ---------------------------------------!
                     temp_700 = tempk(x,y,z-1)
                     tdew_700 = tdewk(x,y,z-1)
                     !---------------------------------------------------------------------!
                  else
                     !------ Interpolate thetae and rtp to 850hPa. ------------------------!
                     epress   = log(p700/press(x,y,z-1))/log(press(x,y,z)/press(x,y,z-1))
                     temp_700 = tempk(x,y,z-1) * ( tempk(x,y,z) / tempk(x,y,z-1) )**epress
                     tdew_700 = tdewk(x,y,z-1) * ( tdewk(x,y,z) / tdewk(x,y,z-1) )**epress
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!

               elseif (press(x,y,z-1) >= p500 .and. press(x,y,z) < p500) then
                  if (press(x,y,z-1) == p500) then
                     !------ Grab data from 850hPa. ---------------------------------------!
                     temp_500 = tempk(x,y,z-1)
                     tdew_500 = tdewk(x,y,z-1)
                     !---------------------------------------------------------------------!
                  else
                     !------ Interpolate thetae and rtp to 850hPa. ------------------------!
                     epress   = log(p500/press(x,y,z-1))/log(press(x,y,z)/press(x,y,z-1))
                     temp_500 = tempk(x,y,z-1) * ( tempk(x,y,z) / tempk(x,y,z-1) )**epress
                     tdew_500 = tdewk(x,y,z-1) * ( tdewk(x,y,z) / tdewk(x,y,z-1) )**epress
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!

                  !----- Everything has been found, leave the zloop. ----------------------!
                  exit zloop
                  !------------------------------------------------------------------------!

               end if
               !---------------------------------------------------------------------------!
            end do zloop
            !------------------------------------------------------------------------------!


            !----- Find the Total totals. -------------------------------------------------!
            kindex(x,y) = temp_850 - temp_500 + tdew_850 - temp_700 + tdew_700 - t00
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!


   return
end subroutine RAMS_comp_kindex
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine RAMS_comp_massint(nx,ny,nz,nc,mass,rho,rmix,beneath,aloft,topt,ngrd)
   use somevars, only : myzmn   & ! intent(in)
                      , mynnzp  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                                 , intent(in)    :: nx
   integer                                 , intent(in)    :: ny
   integer                                 , intent(in)    :: nz
   integer                                 , intent(in)    :: nc
   real            , dimension(nx,ny,nc)   , intent(inout) :: mass
   real            , dimension(nx,ny,nz)   , intent(in)    :: rho
   real            , dimension(nx,ny,nz,nc), intent(in)    :: rmix
   real            , dimension(nx,ny)      , intent(in)    :: beneath
   real            , dimension(nx,ny)      , intent(in)    :: aloft
   real            , dimension(nx,ny)      , intent(in)    :: topt
   integer                                 , intent(in)    :: ngrd
   !----- Local variables. ----------------------------------------------------------------!
   integer                                                 :: x
   integer                                                 :: y
   integer                                                 :: z
   integer                                                 :: c
   real                                                    :: dz
   real                                                    :: ztop
   real                                                    :: rtgt
   real                                                    :: zlow
   real                                                    :: zhigh
   !---------------------------------------------------------------------------------------!




   !----- Define top of the grid. ---------------------------------------------------------!
   ztop = myzmn(mynnzp(1)-1,1)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Loop over x and y.                                                                !
   !---------------------------------------------------------------------------------------!
   yloop: do y=1,ny
      xloop: do x=1,nx
         !------ Initiliase the properties for this layer. --------------------------------!
         rtgt = 1. - topt(x,y) / ztop
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Cloud loop.                                                                 !
         !---------------------------------------------------------------------------------!
         cloop: do c=1,nc

            !------------------------------------------------------------------------------!
            !       Initialise integral, loop over height and add the layers of interest.  !
            !------------------------------------------------------------------------------!
            mass(x,y,c) = 0.
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Loop over this layer.                                                    !
            !------------------------------------------------------------------------------!
            zloop: do z = 2,nz
               zlow  = myzmn(z-1,ngrd) * rtgt
               zhigh = myzmn(  z,ngrd) * rtgt

               !---------------------------------------------------------------------------!
               !     Check whether to add this layer.                                      !
               !---------------------------------------------------------------------------!
               if (zlow >= beneath(x,y) .and. zhigh <= aloft(x,y)) then
                  dz          = min(zhigh,aloft(x,y)) - max(zlow,beneath(x,y))
                  mass(x,y,c) = mass(x,y,c) + rho(x,y,z) * rmix(x,y,z,c) * dz
               else if (zlow > aloft(x,y)) then
                  exit zloop
               end if
               !---------------------------------------------------------------------------!
            end do zloop
            mass(x,y,c) = max(0.,mass(x,y,c))
            !------------------------------------------------------------------------------!
         end do cloop
         !---------------------------------------------------------------------------------!
      end do xloop
      !------------------------------------------------------------------------------------!
   end do yloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine RAMS_comp_massint
!==========================================================================================!
!==========================================================================================!
