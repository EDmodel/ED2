!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!   08/31/08 - MLO - Switching micro by OLAM equivalent, which uses densities rather than  !
!                    mixing ratio.                                                         !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    Microphysics wrap subroutine, which will compute the hydrometeor properties and up-   !
! date the thermodynamic variables.                                                        !
!------------------------------------------------------------------------------------------!
subroutine micro_driver()

   use grid_dims, only:  &
           maxgrds       ! ! intent(in)

   use mem_basic, only : &
           basic_g       ! ! intent(inout)

   use mem_micro, only:  &
           micro_g       ! ! intent(inout)

   use mem_grid, only:   &
           ngrids        & ! intent(in)
          ,ngrid         & ! intent(in)
          ,zm            & ! intent(in)
          ,dzt           & ! intent(in)
          ,dtlt          & ! intent(in)
          ,jdim          & ! intent(in)
          ,maxnzp        & ! intent(in)
          ,time          & ! intent(in)
          ,dtlongn       & ! intent(in)
          ,zt            & ! intent(in)
          ,if_adap       & ! intent(in)
          ,grid_g        ! ! intent(in)

   use node_mod, only :  &
           mmzp          & ! intent(in)
          ,mzp           & ! intent(in)
          ,mxp           & ! intent(in)
          ,myp           & ! intent(in)
          ,ja            & ! intent(in)
          ,jz            & ! intent(in)
          ,ia            & ! intent(in)
          ,iz            & ! intent(in)
          ,mynum         ! ! intent(in)

   use micro_coms, only: &
           pcp_tab       ! ! intent(inout)

   use micphys, only:    &
           nhcat         & ! intent(in)
          ,ch3           & ! intent(out)
          ,pwvt          & ! intent(in)
          ,pwmasi        & ! intent(in)
          ,cdp1          & ! intent(out)
          ,pwvtmasi      & ! intent(out)
          ,icloud        & ! intent(in)
          ,irain         & ! intent(in)
          ,ipris         & ! intent(in)
          ,isnow         & ! intent(in)
          ,iaggr         & ! intent(in)
          ,igraup        & ! intent(in)
          ,ihail         ! ! intent(in)
        
   implicit none

   !----- Local Variables: ----------------------------------------------------------------!
   integer                            :: ngr,lhcat,i,j
   logical, dimension(maxgrds), save  :: ncall2g = .true.
   real                               :: dtlti 
   !---------------------------------------------------------------------------------------!


   !----- Sedimentation and homogeneous freezing table are build only once per grid -------!
   if (ncall2g(ngrid)) then

      call mksedim_tab(mzp,zm,dzt                       , pcp_tab(ngrid)%pcpfillc(1,1,1,1) &
                      ,pcp_tab(ngrid)%pcpfillr(1,1,1,1), pcp_tab(ngrid)%sfcpcp     (1,1,1) )

      call homfrzcl(dtlt,ngrid)

      ncall2g(ngrid) = .false.
   end if
   !---------------------------------------------------------------------------------------!

   call each_call(mzp,dtlt)


   dtlti = 1. / dtlt
   ngr   = ngrid

   do j = ja,jz
      do i = ia,iz

         call range_check(mzp,i,j                        , grid_g(ngr)%flpw       (i,j)    &
                         ,basic_g(ngr)%thp     (1,i,j)   , basic_g(ngr)%theta   (1,i,j)    &
                         ,basic_g(ngr)%pp      (1,i,j)   , basic_g(ngr)%rtp     (1,i,j)    &
                         ,basic_g(ngr)%rv      (1,i,j)   , basic_g(ngr)%wp      (1,i,j)    &
                         ,basic_g(ngr)%dn0     (1,i,j)   , basic_g(ngr)%pi0     (1,i,j)    &
                         ,micro_g(ngr))

         call mcphys_main(mzp,ngr,mynum,if_adap,dtlt,dtlti,time,zm,dzt,zt                  &
                         ,grid_g(ngr)%rtgt          (i,j), micro_g(ngr)%pcpg         (i,j) &
                         ,micro_g(ngr)%qpcpg        (i,j), micro_g(ngr)%dpcpg        (i,j) &
                         ,pcp_tab(ngr)%pcpfillc (1,1,1,1), pcp_tab(ngr)%pcpfillr (1,1,1,1) &
                         ,pcp_tab(ngr)%sfcpcp     (1,1,1))

         call copyback(mzp,i,j                        ,basic_g(ngr)%thp     (1,i,j)        &
                      ,basic_g(ngr)%theta   (1,i,j)   ,basic_g(ngr)%rtp     (1,i,j)        &
                      ,basic_g(ngr)%dn0     (1,i,j)   ,micro_g(ngr)                        )

      end do
   end do

   return
end subroutine micro_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This is the main mcphys subroutine, which will update the hydrometeor mixing ratios   !
! and, if requested, the concentrations. It will also compute the tendency on enthropy     !
! (ice-liquid potential temperature) due to precipitation.                                 !
!------------------------------------------------------------------------------------------!
subroutine mcphys_main(m1,ngr,mynum,if_adap,dtlt,dtlti,time,zm,dzt,zt,rtgt,pcpg,qpcpg      &
                      ,dpcpg,pcpfillc,pcpfillr,sfcpcp)

   use micphys    , only : &
           k1              & ! intent(in)
          ,k2              & ! intent(in)
          ,k3              & ! intent(in)
          ,availcat        & ! intent(in)
          ,progncat        & ! intent(in)
          ,nembfall        & ! intent(in)
          ,maxkfall        & ! intent(in)
          ,lpw             & ! intent(in)
          ,nhcat           & ! intent(in)
          ,jnmb            & ! intent(in)
          ,ncat            & ! intent(in)
          ,thil            & ! intent(out)
          ,ch1             & ! intent(out)
          ,cfvt            & ! intent(in)
          ,dsed_thil       & ! intent(out)
          ,accpx           & ! intent(out)
          ,pcprx           & ! intent(out)
          ,k1cnuc          & ! intent(inout)
          ,k2cnuc          & ! intent(inout)
          ,k1pnuc          & ! intent(inout)
          ,k2pnuc          ! ! intent(inout)
 
    use micro_coms, only : &
            mcats          & ! intent(in)
           ,mivap          & ! intent(in)
           ,mix02          & ! intent(in)
           ,alphasfc       & ! intent(in)
           ,mcat1          & ! intent(in)
           ,mcat2          & ! intent(in)
           ,mcat33         ! ! intent(in)

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer                                    , intent(in)  :: m1,ngr,if_adap,mynum
   real                                       , intent(in)  :: dtlt, dtlti
   real(kind=8)                               , intent(in)  :: time
   real        , dimension(m1)                , intent(in)  :: zm
   real        , dimension(m1)                , intent(in)  :: dzt ! Not used
   real                                       , intent(in)  :: rtgt
   real                                       , intent(out) :: pcpg, qpcpg, dpcpg
   !----- Variables needed for Harrington radiation scheme --------------------------------!
   real, dimension(m1,maxkfall,nembfall,nhcat), intent(in)  :: pcpfillc, pcpfillr
   real, dimension(maxkfall,nembfall,nhcat)   , intent(in)  :: sfcpcp
   real, dimension(m1)                        , intent(in)  :: zt
   !----- Local Variables: ----------------------------------------------------------------!
   integer :: k,jflag,lcat,icv,icx,mc1,mc2,mc3,mc4,mcat
   integer :: lhcat,j1,j2
   !---------------------------------------------------------------------------------------!

   !----- Initialisation of some thermodynamic properties ---------------------------------!
   call thrmstr(m1)
   call each_column(m1,dtlt)

   !----- Diagnose hydrometeor mean mass emb, and if necessary, number concentration. -----!
   jflag = 1
   do lcat = 1,7
      if (availcat(lcat) .and. k2(lcat) >= k1(lcat)) then
         call enemb(lcat,jflag)
      end if
   end do

   !----- Set up matrix for heat/vapor diffusion computation. -----------------------------!
   do lcat = 1,7
      if (availcat(lcat) .and. k2(lcat) >= k1(lcat)) then
         call diffprep(lcat)
      end if
   end do

   !----- Implicit matrix solution of atmospheric vapor density. (1st guess. ) ------------!
   call vapdiff()

   !----- Vapor flux applied to each category.  Do not change the order of these ----------!
   do icv = 1,7
      lcat = mivap(icv)
      if (availcat(lcat) .and. k2(lcat) >= k1(lcat)) call vapflux(lcat)
   end do


   !----- Conversion between pristine ice and snow due to vapor flux ----------------------!
   if (jnmb(4) >= 1) then
      k1(3) = min(k1(3),k1(4))
      k2(3) = max(k2(3),k2(4))
      k1(4) = k1(3)
      k2(4) = k2(3)
      if (k2(3) >= k1(3)) call psxfer()
   end if

   !----- Diagnose hydrometeor mean mass emb, and if necessary, number concentration. -----!
   jflag = 2
   do lcat = 1,7
      if (availcat(lcat) .and. k2(lcat) >= k1(lcat)) then
         call enemb(lcat,jflag)
      end if
   end do

   !----- Adjust temperature --------------------------------------------------------------!
   call newtemp()

   !----- Autoconversion from cloud droplets to rain drops --------------------------------!
   if (availcat(2) .and. k2(1) >= k1(1)) call auto_accret(dtlt)

   !----- This is ???? --------------------------------------------------------------------!
   call effxy(m1)

   !----- Self collection of rain, aggregates, graupel, hail:  number change only. --------!
   selfloop: do lcat = 2, 7
     if (lcat == 3 .or. lcat == 4) cycle selfloop

     mc1 = mcats(lcat)
     if (progncat(lcat) .and. k2(lcat) >= k1(lcat)) call cols(lcat,mc1)
   end do selfloop

   !----- Self collection of pristine ice, snow -------------------------------------------!
   do lcat = 3,4
      mc1 = mcat33(lcat)
      if (availcat(lcat) .and. availcat(mc1) .and. k2(lcat) >= k1(lcat)) then
         call col3344(lcat,5,mc1)
      end if
   end do

   !----- Collection between pristine ice and snow ----------------------------------------!
   if (availcat(5) .and. k2(3) >= k1(3)) call col3443 (3,4,5)

   !----- Ice-ice collisions --------------------------------------------------------------!
   do icx = 1,9
      mc1 = mcat1(icx,1)
      mc2 = mcat1(icx,2)
      mc3 = mcat1(icx,3)
      mc4 = mcat1(icx,4)
      j1  = max(k1(mc1),k1(mc2))
      j2  = min(k2(mc1),k2(mc2))
      if (availcat(mc1) .and. availcat(mc3) .and. j2 >= j1) then
         call col1 (mc1,mc2,mc3,mc4,j1,j2)
      end if
   end do

   !----- Ice-cloud collisions ------------------------------------------------------------!
   do lcat = 4,7
      mc1 = mcat2(lcat,1)
      mc2 = mcat2(lcat,2)
      j1 = max(k1(1),k1(lcat))
      j2 = min(k2(1),k2(lcat))
      if (availcat(lcat) .and. availcat(mc1) .and. j2 >= j1) then
         call col2 (1,lcat,mc1,mc2,j1,j2,dtlt)
      endif
   enddo

   !----- Ice-rain collisions -------------------------------------------------------------!
   do lcat = 3,7
      j1 = max(k1(2),k1(lcat))
      j2 = min(k2(2),k2(lcat))
      if (availcat(lcat) .and. availcat(7) .and. j2 >= j1) then
         call col3(2,lcat,7,j1,j2)
      endif
   enddo

  !----- ????? ----------------------------------------------------------------------------!
  call colxfers()

   !----- ????? ---------------------------------------------------------------------------!
   do mcat = 1,7
      lcat = mix02(mcat)
      if (availcat(lcat)) call x02(lcat)
   end do


   !----- Compute cloud nucleation --------------------------------------------------------!
   if (availcat(1)) call cldnuc(m1)

   !----- ????? ---------------------------------------------------------------------------!
   k1(1) = min(k1(1),k1cnuc)
   k2(1) = max(k2(1),k2cnuc)
   k3(1) = max(k2(1),k3(1))
   jflag = 1
   if (availcat(1)) call enemb(1,jflag)

   !----- Compute ice nucleation ----------------------------------------------------------!
   if (availcat(3)) call icenuc(m1,ngr,dtlt)

   k1(3) = min(k1(3),k1pnuc)
   k2(3) = max(k2(3),k2pnuc)
   k3(3) = max(k2(3),k3(3))

   !----- Do not change the order of these ------------------------------------------------!
   jflag = 1
   if (availcat(3)) call pc03(3,jflag)
   if (availcat(1)) call pc03(1,jflag)   

   !----- Zero out precip arrays. ---------------------------------------------------------!
   pcpg       = 0.
   qpcpg      = 0.
   dpcpg      = 0.
   accpx(:)   = 0.
   pcprx(:)   = 0.

   !----- Fill ch1 array for sedim in this column -----------------------------------------!
   do lhcat = 2,nhcat
      ch1(lhcat) = dtlt * cfvt(lhcat) / rtgt
   enddo

   !----- Zero out array for accumulating sedim changes to thil ---------------------------!
   dsed_thil(:) = 0.

   do lcat = 2,7
      if (availcat(lcat) .and. k2(lcat) >= k1(lcat)) then
         !------ Compute sedimentation for all 6 precipitating categories -----------------!
         call sedim (m1,lcat,if_adap,mynum,pcpg,qpcpg,dpcpg,dtlti,pcpfillc,pcpfillr,sfcpcp &
                    ,dzt)
      end if
   enddo
   !----- Dsed_thil is the tendency on theta_il due to sedimentation ----------------------!
   do k = lpw,m1
      thil(k) = thil(k) + dsed_thil(k)
   end do

  return
end subroutine mcphys_main
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine copyback(m1,i,j,thp,btheta,rtp,dn0,micro)

   use mem_micro, only: &
           micro_vars   ! ! INTENT(IN) ! Only a type structure

   use micphys, only:   &
           k1           & ! intent(in)
          ,k2           & ! intent(in)
          ,k3           & ! intent(in)
          ,lpw          & ! intent(in)
          ,ncat         & ! intent(in)
          ,jnmb         & ! intent(in)
          ,availcat     & ! intent(in)
          ,progncat     & ! intent(in)
          ,rx           & ! intent(in)
          ,cx           & ! intent(in)
          ,qx           & ! intent(in)
          ,pottemp      & ! intent(in)
          ,thil         & ! intent(in)
          ,rtot         & ! intent(in)
          ,rhoa         & ! intent(in)
          ,accpx        & ! intent(in)
          ,pcprx        ! ! intent(in)

   use mem_scratch, only : &
           vctr11        ! ! intent(out)
   
   use rconstants, only: t00,cliqt3,cliq,cice,alli
   use therm_lib , only: qtk

   implicit none

   !----- Arguments: ----------------------------------------------------------------------!
   integer                         , intent(in)    :: m1,i,j
   real             , dimension(m1), intent(inout) :: thp,rtp,dn0,btheta
   type (micro_vars)               , intent(inout) :: micro
   !----- Local variables -----------------------------------------------------------------!
   integer                                         :: k,lcat
   real                                            :: tcoal, fracliq
   !---------------------------------------------------------------------------------------!

  
   do k =lpw,m1-1
      thp(k)    = thil(k)
      !btheta(k) = pottemp(k)
      !dn0(k)    = rhoa(k)
      rtp(k)    = rtot(k)
   end do

   !----- 1. Cloud ------------------------------------------------------------------------!
   if (availcat(1)) then
      do k=lpw,k3(1)
         micro%rcp(k,i,j) = rx(k,1)
         if (progncat(1)) micro%ccp(k,i,j) = cx(k,1)
      end do
   end if
   !---------------------------------------------------------------------------------------!



   !----- 2. Rain -------------------------------------------------------------------------!
   if (availcat(2)) then
      micro%accpr(i,j) = micro%accpr(i,j) + accpx(2)
      micro%pcprr(i,j) = pcprx(2)
      do k=lpw,k2(10)
         micro%q2(k,i,j)  = qx(k,2)
         micro%rrp(k,i,j) = rx(k,2)
         if (progncat(2)) micro%crp(k,i,j) = cx(k,2)
      end do
   end if
   !---------------------------------------------------------------------------------------!



   !----- 3. Pristine ice -----------------------------------------------------------------!
   if (availcat(3)) then
      micro%accpp(i,j) = micro%accpp(i,j) + accpx(3)
      micro%pcprp(i,j) = pcprx(3)
      do k=lpw,k3(3)
         micro%rpp(k,i,j) = rx(k,3)
         micro%cpp(k,i,j) = cx(k,3)
      end do
   end if
   !---------------------------------------------------------------------------------------!



   !----- 4. Snow -------------------------------------------------------------------------!
   if (availcat(4)) then
      micro%accps(i,j) = micro%accps(i,j) + accpx(4)
      micro%pcprs(i,j) = pcprx(4)
      do k=lpw,k2(10)
         micro%rsp(k,i,j) = rx(k,4)
         if (progncat(4)) micro%csp(k,i,j) = cx(k,4)
      end do
   end if
   !---------------------------------------------------------------------------------------!



   !----- 5. Aggregates -------------------------------------------------------------------!
   if (availcat(5)) then
      micro%accpa(i,j) = micro%accpa(i,j) + accpx(5)
      micro%pcpra(i,j) = pcprx(5)
      do k=lpw,k2(10)
         micro%rap(k,i,j) = rx(k,5)
         if (progncat(5)) micro%cap(k,i,j) = cx(k,5)
      end do
   end if
   !---------------------------------------------------------------------------------------!



   !----- 6. Graupel ----------------------------------------------------------------------!
   if (availcat(6)) then
      micro%accpg(i,j) = micro%accpg(i,j) + accpx(6)
      micro%pcprg(i,j) = pcprx(6)
      do k=lpw,k2(10)
         micro%q6(k,i,j)  = qx(k,6)
         micro%rgp(k,i,j) = rx(k,6)
         if (progncat(6)) micro%cgp(k,i,j) = cx(k,6)
      end do
   end if
   !---------------------------------------------------------------------------------------!



   !----- 7. Hail -------------------------------------------------------------------------!
   if (availcat(7)) then
      micro%accph(i,j) = micro%accph(i,j) + accpx(7)
      micro%pcprh(i,j) = pcprx(7)
      do k=lpw,k2(10)
         micro%q7(k,i,j)  = qx(k,7)
         micro%rhp(k,i,j) = rx(k,7)
         if (progncat(7)) micro%chp(k,i,j) = cx(k,7)
      end do
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine copyback
!==========================================================================================!
!==========================================================================================!
