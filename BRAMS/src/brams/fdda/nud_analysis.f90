!===================================== Change Log =========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!     This subroutine will perform the data assimilation during the run time.              !
!------------------------------------------------------------------------------------------!
subroutine datassim()

   use mem_tend    , only : tend      ! ! intent(inout)
   use mem_basic   , only : basic_g   & ! intent(in)
                          , co2_on    ! ! intent(in)
   use mem_grid    , only : ibcon     ! ! intent(in)
   use mem_varinit , only : varinit_g ! ! intent(in)
   use mem_scratch , only : scratch   ! ! intent(in)
   use node_mod    , only : ia        & ! intent(in)
                          , iz        & ! intent(in)
                          , ja        & ! intent(in)
                          , jz        & ! intent(in)
                          , mxp       & ! intent(in)
                          , myp       & ! intent(in)
                          , mzp       ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer :: iwest
   integer :: ieast
   integer :: jsouth
   integer :: jnorth
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Set bounds for nudging this sub-domain.  We first assume the sub-domain is in the !
   ! middle, so we assign the node solvable boundaries.                                    !
   !---------------------------------------------------------------------------------------!
   iwest  = ia
   ieast  = iz
   jsouth = ja
   jnorth = jz
   !----- We now check whether this node subdomain includes edges or corners. -------------!
   if (iand(ibcon,1) /= 0)  iwest = 1   !----- Western edge. ------------------------------!
   if (iand(ibcon,2) /= 0)  ieast = mxp !----- Eastern edge. ------------------------------!
   if(jdim == 1) then
      if (iand(ibcon,4) /= 0) jsouth = 1   !----- Southern edge. --------------------------!
      if (iand(ibcon,8) /= 0) jnorth = myp !----- Northern edge. --------------------------!
   end if
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      We now must copy CO2 to scratch variables so it won't fail in case CO2 is not    !
   ! solved in this run.                                                                   !
   !---------------------------------------------------------------------------------------!
   call azero4(mzp*mxp*myp,scratch%vt3do,scratch%vt3df,scratch%vt3dp,scratch%vt3de)
   if (co2_on) then
      call atob(mzp*mxp*myp,basic_g(ngrid)%co2p   ,scratch%vt3do)
      call atob(mzp*mxp*myp,varinit_g(ngrid)%varop,scratch%vt3dp)
      call atob(mzp*mxp*myp,varinit_g(ngrid)%varof,scratch%vt3df)
      call atob(mzp*mxp*myp,tend%co2t             ,scratch%vt3de)
   end if
   !---------------------------------------------------------------------------------------!

   !----- Now we call the basic boundary and analysis nudging scheme. ---------------------!
   call nudge(mzp, mxp, myp, iwest, ieast, jsouth, jnorth        , varinit_g(ngrid)%varwts &
             , varinit_g(ngrid)%varup  , varinit_g(ngrid)%varvp  , varinit_g(ngrid)%varpp  &
             , varinit_g(ngrid)%vartp  , varinit_g(ngrid)%varrp  , scratch%vt3dp           &
             , varinit_g(ngrid)%varuf  , varinit_g(ngrid)%varvf  , varinit_g(ngrid)%varpf  &
             , varinit_g(ngrid)%vartf  , varinit_g(ngrid)%varrf  , scratch%vt3df           &
             , basic_g(ngrid)%up       , basic_g(ngrid)%vp       , basic_g(ngrid)%theta    &
             , basic_g(ngrid)%rtp      , basic_g(ngrid)%pp       , scratch%vt3do           &
             , tend%ut                 , tend%vt                 , tend%tht                &
             , tend%rtt                , tend%pt                 , scratch%vt3de           )

   !---------------------------------------------------------------------------------------!
   !     If CO2 is solved, then we copy back the tendency.                                 !
   !---------------------------------------------------------------------------------------!
   if (co2_on) then
      call atob(mzp*mxp*myp,scratch%vt3de,tend%co2t)
   end if
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Now we check whether we want to use the condensate nudging scheme.                !
   !---------------------------------------------------------------------------------------!
   if (nud_cond == 1 .and. time >= tcond_beg .and. time <= tcond_end) then
   call nudge_cond(mzp,mxp,myp,iwest,ieast,jsouth,jnorth                                   &
                  ,varinit_g(ngrid)%varwts,varinit_g(ngrid)%varrph,varinit_g(ngrid)%varcph &
                  ,varinit_g(ngrid)%varrfh,varinit_g(ngrid)%varcfh,basic_g(ngrid)%rtp      &
                  ,tend%rtt               )
   end if
   return
end subroutine datassim
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will compute the tendency due to nudging towards the provided data.  !
!------------------------------------------------------------------------------------------! 
subroutine nudge(m1,m2,m3,ia,iz,ja,jz,varwts,varup,varvp,varpp,vartp,varrp,varop           &
                ,varuf,varvf,varpf,vartf,varrf,varof,up,vp,theta,rtp,pp,co2p               &
                ,ut,vt,tht,rtt,pt,co2t)

   use mem_varinit , only : nud_type      & ! intent(in)
                          , wt_nudge_grid & ! intent(in)
                          , wt_nudge_uv   & ! intent(in)
                          , wt_nudge_th   & ! intent(in)
                          , wt_nudge_rt   & ! intent(in)
                          , wt_nudge_co2  ! ! intent(in)
   use mem_scratch , only : vctr1         & ! intent(inout)
                          , vctr2         & ! intent(inout)
                          , vctr3         & ! intent(inout)
                          , vctr4         & ! intent(inout)
                          , vctr5         & ! intent(inout)
                          , vctr6         & ! intent(inout)
                          , vctr11        & ! intent(inout)
                          , vctr12        & ! intent(inout)
                          , vctr13        & ! intent(inout)
                          , vctr14        & ! intent(inout)
                          , vctr15        & ! intent(inout)
                          , vctr16        ! ! intent(inout)
   use mem_grid    , only : ngrid         & ! intent(in)
                          , jdim          & ! intent(in)
                          , time          ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: m1
   integer                     , intent(in)    :: m2
   integer                     , intent(in)    :: m3
   integer                     , intent(in)    :: ia
   integer                     , intent(in)    :: iz
   integer                     , intent(in)    :: ja
   integer                     , intent(in)    :: jz
   real   , dimension(m1,m2,m3), intent(in)    :: varup
   real   , dimension(m1,m2,m3), intent(in)    :: varvp
   real   , dimension(m1,m2,m3), intent(in)    :: vartp
   real   , dimension(m1,m2,m3), intent(in)    :: varrp
   real   , dimension(m1,m2,m3), intent(in)    :: varpp
   real   , dimension(m1,m2,m3), intent(in)    :: varop
   real   , dimension(m1,m2,m3), intent(in)    :: varuf
   real   , dimension(m1,m2,m3), intent(in)    :: varvf
   real   , dimension(m1,m2,m3), intent(in)    :: vartf
   real   , dimension(m1,m2,m3), intent(in)    :: varrf
   real   , dimension(m1,m2,m3), intent(in)    :: varpf
   real   , dimension(m1,m2,m3), intent(in)    :: varof
   real   , dimension(m1,m2,m3), intent(in)    :: varwts
   real   , dimension(m1,m2,m3), intent(in)    :: up
   real   , dimension(m1,m2,m3), intent(in)    :: vp
   real   , dimension(m1,m2,m3), intent(in)    :: theta
   real   , dimension(m1,m2,m3), intent(in)    :: rtp
   real   , dimension(m1,m2,m3), intent(in)    :: pp
   real   , dimension(m1,m2,m3), intent(inout) :: ut
   real   , dimension(m1,m2,m3), intent(inout) :: vt
   real   , dimension(m1,m2,m3), intent(inout) :: tht
   real   , dimension(m1,m2,m3), intent(inout) :: rtt
   real   , dimension(m1,m2,m3), intent(inout) :: pt             
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: i
   integer                                     :: j
   integer                                     :: k
   real                                        :: tfact
   real                                        :: wt_uv
   real                                        :: wt_th
   real                                        :: wt_pi
   real                                        :: wt_rt
   real                                        :: wt_co2
   !---------------------------------------------------------------------------------------!

   !----- Compute the linear time factor. -------------------------------------------------!
   select case (nud_type)
   case (1)
      tfact=sngl((time-htime1)/(htime2-htime1))
   case (2)
      tfact=sngl((time-vtime1)/(vtime2-vtime1))
   end select

   wt_uv  = wt_nudge_grid(ngrid) * wt_nudge_uv
   wt_th  = wt_nudge_grid(ngrid) * wt_nudge_th
   wt_pi  = wt_nudge_grid(ngrid) * wt_nudge_pi
   wt_rt  = wt_nudge_grid(ngrid) * wt_nudge_rt
   wt_co2 = wt_nudge_grid(ngrid) * wt_nudge_co2

   do j = ja, jz
      do i = ia, iz
         do k = 1, m1
            vctr1(k)  = varup(k,i,j)   + (varuf(k,i,j)-varup(k,i,j)) * tfact
            vctr2(k)  = varvp(k,i,j)   + (varvf(k,i,j)-varvp(k,i,j)) * tfact
            vctr3(k)  = vartp(k,i,j)   + (vartf(k,i,j)-vartp(k,i,j)) * tfact
            vctr4(k)  = varrp(k,i,j)   + (varrf(k,i,j)-varrp(k,i,j)) * tfact
            vctr5(k)  = varpp(k,i,j)   + (varpf(k,i,j)-varpp(k,i,j)) * tfact
            vctr6(k)  = varop(k,i,j)   + (varof(k,i,j)-varop(k,i,j)) * tfact

            vctr11(k) = (varwts(k,i,j) + varwts(k,min(m2,i+1),j))    *.5 * wt_uv
            vctr12(k) = (varwts(k,i,j) + varwts(k,i,min(m3,j+jdim))) *.5 * wt_uv
            vctr13(k) =  varwts(k,i,j) * wt_th
            vctr14(k) =  varwts(k,i,j) * wt_rt
            vctr15(k) =  varwts(k,i,j) * wt_pi
            vctr16(k) =  varwts(k,i,j) * wt_co2
         end do

         do k=1,m1
            ut(k,i,j)   = ut(k,i,j)   + vctr11(k) * (vctr1(k) - up   (k,i,j))
            vt(k,i,j)   = vt(k,i,j)   + vctr12(k) * (vctr2(k) - vp   (k,i,j))
            tht(k,i,j)  = tht(k,i,j)  + vctr13(k) * (vctr3(k) - theta(k,i,j))
            rtt(k,i,j)  = rtt(k,i,j)  + vctr14(k) * (vctr4(k) - rtp  (k,i,j))
            pt(k,i,j)   = pt(k,i,j)   + vctr15(k) * (vctr5(k) - pp   (k,i,j))
            co2t(k,i,j) = co2t(k,i,j) + vctr16(k) * (vctr6(k) - co2p (k,i,j))
         end do

      end do
   end do

   return
end subroutine nudge
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will compute the tendency due to nudging towards the provided data   !
! for some condensed water.                                                                !
!------------------------------------------------------------------------------------------! 
subroutine nudge_cond(m1,m2,m3,ia,iz,ja,jz,varwts,varrph,varcph,varrfh,varcfh,rtp,rtt)
   use mem_grid   , only : ngrid     & ! intent(in)
                         , time      ! ! intent(in)
   use mem_varinit, only : condtime1 & ! intent(in)
                         , condtime2 ! ! intent(in)
   use mem_scratch, only : vctr7     & ! intent(inout)
                         , vctr8     & ! intent(inout)
                         , vctr9     ! ! intent(inout)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in)    :: m1
   integer                  , intent(in)    :: m2
   integer                  , intent(in)    :: m3
   integer                  , intent(in)    :: ia
   integer                  , intent(in)    :: iz
   integer                  , intent(in)    :: ja
   integer                  , intent(in)    :: jz
   real, dimension(m1,m2,m3), intent(in)    :: varrph
   real, dimension(m1,m2,m3), intent(in)    :: varcph
   real, dimension(m1,m2,m3), intent(in)    :: varrfh
   real, dimension(m1,m2,m3), intent(in)    :: varcfh
   real, dimension(m1,m2,m3), intent(in)    :: varwts
   real, dimension(m1,m2,m3), intent(in)    :: rtp
   real, dimension(m1,m2,m3), intent(inout) :: rtt             
   !----- Local variables. ----------------------------------------------------------------!
   integer                                  :: i
   integer                                  :: j
   integer                                  :: k
   real                                     :: tfact
   real                                     :: wt_rc
   !---------------------------------------------------------------------------------------!

   !----- Tfact is temporal interpolation weight. -----------------------------------------!
   tfact = sngl((time-condtime1)/(condtime2-condtime1)) * wt_nudgec_grid(ngrid)/t_nudge_rc

   !----- Wt_rc is timescale and grid-dependent weight. -----------------------------------!
   wt_rc = wt_nudgec_grid(ngrid)/t_nudge_rc

   do j=ja,jz
      do i=ia,iz
         do k=1,m1
            vctr7(k)  = varrph(k,i,j) + (varrfh(k,i,j)-varrph(k,i,j)) * tfact
            vctr8(k)  = varcph(k,i,j) + (varcfh(k,i,j)-varcph(k,i,j)) * tfact
            vctr17(k) = varwts(k,i,j) * wt_rc
         end do

         !------ Only nudging total water where condensate exists... ----------------------!
         do k=1,m1
            if (vctr8(k) > 0.) rtt(k,i,j) = rtt(k,i,j) + vctr17(k) * (vctr7(k)-rtp(k,i,j))
         end do
      end do
   end do

   return
end subroutine nudge_cond
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the nudging weight for large scale and model tendencies.    !
!------------------------------------------------------------------------------------------!
subroutine varweight(n1,n2,n3,varwts,topt,rtgt)
   use mem_grid   , only : nzp      ! ! intent(in)
   use mem_varinit, only : tnudcent & ! intent(in)
                         , tnudlat  & ! intent(in)
                         , tnudtop  & ! intent(in)
                         , znudtop  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: n1,n2,n3
   real   , dimension(n1,n2,n3), intent(inout) :: varwts
   real   , dimension(n1,n2)   , intent(in)    :: topt
   real   , dimension(n1,n2)   , intent(in)    :: rtgt
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: i
   integer                                     :: j
   integer                                     :: k
   real                                        :: tnudcenti
   real                                        :: tnudtopi
   real                                        :: tnudlati
   real                                        :: rown
   real                                        :: rows
   real                                        :: rowe
   real                                        :: roww
   real                                        :: zloc
   real                                        :: wttop
   real                                        :: wtlat
   real                                        :: delzi
   !---------------------------------------------------------------------------------------!

   !----- Nothing to be done here if nudging is turned off... -----------------------------!
   if (nudlat <= 0) return


   !---------------------------------------------------------------------------------------!
   !      The weight depends on the inverse of the time scale.  In case the time scale is  !
   ! tiny or zero, we understand that the user doesn't want that kind of nudging.          !
   !---------------------------------------------------------------------------------------!
   !----- Domain centre -------------------------------------------------------------------!
   if (tnudcent > .01) then
      tnudcenti = 1./tnudcent
   else
      tnudcenti = 0.
   end if
   !----- Domain top ----------------------------------------------------------------------!
   if (tnudtop > .01) then
      tnudtopi = 1./tnudtop
   else
      tnudtopi = 0.
   end if
   !----- Domain lateral boundaries. ------------------------------------------------------!
   if (tnudlat > .01) then
      tnudlati = 1./tnudlat
   else
      tnudlati = 0.
   end if
   !---------------------------------------------------------------------------------------!

   if (ztop > znudtop) then
      delzi = 1. / (ztop - znudtop)
   elseif (tnudtop > .01) then
      write (unit=*,fmt='(a)')           'ZNUDTOP can''t be above the model top (ZTOP)...'
      write (unit=*,fmt='(a,1x,es12.5)') ' ZNUDTOP = ',znudtop
      write (unit=*,fmt='(a,1x,es12.5)') '    ZTOP = ',ztop
      call abort_run('Incorrect specification of znudtop!','varweight','nud_analysis.f90')
   end if


   do j=1,n3
      do i=1,n2

         !----- Quadratic weight function for lateral boundaries. -------------------------!
         rown = max(0.,float(j+nudlat-n3))
         rows = max(0.,float(nudlat+1-j))
         rowe = max(0.,float(i+nudlat-n2))
         roww = max(0.,float(nudlat+1-i))
         wtlat = max(rown*rown,rows*rows,rowe*rowe,roww*roww)/real(nudlat*nudlat)

         !----- Linear weight function for top boundary. ----------------------------------!
         do k=1,nzp
            zloc  = zt(k) * rtgt(i,j) + topt(i,j)
            wttop = max(0.,(zloc - znudtop) * delzi)
            !----- Full 3-D weight function. ----------------------------------------------!
            varwts(k,i,j)= tnudcenti + max(tnudlati * wtlat, tnudtopi * wttop)
         end do
      end do
   end do

   return
end subroutine varweight
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subrotuine will interpolate the atmospheric variables between the previous and   !
! the future fields (since meteorological drivers usually have a coarser time resolution   !
! than BRAMS runs).                                                                        !
!------------------------------------------------------------------------------------------!
subroutine vfintrpf(ifm,ifflag)
   use mem_grid   , only : grid_g    & ! intent(in)
                         , nnxp      & ! intent(in)
                         , nnyp      & ! intent(in)
                         , nnzp      & ! intent(in)
                         , maxnxp    & ! intent(in)
                         , maxnyp    ! ! intent(in)
   use mem_scratch, only : scratch   ! ! intent(inout)
   use mem_varinit, only : varinit_g ! ! intent(inout)
   use mem_basic  , only : co2_on    & ! intent(in)
                         , basic_g   ! ! intent(inout)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: ifm
   integer                     , intent(in)    :: ifflag
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: icm
   !---------------------------------------------------------------------------------------!

   !----- Find the parent grid. -----------------------------------------------------------!
   icm = nxtnest(ifm)
   if (icm == 0) return

   !----- Copy topography to vt2da. -------------------------------------------------------!
   call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1,scratch%scr1,grid_g(icm)%topt)
   call eintp(scratch%scr1,scratch%scr2,1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
   call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1,scratch%scr2,scratch%vt2da)

   select case (ifflag)
   case (1)
      !----- Interpolate varwts. ----------------------------------------------------------!
      call fmint4(varinit_g(icm)%varwts,varinit_g(ifm)%varwts,basic_g(icm)%dn0             &
                 ,basic_g(ifm)%dn0,scratch%vt2da,ifm,icm,'t',0)

   case (2)
      !----- Interpolate future level atmospheric variables. ------------------------------!
      call fmint4(varinit_g(icm)%varuf,varinit_g(ifm)%varuf,basic_g(icm)%dn0u              &
                 ,basic_g(ifm)%dn0u,scratch%vt2da,ifm,icm,'u',1)
      call fmint4(varinit_g(icm)%varvf,varinit_g(ifm)%varvf,basic_g(icm)%dn0v              &
                 ,basic_g(ifm)%dn0v,scratch%vt2da,ifm,icm,'v',1)
      call fmint4(varinit_g(icm)%varpf,varinit_g(ifm)%varpf,basic_g(icm)%dn0v              &
                 ,basic_g(ifm)%dn0v,scratch%vt2da,ifm,icm,'t',1)
      call fmint4(varinit_g(icm)%vartf,varinit_g(ifm)%vartf,basic_g(icm)%dn0               &
                 ,basic_g(ifm)%dn0,scratch%vt2da,ifm,icm,'t',1)
      call fmint4(varinit_g(icm)%varrf,varinit_g(ifm)%varrf,basic_g(icm)%dn0               &
                 ,basic_g(ifm)%dn0,scratch%vt2da,ifm,icm,'t',1)
      if (co2_on) then
         call fmint4(varinit_g(icm)%varof,varinit_g(ifm)%varof,basic_g(icm)%dn0            &
                    ,basic_g(ifm)%dn0,scratch%vt2da,ifm,icm,'t',1)
      end if
   end select

   return
end subroutine vfintrpf
!==========================================================================================!
!==========================================================================================!
