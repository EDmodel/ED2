subroutine sw_twostream_clump(salb,   &
     scosz,  &
     scosaoi, &
     sLAIm,   &
     ncoh,  &
     pft,  &
     LAI_in,  &
     PAR_beam_flip,  &
     PAR_diffuse_flip,  &
     SW_abs_beam_flip,   &
     SW_abs_diffuse_flip,   &
     DW_vislo_beam,  &
     DW_vislo_diffuse,  &
     UW_vishi_beam,  &
     UW_vishi_diffuse,  &
     DW_nirlo_beam,  &
     DW_nirlo_diffuse,  &
     UW_nirhi_beam,   &
     UW_nirhi_diffuse)

  use pft_coms, only: clumping_factor, n_pft, phenology
  use canopy_radiation_coms, only: diffuse_backscatter_nir,   &
       diffuse_backscatter_vis, &
       leaf_scatter_nir, leaf_scatter_vis,   &
       visible_fraction_dir, visible_fraction_dif

  implicit none

  integer :: il,ncoh,nfcoh,ipft,ncoh2,iband,i,j
  integer, dimension(2*ncoh) :: indx
  real :: salb,scosz,scosaoi,sLAIm,srshort,srshortd,UW_nirhi_beam,UW_nirhi_diffuse
  real :: DW_vislo_diffuse,UW_vishi_beam,DW_nirlo_beam,DW_nirlo_diffuse
  real :: DW_vislo_beam, UW_vishi_diffuse
  integer, dimension(ncoh) :: pft
  real, dimension(ncoh+1) :: SW_abs_beam_flip,PAR_beam_flip
  real, dimension(ncoh+1) :: SW_abs_diffuse_flip,PAR_diffuse_flip
  real(kind=8) :: alb,cosz,cosaoi,LAIm,rshort,rshortd,lambda,lambda_tot,LAI_reduction
  real(kind=8) :: beam_backscatter,eta,zeta,raddiff,diffuse_band
  real(kind=8), dimension(n_pft) :: leaf_scatter
  real(kind=8) :: exk,exki,zetai
  real(kind=8) :: d,rhoo,sigma,source_bot,source_top
  real(kind=8), dimension(n_pft) :: diffuse_backscatter
  real(kind=8) :: cumulative_lai
  
  real(kind=8), dimension(ncoh) :: LAI_in
  real(kind=8), dimension(ncoh) :: expkl_top,expkl_bot,expamk_top,expamk_bot
  real(kind=8), dimension(ncoh) :: expapk_top,expapk_bot,A_top,A_bot,B_top
  real(kind=8), dimension(ncoh) :: B_bot,C_top,C_bot,F_top,F_bot,G_top,G_bot
  real(kind=8), dimension(ncoh) :: H_top,H_bot,beam_bot
  real(kind=8), dimension(ncoh+1) :: upward_vis_beam,upward_vis_diffuse
  real(kind=8), dimension(ncoh+1) :: upward_nir_beam, upward_nir_diffuse
  real(kind=8), dimension(ncoh+1) :: downward_nir_beam, downward_nir_diffuse
  real(kind=8), dimension(ncoh+1) :: downward_vis_beam,downward_vis_diffuse
  real(kind=8), dimension(2*ncoh) :: mastervec_beam,masveccp_beam
  real(kind=8), dimension(2*ncoh) :: mastervec_diffuse,masveccp_diffuse
  real(kind=8), dimension(2*ncoh,2) :: matal
  real(kind=8), dimension(2*ncoh,5) :: mastermat
  real(kind=8), dimension(2*ncoh,2*ncoh) :: masmatcp  
  integer :: ind
  
  ! Convert input variable to double precision
  alb = dble(salb)
  cosz = max(0.03d0,dble(scosz))
  cosaoi = max(0.03d0,dble(scosaoi))
  LAIm = dble(sLAIm)

  ! Calculate factors common for NIR, PAR
  ncoh2 = 2*ncoh
  lambda = 0.5/cosaoi
  lambda_tot = 0.0
  nfcoh = 0
  do il=1,ncoh
     nfcoh = nfcoh + 1
     ipft = pft(il)
     lambda_tot = lambda_tot + clumping_factor(ipft)
  enddo
  
  lambda_tot = lambda_tot * lambda / float(nfcoh)
  LAI_reduction = dble(min(1.0,sLAIm))
  beam_backscatter = (0.5 + cosz) * (1.0 - cosz*log(1.0+1.0/cosz))
 
  ! Loop over bands
  do iband = 1,2
     if(iband.eq.1)then
        !  First, visible (or PAR)
        do ipft = 1,n_pft
           leaf_scatter(ipft) = dble(leaf_scatter_vis(ipft))
           diffuse_backscatter(ipft) = dble(diffuse_backscatter_vis(ipft))
        enddo
!        raddiff = (rshort -rshortd) * visible_fraction_dir
!        diffuse_band = rshortd * visible_fraction_dif
     elseif(iband.eq.2)then
        !  Then, near infrared (or NIR)
        do ipft = 1,n_pft
           leaf_scatter(ipft) = dble(leaf_scatter_nir)
           diffuse_backscatter(ipft) = dble(diffuse_backscatter_nir)
        enddo
!        raddiff = (rshort -rshortd) * (1.0-visible_fraction_dir)  
!        diffuse_band = rshortd * (1.0-visible_fraction_dif) 
     endif
     
     ! Calculate more factors for this band
     
     ! Calculate the forcings
     beam_bot(ncoh) = exp(-lambda*clumping_factor(pft(ncoh))*LAI_in(ncoh))
     do il=ncoh-1,1,-1
        beam_bot(il) = beam_bot(il+1)  &
             *exp(-lambda*clumping_factor(pft(il))*LAI_in(il))
     enddo
     
     do il=1,ncoh
        ipft = pft(il)
        eta = (1.0 - (1.0-diffuse_backscatter(ipft)) * leaf_scatter(ipft))  &
             * clumping_factor(ipft)
        zeta = leaf_scatter(ipft) * diffuse_backscatter(ipft) *   &
             clumping_factor(ipft)
        exk = sqrt(eta**2 - zeta**2)
        exki = 1.0/exk
        zetai = 1.0/zeta
        ! sources
        source_bot = clumping_factor(ipft)*lambda*leaf_scatter(ipft)  &
             * beam_backscatter * beam_bot(il)

        source_top = source_bot   &
             * exp(lambda*clumping_factor(ipft)*LAI_in(il))
        ! forcing coefficients
        rhoo = - (zeta + eta + clumping_factor(ipft)*lambda)   &
             * clumping_factor(ipft) * lambda    &
             * leaf_scatter(ipft) * beam_backscatter   &
             * beam_bot(il)
        sigma = clumping_factor(ipft) * lambda
        ! calculate exponentials only once
        expkl_bot(il)=1.0
        expkl_top(il)=exp(exk*LAI_in(il))
        expamk_bot(il) = 1.0
        expamk_top(il) = exp((sigma-exk)*LAI_in(il))
        expapk_bot(il) = 1.0
        expapk_top(il) = exp((sigma+exk)*LAI_in(il))
        A_bot(il) = -source_bot*zetai
        A_top(il) = -source_top*zetai &
             +0.5*zetai*(eta*exki-1.0)*expkl_top(il)*rhoo/(sigma-exk) &
             *(expamk_top(il)-expamk_bot(il))  &
             -0.5*zetai*(eta*exki+1.0)/expkl_top(il)*rhoo/(sigma+exk) &
             *(expapk_top(il)-expapk_bot(il))
        B_bot(il) = 0.5*zetai*(eta*exki-1.0)
        B_top(il) = 0.5*zetai*(eta*exki-1.0)*expkl_top(il)
        C_bot(il) = -0.5*zetai*(eta*exki+1.0)
        C_top(il) = -0.5*zetai*(eta*exki+1.0)/expkl_top(il)
        F_bot(il) = 0.0
        F_top(il) = 0.5*exki*expkl_top(il)*rhoo/(sigma-exk)  &
             *(expamk_top(il)-expamk_bot(il))  &
             -0.5*exki/expkl_top(il)*rhoo/(sigma+exk)  &
             *(expapk_top(il)-expapk_bot(il))
        G_bot(il) = 0.5*exki
        G_top(il) = 0.5*exki*expkl_top(il)
        H_bot(il) = -0.5*exki
        H_top(il) = -0.5*exki/expkl_top(il)
     enddo

     
     ! Organize the matrix coefficients
     
     do j=1,ncoh2
        do i=1,ncoh2
           masmatcp(i,j)=0.0
        enddo
        mastervec_beam(j)=0.0
        masveccp_beam(j)=0.0
        mastervec_diffuse(j)=0.0
        masveccp_diffuse(j)=0.0
     enddo
     
     masmatcp(1,1)=G_top(ncoh)
     masmatcp(1,2)=H_top(ncoh)
     mastervec_beam(1)=-F_top(ncoh)
     mastervec_diffuse(1)=1.0
     masveccp_beam(1)=mastervec_beam(1)
     masveccp_diffuse(1)=mastervec_diffuse(1)

     do i=2,ncoh2-2,2
        masmatcp(i,i-1)=G_bot(nint((ncoh2-i+2)*0.5))
        masmatcp(i,i)=H_bot(nint((ncoh2-i+2)*0.5))
        masmatcp(i,i+1)=-G_top(nint((ncoh2-i)*0.5))
        masmatcp(i,i+2)=-H_top(nint((ncoh2-i)*0.5))
        mastervec_beam(i)=-F_bot(nint((ncoh2-i+2)*0.5))+F_top(nint((ncoh2-i)*0.5))
        mastervec_diffuse(i)=0.0
        masveccp_beam(i)=mastervec_beam(i)
        masveccp_diffuse(i)=mastervec_diffuse(i)
     enddo
     do i=3,ncoh2-1,2
        masmatcp(i,i-2)=B_bot(nint((ncoh2-i+3)*0.5))
        masmatcp(i,i-1)=C_bot(nint((ncoh2-i+3)*0.5))
        masmatcp(i,i)=-B_top(nint((ncoh2-i+1)*0.5))
        masmatcp(i,i+1)=-C_top(nint((ncoh2-i+1)*0.5))
        mastervec_beam(i)=-A_bot(nint((ncoh2-i+3)*0.5))+A_top(nint((ncoh2-i+1)*0.5))
        masveccp_beam(i)=mastervec_beam(i)
        mastervec_diffuse(i)=0.0
        masveccp_diffuse(i)=mastervec_diffuse(i)
     enddo
     masmatcp(ncoh2,ncoh2-1)=B_bot(1)-alb*G_bot(1)
     masmatcp(ncoh2,ncoh2)=C_bot(1)-alb*H_bot(1)
     mastervec_beam(ncoh2)= -A_bot(1)+alb*beam_bot(1)
     masveccp_beam(ncoh2)= mastervec_beam(ncoh2)
     mastervec_diffuse(ncoh2)= 0.0
     masveccp_diffuse(ncoh2)= mastervec_diffuse(ncoh2)
     
     ! Prep for inversion
     
     mastermat(1,1)=0.
     mastermat(1,2)=0.
     mastermat(1,3)=masmatcp(1,1)
     mastermat(1,4)=masmatcp(1,2)
     mastermat(1,5)=0.
     do i=2,ncoh2-2,2
        mastermat(i,1)=0.
        mastermat(i,2)=masmatcp(i,i-1)
        mastermat(i,3)=masmatcp(i,i)
        mastermat(i,4)=masmatcp(i,i+1)
        mastermat(i,5)=masmatcp(i,i+2)
     enddo
     do i=3,ncoh2-1,2
        mastermat(i,1)=masmatcp(i,i-2)
        mastermat(i,2)=masmatcp(i,i-1)
        mastermat(i,3)=masmatcp(i,i)
        mastermat(i,4)=masmatcp(i,i+1)
        mastermat(i,5)=0.
     enddo
     mastermat(ncoh2,1)=0.
     mastermat(ncoh2,2)=masmatcp(ncoh2,ncoh2-1)
     mastermat(ncoh2,3)=masmatcp(ncoh2,ncoh2)
     mastermat(ncoh2,4)=0.
     mastermat(ncoh2,5)=0.
     
     ! Invert the matrix
     call bandec(mastermat,ncoh2,2,2,matal,indx,d)

     ! Backsubstitute for beam and diffuse
     call banbks(mastermat,ncoh2,2,2,matal,indx,mastervec_beam)
     call banbks(mastermat,ncoh2,2,2,matal,indx,mastervec_diffuse)
     
     ! Improve the solution
     call mprove(masmatcp,mastermat,matal,ncoh2,5,2,indx,masveccp_beam,mastervec_beam)
     call mprove(masmatcp,mastermat,matal,ncoh2,5,2,indx,masveccp_diffuse,mastervec_diffuse)
     

     if(iband.eq.1)then
        do i=3,ncoh2-1,2
           ind = nint((ncoh2-i+1)*0.5)+1
           upward_vis_beam(ind) = A_bot(ind) + B_bot(ind) *  &
                mastervec_beam(i-2) + C_bot(ind) * mastervec_beam(i-1)
           upward_vis_diffuse(ind) = B_bot(ind) * mastervec_diffuse(i-2)   &
                + C_bot(ind) * mastervec_diffuse(i-1)
        enddo
        do i=2,ncoh2-2,2
           ind = nint((ncoh2-i)*0.5)+1
           downward_vis_beam(ind) = beam_bot(ind) + F_bot(ind)  &
                + H_bot(ind) * mastervec_beam(i)  &
                + G_bot(ind) * mastervec_beam(i-1)
           downward_vis_diffuse(ind) = H_bot(ind) * mastervec_diffuse(i)  &
                + G_bot(ind) * mastervec_diffuse(i-1)
        enddo
        upward_vis_beam(ncoh+1) = B_top(ncoh) * mastervec_beam(1)  &
             + C_top(ncoh) * mastervec_beam(2) + A_top(ncoh)
        upward_vis_diffuse(ncoh+1) = B_top(ncoh) * mastervec_diffuse(1)  &
             + C_top(ncoh) * mastervec_diffuse(2)
        downward_vis_beam(ncoh+1) = 1.0
        downward_vis_diffuse(ncoh+1) = 1.0
        downward_vis_beam(1) = G_bot(1) * mastervec_beam(ncoh2-1) +   &
             H_bot(1) * mastervec_beam(ncoh2) + F_bot(1) + beam_bot(1)
        downward_vis_diffuse(1) = G_bot(1) * mastervec_diffuse(ncoh2-1) +  &
             H_bot(1) * mastervec_diffuse(ncoh2)
        upward_vis_beam(1) = alb * downward_vis_beam(1)
        upward_vis_diffuse(1) = alb * downward_vis_diffuse(1)
     elseif(iband.eq.2)then
        do i=3,ncoh2-1,2
           ind = nint((ncoh2-i+1)*0.5)+1
           upward_nir_beam(ind) = A_bot(ind) + B_bot(ind) *   &
                mastervec_beam(i-2) + C_bot(ind) * mastervec_beam(i-1)
           upward_nir_diffuse(ind) = B_bot(ind) * mastervec_diffuse(i-2) +  &
                C_bot(ind) * mastervec_diffuse(i-1)
        enddo
        do i=2,ncoh2-2,2
           ind = nint((ncoh2-i)*0.5)+1
           downward_nir_beam(ind) = beam_bot(ind) + F_bot(ind) +   &
                H_bot(ind) * mastervec_beam(i) + G_bot(ind) *   &
                mastervec_beam(i-1)
           downward_nir_diffuse(ind) = H_bot(ind) * mastervec_diffuse(i)  +  &
                G_bot(ind) * mastervec_diffuse(i-1)
        enddo
        upward_nir_beam(ncoh+1) = B_top(ncoh) * mastervec_beam(1) +  &
             C_top(ncoh) * mastervec_beam(2) + A_top(ncoh)
        upward_nir_diffuse(ncoh+1) = B_top(ncoh) * mastervec_diffuse(1) +  &
             C_top(ncoh) * mastervec_diffuse(2)
        downward_nir_beam(ncoh+1) = 1.0
        downward_nir_diffuse(ncoh+1) = 1.0
        downward_nir_beam(1) = G_bot(1) * mastervec_beam(ncoh2-1) +  &
             H_bot(1) * mastervec_beam(ncoh2) + F_bot(1) + beam_bot(1)
        downward_nir_diffuse(1) = G_bot(1) * mastervec_diffuse(ncoh2-1) +  &
             H_bot(1) * mastervec_diffuse(ncoh2)
        upward_nir_beam(1) = alb*downward_nir_beam(1)
        upward_nir_diffuse(1) = alb*downward_nir_diffuse(1)
     endif
  enddo
  
  do il=1,ncoh
     PAR_beam_flip(il) = visible_fraction_dir * sngl(    &
          downward_vis_beam(il+1) - downward_vis_beam(il)  +   &
          upward_vis_beam(il) - upward_vis_beam(il+1)  )
     PAR_diffuse_flip(il) = visible_fraction_dif * sngl(    &
          downward_vis_diffuse(il+1) - downward_vis_diffuse(il)  +   &
          upward_vis_diffuse(il) - upward_vis_diffuse(il+1)  )
     SW_abs_beam_flip(il) = PAR_beam_flip(il) +   &
          (1.0 - visible_fraction_dir) * sngl(  &
          downward_nir_beam(il+1)-downward_nir_beam(il)  &
          +upward_nir_beam(il)-upward_nir_beam(il+1))
     SW_abs_diffuse_flip(il) = PAR_diffuse_flip(il) +   &
          (1.0 - visible_fraction_dif) * sngl(  &
          downward_nir_diffuse(il+1)-downward_nir_diffuse(il)  &
          +upward_nir_diffuse(il)-upward_nir_diffuse(il+1))
     if(PAR_beam_flip(il).lt.0.0)PAR_beam_flip(il)=0.0
     if(PAR_diffuse_flip(il).lt.0.0)PAR_diffuse_flip(il)=0.0
     if(SW_abs_beam_flip(il).lt.0.0)SW_abs_beam_flip(il)=0.0
     if(SW_abs_diffuse_flip(il).lt.0.0)SW_abs_diffuse_flip(il)=0.0
  enddo
  
  DW_vislo_beam = max(0.0,sngl(downward_vis_beam(1))) * visible_fraction_dir
  DW_vislo_diffuse = max(0.0,sngl(downward_vis_diffuse(1))) * visible_fraction_dif
  UW_vishi_beam = max(0.0,sngl(upward_vis_beam(ncoh+1))) * visible_fraction_dir
  UW_vishi_diffuse = max(0.0,sngl(upward_vis_diffuse(ncoh+1))) *visible_fraction_dif
  DW_nirlo_beam = max(0.0,sngl(downward_nir_beam(1))) * (1.0-visible_fraction_dir)
  DW_nirlo_diffuse = max(0.0,sngl(downward_nir_diffuse(1))) * (1.0-visible_fraction_dif)
  UW_nirhi_beam = max(0.0,sngl(upward_nir_beam(ncoh+1))) * (1.0-visible_fraction_dir)
  UW_nirhi_diffuse = max(0.0,sngl(upward_nir_diffuse(ncoh+1))) * (1.0-visible_fraction_dif)
  
  return
end subroutine sw_twostream_clump

!==============================================================
subroutine bandec(a,n,m1,m2,al,indx,d)
implicit none

real(kind=8), parameter :: tiny=1.0e-20
integer :: n,m1,m2
real(kind=8), dimension(n,m1+m2+1) :: a
real(kind=8), dimension(n,m1) :: al
integer, dimension(n) :: indx
real(kind=8) :: d,tvar
integer :: i,j,k,l,mm
real(kind=8) :: dum

!------------------

mm=m1+m2+1
l=m1
do i=1,m1
   do j=m1+2-i,mm
      a(i,j-l)=a(i,j)
   enddo
   l=l-1
   do j=mm-l,mm
      a(i,j)=0.
   enddo
enddo
d=1.0
l=m1
do k=1,n
   dum=a(k,1)
   i=k
   if(l.lt.n)l=l+1
   do j=k+1,l
      if(abs(a(j,1)).gt.abs(dum))then
         dum=a(j,1)
         i=j
      endif
   enddo
   indx(k)=i
   if(dum.eq.0.0)a(k,1)=tiny
   if(i.ne.k)then
      d=-d
      do j=1,mm
         tvar=a(k,j)
         a(k,j)=a(i,j)
         a(i,j)=tvar
      enddo
   endif
   do i=k+1,l
      dum=a(i,1)/a(k,1)
      al(k,i-k)=dum
      do j=2,mm
         a(i,j-1)=a(i,j)-dum*a(k,j)
      enddo
      a(i,mm)=0.0
   enddo
enddo
return
end subroutine bandec


!==========================================================================
subroutine banbks(a,n,m1,m2,al,indx,b)
implicit none

integer :: n,m1,m2,i,k,l,mm
real(kind=8), dimension(n,m1+m2+1) :: a
real(kind=8), dimension(n,m1) :: al
integer, dimension(n) :: indx
real(kind=8), dimension(n) :: b
real(kind=8) :: dum,tvar
!-------------------------

mm=m1+m2+1
l=m1
do k=1,n
   i=indx(k)
   if(i.ne.k)then
      tvar=b(k)
      b(k)=b(i)
      b(i)=tvar
   endif
   if(l.lt.n)l=l+1
   do i=k+1,l
      b(i)=b(i)-al(k,i-k)*b(k)
   enddo
enddo
l=1




do i=n,1,-1
   dum=b(i)
   do k=2,l
      dum=dum-a(i,k)*b(k+i-1)
   enddo
   b(i)=dum/a(i,1)
   if(l.lt.mm)l=l+1
enddo



return
end subroutine banbks

!=======================================================================
subroutine mprove(a,alud,matal,n,np,npp,indx,b,x)
implicit none

integer, parameter :: nmax=100
integer :: n,np,i,j,npp
real(kind=8), dimension(n,n) :: a
real(kind=8), dimension(n,np) :: alud
real(kind=8), dimension(n,npp) :: matal
integer, dimension(n) :: indx
real(kind=8), dimension(n) :: b,x
real(kind=8), dimension(n) :: r
real(kind=8) :: sdp

!-------------------------------------
do i=1,n
   sdp=-b(i)
   do j=1,n
      sdp=sdp+a(i,j)*x(j)
   enddo
   r(i)=sdp
enddo

call banbks(alud,n,npp,npp,matal,indx,r)
do i=1,n
   x(i)=x(i)-r(i)
enddo
return
end subroutine mprove

!=====================================================================
subroutine lw_twostream(ncoh, semgs, sT_grnd, pft, LAI, T_veg,   &
     lw_v_surf, lw_v_incid, downward_lw_below_surf, downward_lw_below_incid,  &
     upward_lw_below_surf, upward_lw_below_incid, upward_lw_above_surf,  &
     upward_lw_above_incid)

  use canopy_radiation_coms, only: emis_v, mubar
  use consts_coms, only: stefan

  implicit none
  
  integer :: ncoh
  real :: semgs
  real :: st_grnd
  integer, dimension(ncoh) :: pft
  real(kind=8), dimension(ncoh) :: LAI
  real(kind=8), dimension(ncoh) :: T_veg
  real, dimension(ncoh) :: lw_v_surf
  real, dimension(ncoh) :: lw_v_incid
  real :: downward_lw_below_surf
  real :: upward_lw_above_surf
  real :: upward_lw_below_surf
  real :: downward_lw_below_incid
  real :: upward_lw_above_incid
  real ::  upward_lw_below_incid
  real(kind=8) :: emgs
  real(kind=8) :: T_grnd
  integer :: ncoh2
  integer :: il
  real(kind=8) :: zeta
  real(kind=8) :: eta
  real(kind=8) :: exk
  real(kind=8) :: zetai
  real(kind=8) :: exki
  real(kind=8), dimension(ncoh) :: source
  real(kind=8), dimension(ncoh) :: forcing
  real(kind=8), dimension(ncoh+1) :: explai
  real(kind=8), dimension(ncoh+1) :: exmlai
  real(kind=8), dimension(ncoh) :: A_dw
  real(kind=8), dimension(ncoh) :: B_dw
  real(kind=8), dimension(ncoh) :: C_dw
  real(kind=8), dimension(ncoh) :: D_dw
  real(kind=8), dimension(ncoh) :: A_uw
  real(kind=8), dimension(ncoh) :: B_uw
  real(kind=8), dimension(ncoh) :: C_uw
  real(kind=8), dimension(ncoh) :: D_uw
  real(kind=8), dimension(ncoh) :: E_uw
  real(kind=8), dimension(ncoh) :: F_uw
  integer :: i
  integer :: j
  real(kind=8), dimension(2*ncoh,2*ncoh) :: masmatcp  
  real(kind=8), dimension(2*ncoh) :: mastervec_surf
  real(kind=8), dimension(2*ncoh) :: mastervec_incid
  real(kind=8), dimension(ncoh+1) :: DW_incid
  real(kind=8), dimension(ncoh+1) :: DW_surf
  real(kind=8), dimension(ncoh+1) :: UW_incid
  real(kind=8), dimension(ncoh+1) :: UW_surf
  real(kind=8), dimension(2*ncoh,5) :: mastermat
  real(kind=8), dimension(2*ncoh,2) :: matal
  integer :: ind
  integer, dimension(2*ncoh) :: indx
  real(kind=8) :: d

  emgs = dble(semgs)
  t_grnd = dble(st_grnd)

  ncoh2 = 2*ncoh

  do il=1,ncoh
     zeta=2.0*(1.0-emis_v(pft(il)))/(3.0*mubar)
     eta=(2.0+emis_v(pft(il)))/(3.0*mubar)
     exk=sqrt(eta**2-zeta**2)
     exki=1.0/exk
     zetai=1.0/zeta
     source(il) = emis_v(pft(il))*stefan*T_veg(il)**4
     forcing(il) = -(zeta+eta)*source(il)
     explai(il)=exp(exk*LAI(il))
     exmlai(il)=exp(-exk*LAI(il))

     ! coefficient of lambda1 (and minus the coefficient of lambda2) for the bottom of a layer, downwelling radiation.
     A_dw(il) = 0.5 * exki

     ! coefficient of lambda1, top of layer, downwelling radiation
     B_dw(il) = 0.5*exki*explai(il)

     ! coefficient of lambda2, top of layer, downwelling radiation
     C_dw(il) = -0.5*exki*exmlai(il)

     ! term of downwelling radiation not multiplying a lambda
     D_dw(il) = 0.5*(exki**2)*forcing(il) * (explai(il) + exmlai(il) - 2.0)

     A_uw(il) = 0.5*zetai*(eta*exki-1.0)
     B_uw(il) = -0.5*zetai*(eta*exki+1.0)
     C_uw(il) =  -source(il)*zetai
     D_uw(il) = A_uw(il) * explai(il)
     E_uw(il) = B_uw(il) * exmlai(il)
     F_uw(il) = -source(il)*zetai  &
          +0.5*zetai*(eta*exki-1.0)*explai(il)  &
          *(forcing(il)*exki*(1.0           &
          -exmlai(il)))  &
          -0.5*zetai*(eta*exki+1.0)*exmlai(il)  &
          *(forcing(il)*exki*(explai(il)-1.0))
  enddo

!  print*,"ncoh:",ncoh
!  print*,"A_dw",A_dw(1:ncoh)
!  print*,"B_dw",B_dw(1:ncoh)
!  print*,"C_dw",C_dw(1:ncoh)
!  print*,"D_dw",D_dw(1:ncoh)!

!  print*,"A_uw",A_uw(1:ncoh)
!  print*,"B_uw",B_uw(1:ncoh)
!  print*,"C_uw",C_uw(1:ncoh)
 ! print*,"D_uw",D_uw(1:ncoh)


  do j=1,ncoh2
     do i=1,ncoh2
        masmatcp(i,j)=0.0
     enddo
     mastervec_surf(j)=0.0
     mastervec_incid(j)=0.0
  enddo

  ! Vector is of the form: (lambda_N, lambda_{N-1},...,lambda_1)

  masmatcp(1,1)=B_dw(ncoh)
  masmatcp(1,2)=C_dw(ncoh)
  mastervec_surf(1)=-D_dw(ncoh)
  mastervec_incid(1)=1.0

  do i=2,ncoh2-2,2
     ind = nint((ncoh2-i)*0.5)
     masmatcp(i,i-1)=-A_dw(ind+1)
     masmatcp(i,i)=A_dw(ind+1)
     masmatcp(i,i+1)=B_dw(ind)
     masmatcp(i,i+2)=C_dw(ind)
     mastervec_surf(i)=-D_dw(ind)
     mastervec_incid(i)=0.0
  enddo
  do i=3,ncoh2-1,2
     ind = nint((ncoh2-i+1)*0.5)
     masmatcp(i,i-2)=-A_uw(ind+1)
     masmatcp(i,i-1)=-B_uw(ind+1)
     masmatcp(i,i)=D_uw(ind)
     masmatcp(i,i+1)=E_uw(ind)
     mastervec_surf(i)=C_uw(ind+1)-F_uw(ind)
     mastervec_incid(i)=0.0
  enddo
  masmatcp(ncoh2,ncoh2-1)=A_uw(1)-(1.0-emgs)*A_dw(1)
  masmatcp(ncoh2,ncoh2)=B_uw(1)+(1.0-emgs)*A_dw(1)
  mastervec_surf(ncoh2)=emgs*stefan*T_grnd**4 - C_uw(1)
  mastervec_incid(ncoh2)=0.0

  mastermat(1,1)=0.
  mastermat(1,2)=0.
  mastermat(1,3)=masmatcp(1,1)
  mastermat(1,4)=masmatcp(1,2)
  mastermat(1,5)=0.
  do i=2,ncoh2-2,2
     mastermat(i,1)=0.
     mastermat(i,2)=masmatcp(i,i-1)
     mastermat(i,3)=masmatcp(i,i)
     mastermat(i,4)=masmatcp(i,i+1)
     mastermat(i,5)=masmatcp(i,i+2)
  enddo
  do i=3,ncoh2-1,2
     mastermat(i,1)=masmatcp(i,i-2)
     mastermat(i,2)=masmatcp(i,i-1)
     mastermat(i,3)=masmatcp(i,i)
     mastermat(i,4)=masmatcp(i,i+1)
     mastermat(i,5)=0.
  enddo
  mastermat(ncoh2,1)=0.
  mastermat(ncoh2,2)=masmatcp(ncoh2,ncoh2-1)
  mastermat(ncoh2,3)=masmatcp(ncoh2,ncoh2)
  mastermat(ncoh2,4)=0.
  mastermat(ncoh2,5)=0.
  
  ! Invert matrix
  call bandec(mastermat,ncoh2,2,2,matal,indx,d)

  ! Backsubstitute for contributions of ground and vegetation
  call banbks(mastermat,ncoh2,2,2,matal,indx,mastervec_surf)

  ! Backsubstitute for contribution of incident longwave at canopy top
  call banbks(mastermat,ncoh2,2,2,matal,indx,mastervec_incid)

  do i=3,ncoh2-1,2
     ind = nint((ncoh2-i+1)*0.5)
     UW_surf(ind+1)=masmatcp(i,i)*mastervec_surf(i)   &
          +masmatcp(i,i+1)*mastervec_surf(i+1)+F_uw(ind)
     UW_incid(ind+1)=masmatcp(i,i)*mastervec_incid(i)   &
          +masmatcp(i,i+1)*mastervec_incid(i+1)
  enddo
  do i=2,ncoh2-2,2
     ind = nint((ncoh2-i)*0.5)
     DW_surf(ind+1)=masmatcp(i,i+1)*mastervec_surf(i+1)   &
          +masmatcp(i,i+2)*mastervec_surf(i+2)+D_dw(ind)
     DW_incid(ind+1)=masmatcp(i,i+1)*mastervec_incid(i+1)   &
          +masmatcp(i,i+2)*mastervec_incid(i+2)
  enddo

  UW_surf(ncoh+1)=D_uw(ncoh)*mastervec_surf(1)+E_uw(ncoh)*mastervec_surf(2)  &
       + F_uw(ncoh)
  UW_incid(ncoh+1)=D_uw(ncoh)*mastervec_incid(1)+E_uw(ncoh)*mastervec_incid(2)
  DW_surf(ncoh+1)=0.0
  DW_incid(ncoh+1)=1.0
  DW_surf(1)=A_dw(1)*(mastervec_surf(ncoh2-1)-mastervec_surf(ncoh2))
  DW_incid(1)=A_dw(1)*(mastervec_incid(ncoh2-1)-mastervec_incid(ncoh2))
  UW_surf(1)=(1.0-emgs)*DW_surf(1)+emgs*stefan*T_grnd**4
  UW_incid(1)=(1.0-emgs)*DW_incid(1)

  do il = 1,ncoh
     lw_v_surf(il) = sngl(DW_surf(il+1) - DW_surf(il) + UW_surf(il) -   &
          UW_surf(il+1))
     lw_v_incid(il) = sngl(DW_incid(il+1) - DW_incid(il) + UW_incid(il) -   &
          UW_incid(il+1))
  enddo

  downward_lw_below_surf = sngl(DW_surf(1))
  downward_lw_below_incid = sngl(DW_incid(1))
  upward_lw_below_surf = sngl(UW_surf(1))
  upward_lw_below_incid = sngl(UW_incid(1))
  upward_lw_above_surf = sngl(UW_surf(ncoh+1))
  upward_lw_above_incid = sngl(UW_incid(ncoh+1))

  return
end subroutine lw_twostream

!================================================================
subroutine sw_beers_law(ncoh,rshort,SW_abs_flip,PAR_flip,DW_vislo,DW_nirlo  &
     ,UW_vishi,UW_nirhi,lai_in)

  use canopy_radiation_coms, only: visible_fraction
  implicit none
  integer :: ncoh,il
  real :: rshort,DW_vislo,DW_nirlo,UW_vishi,UW_nirhi,beam
  real, dimension(ncoh) :: SW_abs_flip,PAR_flip,lai_in
  real, dimension(ncoh+1) :: downward
  
  !=========================================================
  UW_vishi = 0.0
  UW_nirhi = 0.0

  beam = 1.0 * rshort
  downward(ncoh+1) = beam
  do il=ncoh,1,-1
     downward(il) = downward(il+1)*exp(-0.5*lai_in(il))
     SW_abs_flip(il) = downward(il+1) - downward(il)
     PAR_flip(il) = visible_fraction * SW_abs_flip(il)
  enddo
  DW_vislo = downward(1)
  DW_nirlo = 0.0


  return
end subroutine sw_beers_law


!================================================================
subroutine sw_beers_law_crown(ncoh,rshort,SW_abs_flip,PAR_flip,DW_vislo,DW_nirlo  &
     ,UW_vishi,UW_nirhi,lai_in,pft,CA)
  use canopy_radiation_coms, only: visible_fraction
  use pft_coms, only: clumping_factor
  implicit none

  integer :: ncoh,il
  real :: rshort,DW_vislo,DW_nirlo,UW_vishi,UW_nirhi,beam
  real, dimension(ncoh) :: SW_abs_flip,PAR_flip,lai_in,CA
  integer, dimension(ncoh) :: pft
  real, dimension(ncoh+1) :: downward
  
  !=========================================================
  UW_vishi = 0.0
  UW_nirhi = 0.0

  beam = 1.0 * rshort
  downward(ncoh+1) = beam
  do il=ncoh,1,-1
     downward(il) = downward(il+1)*(1-CA(il)*(1-exp(-0.5*clumping_factor(pft(il))*lai_in(il)/CA(il))))
     SW_abs_flip(il) = downward(il+1) - downward(il)
     PAR_flip(il) = visible_fraction * SW_abs_flip(il)
  enddo
  DW_vislo = downward(1)
  DW_nirlo = 0.0


  return
end subroutine sw_beers_law_crown
