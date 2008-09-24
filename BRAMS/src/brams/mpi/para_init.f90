!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine node_decomp(init)

  use mem_grid
  use var_tables
  use mem_scratch
  use mem_basic
  use rpara

  implicit none

  integer :: init

  integer :: ngr,idn,isn,ng,mp_nzp,numbuff,icm,nsiz,jnode,ncols,ifm &
       ,nestvar,nc,nf,npvar2,npvar3,nv,num_lbc_buff,num_nest_buff &
       ,num_feed_buff,itype,i1,i2,j1,j2,ixy,ixyz,memf

  logical :: failed

  ! try to read domain decomposition data from input file

  call PAR_decomp_input (maxmach, maxgrds, nmachs, ngrids, nnxp, nnyp, &
       ixb, ixe, iyb, iye, failed)


  !      Decompose all grids into subdomains

  do ngr = 1,ngrids


     ! escape if domain decomposition specified by input file

     if (failed) then

        ! Obtain estimates of the fraction of computational time (work) required
        ! for each column in the region of the domain.
        
        
        call PAR_est_time(nnxp(ngr),nnyp(ngr),scratch%scr1(1)  &
             ,basic_g(ngr)%cputime(1,1),init)
        
        ! Decompose the grid taking into account the work numbers.
        
        nsiz = nnxp(ngr)+nnyp(ngr)
        
        call PAR_decomp(nnxp(ngr),nnyp(ngr),nsiz,nmachs,scratch%scr1(1)  &
             ,scratch%scr2(1),scratch%scr2(1+nsiz),scratch%scr2(1+2*nsiz)  &
             ,scratch%scr2(1+3*nsiz),scratch%scr2(1+4*nsiz)  &
             ,scratch%scr2(1+5*nsiz),scratch%scr2(1+6*nsiz)  &
             ,ixb(1,ngr),ixe(1,ngr),iyb(1,ngr),iye(1,ngr))
     end if
        
  enddo

  ! Compute various bounds for the subdomains

  call PAR_decomp_bounds(ngrids,nnxp,nnyp,1,1)

  ! Determine node sending paths and numbers of receiving nodes

  call PAR_node_paths(maxgrds,ngrids,nxpmax,nypmax,nnxp,nnyp  &
       ,nxtnest,nmachs,machnum,ibcflg,nxbeg,nxend,nybeg,nyend  &
       ,ipm,jpm,ixb,ixe,iyb,iye  &
       ,ibounds                  &        ! For reproducibility - Saulo Barros
       ,inode_paths_master,iget_paths_master,2,ibnd,jbnd)


  !     Compute  send and receive buffer sizes. These will be maximum of
  !       long timestep, turbulence, nest boundaries, and nest feedback.
  !       Small timestep will use same buffers as they are always smaller.

  do idn=1,nmachs
     do isn=1,nmachs
        lbc_buffs(1,idn,isn)=0
        lbc_buffs(2,idn,isn)=0
     enddo
  enddo


  mp_nzp=0
  do ng=1,ngrids
     mp_nzp=max(mp_nzp,nnzp(ng))
  enddo

  do ng=1,ngrids

     !          Find number of nested variables to be communicated.
     icm=ng
     ifm=ng
     if(ng /= 1) icm=nxtnest(ifm)
     nestvar=4
     do nf=1,num_scalar(ifm)
        do nc=1,num_scalar(icm)
           if(scalar_tab(nf,ifm)%name==scalar_tab(nc,icm)%name)  &
                nestvar=nestvar+1
        enddo
     enddo

     !  Find number of lbc variables to be communicated.
     npvar3=0 ; npvar2=0
     do nv = 1,num_var(ng)
        if(vtab_r(nv,ng)%impt1 == 1 ) then
           if (vtab_r(nv,ng)%idim_type==2) npvar2=npvar2+1 
           if (vtab_r(nv,ng)%idim_type==3) npvar3=npvar3+1 
        endif
     enddo


     do isn=1,nmachs
        do idn=1,nmachs
           num_lbc_buff=0
           num_nest_buff=0
           num_feed_buff=0

           itype=1
           i1=inode_paths_master(1,itype,ng,idn,isn)
           i2=inode_paths_master(2,itype,ng,idn,isn)
           j1=inode_paths_master(3,itype,ng,idn,isn)
           j2=inode_paths_master(4,itype,ng,idn,isn)
           if(i1.ne.0) then
              ixy=(i2-i1+1)*(j2-j1+1)
              ixyz=(i2-i1+1)*(j2-j1+1)*(mp_nzp)
              num_lbc_buff=ixyz*npvar3+ixy*npvar2  &
                   +2*(npvar3+npvar2+100)
           endif

           itype=5
           i1=inode_paths_master(1,itype,ng,idn,isn)
           i2=inode_paths_master(2,itype,ng,idn,isn)
           j1=inode_paths_master(3,itype,ng,idn,isn)
           j2=inode_paths_master(4,itype,ng,idn,isn)
           if(i1.ne.0) then
              ixyz=(i2-i1+1)*(j2-j1+1)*(mp_nzp)
              num_nest_buff=ixyz*nestvar+2*(nestvar+100)
           endif

           !  itype=6
           itype=7                  !For reproducibility - Saulo Barros
           i1=inode_paths_master(1,itype,ng,idn,isn)
           i2=inode_paths_master(2,itype,ng,idn,isn)
           j1=inode_paths_master(3,itype,ng,idn,isn)
           j2=inode_paths_master(4,itype,ng,idn,isn)
           if(i1.ne.0) then
              ixyz=(i2-i1+1)*(j2-j1+1)*(mp_nzp)
              num_feed_buff=ixyz*nestvar+2*(nestvar+100)
           endif
           lbc_buffs(1,idn,isn)= max(lbc_buffs(1,idn,isn)  &
                ,num_lbc_buff,num_nest_buff,num_feed_buff)
        enddo
     enddo

     do isn=1,nmachs
        do idn=1,nmachs
           num_lbc_buff=0
           num_nest_buff=0
           num_feed_buff=0

           itype=1
           i1=inode_paths_master(1,itype,ng,isn,idn)
           i2=inode_paths_master(2,itype,ng,isn,idn)
           j1=inode_paths_master(3,itype,ng,isn,idn)
           j2=inode_paths_master(4,itype,ng,isn,idn)
           if(i1.ne.0) then
              ixy=(i2-i1+1)*(j2-j1+1)
              ixyz=(i2-i1+1)*(j2-j1+1)*(mp_nzp)
              num_lbc_buff=ixyz*npvar3+ixy*npvar2  &
                   +2*(npvar3+npvar2+100)
           endif

           itype=5
           i1=inode_paths_master(1,itype,ng,isn,idn)
           i2=inode_paths_master(2,itype,ng,isn,idn)
           j1=inode_paths_master(3,itype,ng,isn,idn)
           j2=inode_paths_master(4,itype,ng,isn,idn)
           if(i1.ne.0) then
              ixyz=(i2-i1+1)*(j2-j1+1)*(mp_nzp)
              num_nest_buff=ixyz*nestvar+2*(nestvar+100)
           endif

           !  itype=6
           itype=7                  !For reproducibility - Saulo Barros
           i1=inode_paths_master(1,itype,ng,isn,idn)
           i2=inode_paths_master(2,itype,ng,isn,idn)
           j1=inode_paths_master(3,itype,ng,isn,idn)
           j2=inode_paths_master(4,itype,ng,isn,idn)
           if(i1.ne.0) then
              ixyz=(i2-i1+1)*(j2-j1+1)*(mp_nzp)
              num_feed_buff=ixyz*nestvar+2*(nestvar+100)
           endif
           lbc_buffs(2,idn,isn)= max(lbc_buffs(2,idn,isn)  &
                ,num_lbc_buff,num_nest_buff,num_feed_buff)
        enddo
     enddo
  enddo

  !       Check nest boundary receive buffer size

  itype=5

  do idn=1,nmachs
     newbuff_nest1(idn)=1
     nbuff_nest1(idn)=0

     do ng=1,ngrids
        numbuff=0
        icm=nxtnest(ng)
        do isn=1,nmachs
           i1=inode_paths_master(1,itype,ng,idn,isn)
           i2=inode_paths_master(2,itype,ng,idn,isn)
           j1=inode_paths_master(3,itype,ng,idn,isn)
           j2=inode_paths_master(4,itype,ng,idn,isn)
           memf=(i2-i1+1)*(j2-j1+1)*(mp_nzp)*nestvar
           numbuff=numbuff+memf
        enddo
        nbuff_nest1(idn)=max(nbuff_nest1(idn),numbuff)
     enddo

  enddo


  return
end subroutine node_decomp

!     ****************************************************************

subroutine onenode()

  use mem_grid
  use node_mod

  implicit none

  mmzp=nnzp
  mmxp=nnxp
  mmyp=nnyp

  mxp=nnxp(ngrid)
  myp=nnyp(ngrid)
  mzp=nnzp(ngrid)

  ia=2
  iz=mxp-1
  ia_1=max(ia-1,1)
  ia_2=max(ia-2,1)
  ia_3=max(ia-3,1)
  ia1=ia+1
  ia2=ia+2
  ia3=ia+3
  iz_1=iz-1
  iz_2=iz-2
  iz_3=iz-3
  iz1=min(iz+1,mxp)
  iz2=min(iz+2,mxp)
  iz3=min(iz+3,mxp)
  izu=iz-1

  if(nnyp(ngrid).gt.1) then
     ja=2
     jz=myp-1
     ja_1=max(ja-1,1)
     ja_2=max(ja-2,1)
     ja_3=max(ja-3,1)
     ja1=ja+1
     ja2=ja+2
     ja3=ja+3
     jz_1=jz-1
     jz_2=jz-2
     jz_3=jz-3
     jz1=min(jz+1,myp)
     jz2=min(jz+2,myp)
     jz3=min(jz+3,myp)
     jzv=jz-1
  else
     ja=1
     jz=1
     ja_1=1
     ja_2=1
     ja_3=1
     ja1=1
     ja2=1
     ja3=1
     jz_1=1
     jz_2=1
     jz_3=1
     jz1=1
     jz2=1
     jz3=1
     jzv=1
  endif

  i0=0
  j0=0
  ibcon=1+2+4+8
  ipara=0

  mi0(1:ngrids)=i0
  mj0(1:ngrids)=j0
  mia(1:ngrids)=ia
  miz(1:ngrids)=mmxp(1:ngrids)-1
  mja(1:ngrids)=ja
  mjz(1:ngrids)=mmyp(1:ngrids)-1



  return
end subroutine onenode
!
!     ****************************************************************
!
subroutine PAR_decomp_bounds(ngrids,nnxp,nnyp,nbndx,nbndy)

  use rpara

  implicit none

  integer :: ngrids,nnxp(*),nnyp(*),nbndx,nbndy

  integer :: ng,nm,nx,ny

  !        Compute various subdomain boundary numbers for the nodes
  !             nxbeg,nybeg,nxend,nyend - portions of full domain that nodes will have
  !                                     - includes overlap region
  !             ixoff,iyoff  - subdomain offsets relative to full domain
  !             nxbegc,nxendc,nybegc,nyendc - subdomain "compute" points,
  !                         or normal thermodynamic tendency points (2-nx, 2-ny for
  !                          non-parallel run
  !             ibcflag - flag denoting if real boundary is on subdomain
  !                       bit 1=west, bit 2=east, bit 3=south, bit 4=north

  do ng=1,ngrids

     nx=nnxp(ng)
     ny=nnyp(ng)

     do nm=1,nmachs

        ibcflg(nm,ng)=0

        if(ixb(nm,ng).eq.2) then
           nxbeg(nm,ng)=1
           ixoff(nm,ng)=0
           nxbegc(nm,ng)=2
           ibcflg(nm,ng)=ibcflg(nm,ng)+1
        else
           nxbeg(nm,ng)=ixb(nm,ng)-nbndx
           ixoff(nm,ng)=nxbeg(nm,ng)-1
           nxbegc(nm,ng)=nbndx+1
        endif

        if(ixe(nm,ng).eq.nx-1) then
           nxend(nm,ng)=nx
           ibcflg(nm,ng)=ibcflg(nm,ng)+2
           nxendc(nm,ng)=(1+ixe(nm,ng)-ixb(nm,ng))+nxbegc(nm,ng)-1
        else
           nxend(nm,ng)=ixe(nm,ng)+nbndx
           nxendc(nm,ng)=(ixe(nm,ng)-ixb(nm,ng))+nxbegc(nm,ng)
        endif

        if(iyb(nm,ng).eq.2) then
           nybeg(nm,ng)=1
           iyoff(nm,ng)=0
           nybegc(nm,ng)=2
           ibcflg(nm,ng)=ibcflg(nm,ng)+4
        else
           nybeg(nm,ng)=iyb(nm,ng)-nbndy
           iyoff(nm,ng)=nybeg(nm,ng)-1
           nybegc(nm,ng)=nbndy+1
        endif

        if(iye(nm,ng).eq.ny-1) then
           nyend(nm,ng)=ny
           ibcflg(nm,ng)=ibcflg(nm,ng)+8
           nyendc(nm,ng)=(1+iye(nm,ng)-iyb(nm,ng))+nybegc(nm,ng)-1
        else
           nyend(nm,ng)=iye(nm,ng)+nbndy
           nyendc(nm,ng)=(iye(nm,ng)-iyb(nm,ng))+nybegc(nm,ng)
        endif

     enddo



     !   print*,'Grid:',ng
     !   print *,'   nm nxbeg nxend nybeg nyend'  &
     !        ,' nxbegc nxendc nybegc nyendc ixoff iyoff ibcflg'

     do nm=1,nmachs
        npxy(nm,ng)=(nxend(nm,ng)-nxbeg(nm,ng)+1)  &
             *(nyend(nm,ng)-nybeg(nm,ng)+1)

        !            print '(16i6)',nm
        !     +           ,nxbeg(nm,ng),nxend(nm,ng),nybeg(nm,ng),nyend(nm,ng)
        !     +           ,nxbegc(nm,ng),nxendc(nm,ng),nybegc(nm,ng)
        !     +           ,nyendc(nm,ng)
        !     +           ,ixoff(nm,ng),iyoff(nm,ng),ibcflg(nm,ng)
     enddo

  enddo

  return
end subroutine PAR_decomp_bounds




subroutine dump_Domain_Decomposition()
  use mem_grid, only: &
       ngrids,grid_g

  use rpara, only: &
       nmachs,     &
       ixb,        &
       ixe,        &
       iyb,        &
       iye

  implicit none
  integer :: ngr
  integer :: jnode
  integer :: ncols
  integer :: i,j

  write(*,"(a)") ' === Domain decomposition ==='
  write(*,"(a)") '    grid  node x-beg x-end y-beg y-end  cols'
  do ngr = 1, ngrids
     do jnode = 1,nmachs
        do i=ixb(jnode,ngr),ixe(jnode,ngr)
          do j=iyb(jnode,ngr),iye(jnode,ngr)
            grid_g(ngr)%fmynum(i,j)=real(jnode)
          end do
        end do
        ncols= (1+ixe(jnode,ngr)-ixb(jnode,ngr))  &
             *(1+iye(jnode,ngr)-iyb(jnode,ngr))
        write(*,"('  ',7i6,f12.4)") ngr,jnode,ixb(jnode,ngr),ixe(jnode,ngr)  &
             ,iyb(jnode,ngr),iye(jnode,ngr),ncols
     enddo
     write(*,"(a)")
  end do
end subroutine dump_Domain_Decomposition
  


