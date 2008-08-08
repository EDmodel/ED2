!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine master_sendinit()

  use var_tables
  use mem_scratch
  use mem_cuparm
  use mem_varinit
  use mem_grid
  use io_params
  use rpara
  ! For CATT
  use mem_aerad, only: nwave

  implicit none

  integer :: ierr
  include 'mpif.h'
  include 'interface.h'
  integer :: nm,nmp,ind,nwords,mxp,myp,mxyp,mxyzp,nv,npts,ng,i

  do nm = 1,nmachs
     call MPI_Send(vtime1,1,MPI_DOUBLE_PRECISION,machnum(nm),110,MPI_COMM_WORLD,ierr)
     call MPI_Send(vtime2,1,MPI_DOUBLE_PRECISION,machnum(nm),111,MPI_COMM_WORLD,ierr)
     call MPI_Send(htime1,1,MPI_DOUBLE_PRECISION,machnum(nm),112,MPI_COMM_WORLD,ierr)
     call MPI_Send(htime2,1,MPI_DOUBLE_PRECISION,machnum(nm),113,MPI_COMM_WORLD,ierr)
     call MPI_Send(condtime1,1,MPI_DOUBLE_PRECISION,machnum(nm),114,MPI_COMM_WORLD,ierr)
     call MPI_Send(condtime2,1,MPI_DOUBLE_PRECISION,machnum(nm),115,MPI_COMM_WORLD,ierr)
     call MPI_Send(cutime1,1,MPI_DOUBLE_PRECISION,machnum(nm),116,MPI_COMM_WORLD,ierr)
     call MPI_Send(cutime2,1,MPI_DOUBLE_PRECISION,machnum(nm),117,MPI_COMM_WORLD,ierr)
     call MPI_Send(ssttime1,maxgrds,MPI_DOUBLE_PRECISION,machnum(nm),118,MPI_COMM_WORLD,ierr)
     call MPI_Send(ssttime2,maxgrds,MPI_DOUBLE_PRECISION,machnum(nm),119,MPI_COMM_WORLD,ierr)
     call MPI_Send(ndvitime1,maxgrds,MPI_DOUBLE_PRECISION,machnum(nm),120,MPI_COMM_WORLD,ierr)
     call MPI_Send(ndvitime2,maxgrds,MPI_DOUBLE_PRECISION,machnum(nm),121,MPI_COMM_WORLD,ierr)

     do ng = 1,ngrids

        mxp = nxend(nm,ng) - nxbeg(nm,ng) + 1
        myp = nyend(nm,ng) - nybeg(nm,ng) + 1
        mxyp = mxp * myp
        mxyzp = mxyp * nnzp(ng)

        !  Send variables that the nodes will only have the subdomain portion of.

        do nv = 1,num_var(ng)

           if ( vtab_r(nv,ng)%impti == 1) then

              if ( vtab_r(nv,ng)%idim_type == 2) then

                 call mk_2_buff(vtab_r(nv,ng)%var_p  &
                      ,scratch%scr2(1),nnxp(ng),nnyp(ng),mxp,myp  &
                      ,nxbeg(nm,ng),nxend(nm,ng)  &
                      ,nybeg(nm,ng),nyend(nm,ng))
                 npts=mxyp

              elseif ( vtab_r(nv,ng)%idim_type == 3) then

                 call mk_3_buff(vtab_r(nv,ng)%var_p  &
                      ,scratch%scr2(1),nnzp(ng),nnxp(ng),nnyp(ng)  &
                      ,nnzp(ng),mxp,myp  &
                      ,nxbeg(nm,ng),nxend(nm,ng)  &
                      ,nybeg(nm,ng),nyend(nm,ng))
                 npts=mxyzp

              elseif ( vtab_r(nv,ng)%idim_type == 4) then

                 call mk_4_buff(vtab_r(nv,ng)%var_p  &
                      ,scratch%scr2(1),nzg,nnxp(ng),nnyp(ng)  &
                      ,npatch,nzg,mxp,myp,npatch &
                      ,nxbeg(nm,ng),nxend(nm,ng)  &
                      ,nybeg(nm,ng),nyend(nm,ng))
                 npts=mxyp*nzg*npatch

              elseif ( vtab_r(nv,ng)%idim_type == 5) then

                 call mk_4_buff(vtab_r(nv,ng)%var_p  &
                      ,scratch%scr2(1),nzs,nnxp(ng),nnyp(ng)  &
                      ,npatch,nzs,mxp,myp,npatch &
                      ,nxbeg(nm,ng),nxend(nm,ng)  &
                      ,nybeg(nm,ng),nyend(nm,ng))
                 npts=mxyp*nzs*npatch

              elseif ( vtab_r(nv,ng)%idim_type == 6) then

                 call mk_2p_buff(vtab_r(nv,ng)%var_p  &
                      ,scratch%scr2(1),nnxp(ng),nnyp(ng)  &
                      ,npatch,mxp,myp,npatch &
                      ,nxbeg(nm,ng),nxend(nm,ng)  &
                      ,nybeg(nm,ng),nyend(nm,ng))
                 npts=mxyp*npatch
              elseif ( vtab_r(nv,ng)%idim_type == 7) then

                 call mk_buff_carma(vtab_r(nv,ng)%var_p  &
                      ,scratch%scr2(1) &
                      ,nnxp(ng),nnyp(ng),nwave &
                      ,mxp,myp,nwave &
                      ,nxbeg(nm,ng),nxend(nm,ng)  &
                      ,nybeg(nm,ng),nyend(nm,ng))
                 npts=mxyp*nwave

              elseif ( vtab_r(nv,ng)%idim_type == 8) then

                 call mk_4_buff(vtab_r(nv,ng)%var_p  &
                      ,scratch%scr2(1),nnzp(ng),nnxp(ng),nnyp(ng)  &
                      ,nclouds,nnzp(ng),mxp,myp,nclouds &
                      ,nxbeg(nm,ng),nxend(nm,ng)  &
                      ,nybeg(nm,ng),nyend(nm,ng))
                 npts=mxyzp*nclouds

              elseif ( vtab_r(nv,ng)%idim_type == 9 ) then

                 call mk_2p_buff(vtab_r(nv,ng)%var_p  &
                      ,scratch%scr2(1),nnxp(ng),nnyp(ng)  &
                      ,nclouds,mxp,myp,nclouds &
                      ,nxbeg(nm,ng),nxend(nm,ng)  &
                      ,nybeg(nm,ng),nyend(nm,ng))
                 npts=mxyp*nclouds

              endif

              call MPI_Send(scratch%scr2(1),npts,MPI_REAL,machnum(nm),  &
                            110000+ng*1000+nv,MPI_COMM_WORLD,ierr)

           endif

        enddo

     enddo
  enddo

  return
end subroutine master_sendinit

!---------------------------------------------------------------------------

!For Carma - CATT
subroutine mk_buff_carma(a,b,n1,n2,n3,m1,m2,m3,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,m1,m2,m3,i1,i2,j1,j2
  real :: a(n1,n2,n3),b(m1,m2,m3)

  b(1:m1,1:m2,1:m3)=a(i1:i2,j1:j2,1:n3)

  return
end subroutine mk_buff_carma

!---------------------------------------------------------------------------

subroutine mk_2_buff(a,b,n1,n2,m1,m2,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,m1,m2,i1,i2,j1,j2
  real :: a(n1,n2),b(m1,m2)

  b(1:m1,1:m2)=a(i1:i2,j1:j2)

  return
end subroutine mk_2_buff

!---------------------------------------------------------------------------

subroutine mk_2p_buff(a,b,n1,n2,n3,m1,m2,m3,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,m1,m2,m3,i1,i2,j1,j2
  real :: a(n1,n2,n3),b(m1,m2,m3)

  b(1:m1,1:m2,1:m3)=a(i1:i2,j1:j2,1:n3)

  return
end subroutine mk_2p_buff

!---------------------------------------------------------------------------

subroutine mk_3_buff(a,b,n1,n2,n3,m1,m2,m3,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,m1,m2,m3,i1,i2,j1,j2
  real :: a(n1,n2,n3),b(m1,m2,m3)

  b(1:m1,1:m2,1:m3)=a(1:n1,i1:i2,j1:j2)

  return
end subroutine mk_3_buff

!---------------------------------------------------------------------------

subroutine mk_4_buff(a,b,n1,n2,n3,n4,m1,m2,m3,m4,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,n4,m1,m2,m3,m4,i1,i2,j1,j2
  real :: a(n1,n2,n3,m4),b(m1,m2,m3,m4)

  b(1:m1,1:m2,1:m3,1:m4)=a(1:n1,i1:i2,j1:j2,1:n4)

  return
end subroutine mk_4_buff

!     ****************************************************************

subroutine node_getinit()

  use var_tables
  use mem_scratch
  use node_mod
  use mem_cuparm
  use mem_varinit
  use mem_grid
  use io_params

  implicit none

  include 'mpif.h'
  include 'interface.h'
  integer :: ierr
  integer :: mxyp,mxyzp,nv,nmp
  integer :: ng


  call MPI_Recv(vtime1,1,MPI_DOUBLE_PRECISION,master_num,110,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call MPI_Recv(vtime2,1,MPI_DOUBLE_PRECISION,master_num,111,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call MPI_Recv(htime1,1,MPI_DOUBLE_PRECISION,master_num,112,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call MPI_Recv(htime2,1,MPI_DOUBLE_PRECISION,master_num,113,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call MPI_Recv(condtime1,1,MPI_DOUBLE_PRECISION,master_num,114,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call MPI_Recv(condtime2,1,MPI_DOUBLE_PRECISION,master_num,115,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call MPI_Recv(cutime1,1,MPI_DOUBLE_PRECISION,master_num,116,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call MPI_Recv(cutime2,1,MPI_DOUBLE_PRECISION,master_num,117,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call MPI_Recv(ssttime1,maxgrds,MPI_DOUBLE_PRECISION,master_num,118,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call MPI_Recv(ssttime2,maxgrds,MPI_DOUBLE_PRECISION,master_num,119,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call MPI_Recv(ndvitime1,maxgrds,MPI_DOUBLE_PRECISION,master_num,120,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call MPI_Recv(ndvitime2,maxgrds,MPI_DOUBLE_PRECISION,master_num,121,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)


  do ng = 1,ngrids

     do nv=1,num_var(ng)

        if (vtab_r(nv,ng)%impti == 1) then
           call MPI_Recv(vtab_r(nv,ng)%var_p,vtab_r(nv,ng)%npts,MPI_REAL, &
                         master_num,110000+ng*1000+nv,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
        endif

     enddo

  enddo

  return
end subroutine node_getinit

!     *****************************************************************

subroutine node_sendall()

  use node_mod
  use mem_grid
  use var_tables
  use mem_scratch
  implicit none

  include 'interface.h'
  include 'mpif.h'

  integer :: ierr,msiii
  integer, parameter :: ntags=5
  integer :: msgtags(ntags),mxyp,mxyzp,nmp,nv,ind,ng

  msiii = 0
  do ng=1,ngrids

     msgtags(1)=mynum
     msgtags(2)=ng

     do nv=1,num_var(ng)

        if (vtab_r(nv,ng)%impt3 == 1) then
           msgtags(3)= vtab_r(nv,ng)%npts
           msgtags(4)= vtab_r(nv,ng)%idim_type
           msgtags(5)= nv

           msiii = msiii +1
           call MPI_Send(msgtags,ntags,MPI_INTEGER,master_num,            &
                         210000+mynum*100+msiii,MPI_COMM_WORLD,ierr)
           call MPI_Send(vtab_r(nv,ng)%var_p,vtab_r(nv,ng)%npts,MPI_REAL, &
                         master_num,220000+mynum*100+msiii,MPI_COMM_WORLD,ierr)
        endif
        
     enddo

  enddo


  return
end subroutine node_sendall

!     ****************************************************************

subroutine master_getall()

  use mem_grid
  use rpara
  use var_tables
  use mem_scratch
  use mem_aerad, only: nwave
  use mem_cuparm, only: nclouds
  
  implicit none

  include 'interface.h'
  include 'mpif.h'

  integer, parameter :: ntags=5
  integer :: msgtags(ntags),numvars,nvvv,mxp,myp,mxyp,mxyzp
  integer :: nm,il1,ir2,jb1,jt2,ind,nv,idim_type,npts,ng
  integer :: ierr,msiii,nmi

  numvars=0
  do ng=1,ngrids
     do nv=1,num_var(ng)
        if (vtab_r(nv,ng)%impt3 == 1) numvars=numvars+1
     enddo
  enddo

  do nmi = 1,nmachs
     msiii = 0
     do nvvv=1,numvars
        msiii = msiii + 1
        call MPI_Recv(msgtags,ntags,MPI_INTEGER,machnum(nmi),   &
                      210000+machnum(nmi)*100+msiii,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
        nm=msgtags(1)
        ng=msgtags(2)
        npts=msgtags(3)
        idim_type=msgtags(4)
        nv=msgtags(5)

        il1=nxbegc(nm,ng)
        ir2=nxendc(nm,ng)
        jb1=nybegc(nm,ng)
        jt2=nyendc(nm,ng)
        call MPI_Recv(scratch%scr1(1),npts,MPI_REAL,nm,         &
                      220000+nm*100+msiii,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
        if(iand(ibcflg(nm,ng),1).ne.0) il1=il1-1
        if(iand(ibcflg(nm,ng),2).ne.0) ir2=ir2+1
        if(iand(ibcflg(nm,ng),4).ne.0) jb1=jb1-1
        if(iand(ibcflg(nm,ng),8).ne.0) jt2=jt2+1

        mxp = nxend(nm,ng) - nxbeg(nm,ng) + 1
        myp = nyend(nm,ng) - nybeg(nm,ng) + 1
        mxyp = mxp * myp

        if (idim_type == 2) then
           call ex_2_buff(vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
                ,nnxp(ng),nnyp(ng) ,mxp,myp  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
        elseif (idim_type == 3) then
           call ex_3_buff(vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
                ,nnzp(ng),nnxp(ng),nnyp(ng),nnzp(ng),mxp,myp  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
        elseif (idim_type == 4) then
           call ex_4_buff(vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
                ,nzg,nnxp(ng),nnyp(ng),npatch,nzg,mxp,myp,npatch  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
        elseif (idim_type == 5) then
           call ex_4_buff(vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
                ,nzs,nnxp(ng),nnyp(ng),npatch,nzs,mxp,myp,npatch  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
        elseif (idim_type == 6) then
           call ex_2p_buff(vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
                ,nnxp(ng),nnyp(ng),npatch,mxp,myp,npatch  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
       elseif (idim_type == 7) then
          call ex_buff_carma(vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
               ,nnxp(ng),nnyp(ng),nwave,mxp,myp,nwave  &
               ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
       elseif (idim_type == 8) then
          call ex_4_buff(vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
               ,nnzp(ng),nnxp(ng),nnyp(ng),nclouds,nnzp(ng),mxp,myp,nclouds  &
               ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
        elseif (idim_type == 9) then
           call ex_2p_buff(vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
                ,nnxp(ng),nnyp(ng),nclouds,mxp,myp,nclouds      &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
       endif

     enddo
  enddo
  return
end subroutine master_getall

!---------------------------------------------------------------------------

subroutine ex_buff_carma(a,b,n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2)
  !i0 = ixoff(nm,ng)
  !j0 = iyoff(nm,ng)
  !i1 = il1
  !i2 = ir2
  !j1 = jb1
  !j2 = jt2
  implicit none
  integer :: n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2
  real :: a(n1,n2,n3),b(m1,m2,m3)

  a(i1+i0:i2+i0,j1+j0:j2+j0,1:n3) = b(i1:i2,j1:j2,1:m3)

  return
end subroutine ex_buff_carma

!---------------------------------------------------------------------------

subroutine ex_2_buff(a,b,n1,n2,m1,m2,i0,j0,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,m1,m2,i0,j0,i1,i2,j1,j2
  real :: a(n1,n2),b(m1,m2)

  a(i1+i0:i2+i0,j1+j0:j2+j0) = b(i1:i2,j1:j2)

  return
end subroutine ex_2_buff

!---------------------------------------------------------------------------

subroutine ex_3_buff(a,b,n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2
  real :: a(n1,n2,n3),b(m1,m2,m3)

  a(1:n1,i1+i0:i2+i0,j1+j0:j2+j0) = b(1:m1,i1:i2,j1:j2)

  return
end subroutine ex_3_buff

!---------------------------------------------------------------------------

subroutine ex_4_buff(a,b,n1,n2,n3,n4,m1,m2,m3,m4,i0,j0,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,n4,m1,m2,m3,m4,i0,j0,i1,i2,j1,j2
  real :: a(n1,n2,n3,n4),b(m1,m2,m3,m4)

  a(1:n1,i1+i0:i2+i0,j1+j0:j2+j0,1:n4) = b(1:m1,i1:i2,j1:j2,1:m4)

  return
end subroutine ex_4_buff

!---------------------------------------------------------------------------

subroutine ex_2p_buff(a,b,n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2)
  implicit none
  integer :: n1,n2,n3,m1,m2,m3,i0,j0,i1,i2,j1,j2
  real :: a(n1,n2,n3),b(m1,m2,m3)

  a(i1+i0:i2+i0,j1+j0:j2+j0,1:n3) = b(i1:i2,j1:j2,1:m3)

  return
end subroutine ex_2p_buff
