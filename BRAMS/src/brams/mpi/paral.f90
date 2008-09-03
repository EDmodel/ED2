!---------------------------------------------------------------------------


subroutine node_sendanl(vtype)
  !     Send just the 'type' variables back to the master process

  use node_mod
  use mem_grid
  use var_tables
  use mem_scratch

  implicit none

  character*(*) vtype
  include 'interface.h'
  include 'mpif.h'

  integer :: ierr,msiii

  integer, parameter :: ntags=5  !4 ! Changed by ALF
  integer :: msgtags(ntags)
  integer :: msgnum,nv,ng
  integer :: mxyp,mxyzp,nmp,ind,isend!,msgnum

  msiii = 0
  msgtags(1)=mynum

  do ng=1,ngrids

     msgtags(2)=ng

     do nv=1, num_var(ng)

        if(vtype == 'LITE' .and. vtab_r(nv,ng)%ilite == 1) then
           msiii = msiii + 1

           msgtags(3)= vtab_r(nv,ng)%npts
           msgtags(4)= vtab_r(nv,ng)%idim_type
           msgtags(5)= nv

           call MPI_Send(msgtags,ntags,MPI_INTEGER,master_num,            &
                         230000+mynum*100+msiii,MPI_COMM_WORLD,ierr)
           call MPI_Send(vtab_r(nv,ng)%var_p, vtab_r(nv,ng)%npts,MPI_REAL, &
                         master_num,240000+mynum*100+msiii,MPI_COMM_WORLD,ierr)

        elseif( (vtype == 'MEAN' .or. vtype == 'BOTH')  &
             .and. vtab_r(nv,ng)%imean == 1) then
           msiii = msiii + 1

           msgtags(3)= vtab_r(nv,ng)%npts
           msgtags(4)= vtab_r(nv,ng)%idim_type
           msgtags(5)= nv

           call MPI_Send(msgtags,ntags,MPI_INTEGER,master_num,            &
                         230000+mynum*100+msiii,MPI_COMM_WORLD,ierr)
           call MPI_Send(vtab_r(nv,ng)%var_m, vtab_r(nv,ng)%npts,MPI_REAL, &
                         master_num,240000+mynum*100+msiii,MPI_COMM_WORLD,ierr)
        endif
     enddo
  enddo

  return
end subroutine node_sendanl

!     ****************************************************************


subroutine master_getanl(vtype)
  !     Get just the 'type' variables from the nodes

  use mem_grid
  use rpara
  use mem_scratch
  use var_tables
  use mem_aerad, only: nwave

  implicit none

  character(len=*) :: vtype

  include 'interface.h'
  include 'mpif.h'

  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: ierr,msiii,nmi


  integer, parameter :: ntags=5  !4 ! Changed by ALF
  integer :: msgtags(ntags)
  integer :: numvars,ng,nummess,nvvv,ibytes,msgtyp,ihostnum
  integer :: nm,il1,ir2,jb1,jt2,mxp,myp,mxyp
  integer :: nv,npts,idim_type

  ! Checking the amount of information (messages) to be received
  numvars = 0
  do ng = 1, ngrids
     do nv = 1, num_var(ng)
        select case (trim(vtype))
        case ('MEAN','BOTH')
           if (vtab_r(nv,ng)%imean == 1) numvars=numvars+1
        case ('LITE')
           if (vtab_r(nv,ng)%ilite == 1) numvars=numvars+1
        end select
     enddo
  enddo

  do nvvv =1,numvars
     msiii=msiii+1
     call MPI_Recv(msgtags,ntags,MPI_INTEGER,machnum(nmi),    &
                   230000+machnum(nmi)*100+msiii,MPI_COMM_WORLD,status,ierr)
     nm = msgtags(1)        ! Number of process for the message received
     ng=msgtags(2)          ! Number of the Grid
     npts=msgtags(3)        ! Total number of points for the real array
     idim_type=msgtags(4)   ! Dimension of the real array
     nv=msgtags(5)          ! Number of the variable in: vtab_r(nv,ng)

     ! Calculating the position (plane coordinates) to fit on grid
     il1=nxbegc(nm,ng)
     ir2=nxendc(nm,ng)
     jb1=nybegc(nm,ng)
     jt2=nyendc(nm,ng)

     if(iand(ibcflg(nm,ng),1).ne.0) il1=il1-1
     if(iand(ibcflg(nm,ng),2).ne.0) ir2=ir2+1
     if(iand(ibcflg(nm,ng),4).ne.0) jb1=jb1-1
     if(iand(ibcflg(nm,ng),8).ne.0) jt2=jt2+1

     mxp = nxend(nm,ng) - nxbeg(nm,ng) + 1
     myp = nyend(nm,ng) - nybeg(nm,ng) + 1
     mxyp = mxp * myp
     
     call MPI_Recv(scratch%scr1(1),npts,MPI_REAL,nm,    &
                   240000+nm*100+msiii,MPI_COMM_WORLD,status,ierr)


     select case (trim(vtype))
     case ('MEAN','BOTH')
        ! Copying the buffer scratch in the real array
        select case (idim_type)
        case (2) 
           call ex_2_buff(vtab_r(nv,ng)%var_m,scratch%scr1(1)  &
                ,nnxp(ng),nnyp(ng)  &
                ,mxp,myp  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

        case (3)
           call ex_3_buff(vtab_r(nv,ng)%var_m,scratch%scr1(1)  &
                ,nnzp(ng),nnxp(ng),nnyp(ng)  &
                ,nnzp(ng),mxp,myp  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

        case (4)
           call ex_4_buff(vtab_r(nv,ng)%var_m,scratch%scr1(1)  &
                ,nzg,nnxp(ng),nnyp(ng),npatch  &
                ,nzg,mxp,myp,npatch  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

        case (5)
           call ex_4_buff(vtab_r(nv,ng)%var_m,scratch%scr1(1)  &
                ,nzs,nnxp(ng),nnyp(ng),npatch  &
                ,nzs,mxp,myp,npatch  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

        case (6)
           call ex_2p_buff(vtab_r(nv,ng)%var_m,scratch%scr1(1)  &
                ,nnxp(ng),nnyp(ng),npatch  &
                ,mxp,myp,npatch  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

        case (7)
           call ex_buff_carma(vtab_r(nv,ng)%var_m,scratch%scr1(1)  &
                ,nnxp(ng),nnyp(ng),nwave,mxp,myp,nwave  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

        case (8)
           !tridimensional variables, with vertical dimension 
           !being the number soil levels
           call ex_2p_buff(vtab_r(nv,ng)%var_m,scratch%scr1(1)  &
                ,nnxp(ng),nnyp(ng),nzg,mxp,myp,nzg  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
        end select

     case ('LITE')
        ! Copying the buffer scratch in the real array
        select case (idim_type)
        case (2)
           call ex_2_buff(vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
                ,nnxp(ng),nnyp(ng)  &
                ,mxp,myp  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

        case (3)
           call ex_3_buff(vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
                ,nnzp(ng),nnxp(ng),nnyp(ng)  &
                ,nnzp(ng),mxp,myp  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

        case (4)
           call ex_4_buff(vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
                ,nzg,nnxp(ng),nnyp(ng),npatch  &
                ,nzg,mxp,myp,npatch  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

        case (5)
           call ex_4_buff(vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
                ,nzs,nnxp(ng),nnyp(ng),npatch  &
                ,nzs,mxp,myp,npatch  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

        case (6)
           call ex_2p_buff(vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
                ,nnxp(ng),nnyp(ng),npatch  &
                ,mxp,myp,npatch  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

        case (7)
           call ex_buff_carma(vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
                ,nnxp(ng),nnyp(ng),nwave,mxp,myp,nwave  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

        case (8)
           ![ED2-MLO: tridimensional variables, with vertical dimension 
           !being the number soil levels
           call ex_2p_buff(vtab_r(nv,ng)%var_p,scratch%scr1(1)  &
                ,nnxp(ng),nnyp(ng),nzg,mxp,myp,nzg  &
                ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)
        end select
     end select
  
  end do

  return
end subroutine master_getanl

