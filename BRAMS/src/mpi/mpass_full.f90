!===================================== Change Log =========================================!
! 5.0.0                                                                                    !
!                                                                                          !
! MLO - 09/29/08 Adding the cloud spectral size related dimensions (8 and 9). Eliminated   !
!       the carma subroutine, since they did the exact same thing as the "2p" ones. Also   !
!       eliminated the mz and mpatch (and related) arguments from the packing and unpack-  !
!       ing subroutines since they were always the same as the full domain.                !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine is called by the head node, and it will wrap all variables into       !
! smaller subdomains and send them to the slave nodes.                                     !
!------------------------------------------------------------------------------------------!
subroutine master_sendinit()

   use var_tables
   use mem_scratch
   use mem_cuparm
   use mem_varinit
   use mem_grid
   use io_params
   use rpara
   use mem_aerad, only: nwave

   implicit none

   !----- External variable declaration ---------------------------------------------------!
   include 'mpif.h'
   include 'interface.h'
   !----- Local variables -----------------------------------------------------------------!
   integer :: ierr
   integer :: nm,nmp,ind,nwords,mxp,myp,mxyp,mxyzp,nv,npts,ng,i
   integer :: mpiid
   !---------------------------------------------------------------------------------------!

   machloop: do nm = 1,nmachs
      !----- Sending time information -----------------------------------------------------!
      call MPI_Send(vtime1,1,MPI_DOUBLE_PRECISION,machnum(nm),200,MPI_COMM_WORLD,ierr)
      call MPI_Send(vtime2,1,MPI_DOUBLE_PRECISION,machnum(nm),201,MPI_COMM_WORLD,ierr)
      call MPI_Send(htime1,1,MPI_DOUBLE_PRECISION,machnum(nm),202,MPI_COMM_WORLD,ierr)
      call MPI_Send(htime2,1,MPI_DOUBLE_PRECISION,machnum(nm),203,MPI_COMM_WORLD,ierr)
      call MPI_Send(condtime1,1,MPI_DOUBLE_PRECISION,machnum(nm),204,MPI_COMM_WORLD,ierr)
      call MPI_Send(condtime2,1,MPI_DOUBLE_PRECISION,machnum(nm),205,MPI_COMM_WORLD,ierr)
      call MPI_Send(cutime1,1,MPI_DOUBLE_PRECISION,machnum(nm),206,MPI_COMM_WORLD,ierr)
      call MPI_Send(cutime2,1,MPI_DOUBLE_PRECISION,machnum(nm),207,MPI_COMM_WORLD,ierr)
      call MPI_Send(ssttime1,maxgrds,MPI_DOUBLE_PRECISION,machnum(nm),208,MPI_COMM_WORLD   &
                   ,ierr)
      call MPI_Send(ssttime2,maxgrds,MPI_DOUBLE_PRECISION,machnum(nm),209,MPI_COMM_WORLD   &
                   ,ierr)
      call MPI_Send(ndvitime1,maxgrds,MPI_DOUBLE_PRECISION,machnum(nm),210,MPI_COMM_WORLD  &
                   ,ierr)
      call MPI_Send(ndvitime2,maxgrds,MPI_DOUBLE_PRECISION,machnum(nm),211,MPI_COMM_WORLD  &
                   ,ierr)

      mpiid = 100000
      gridloop: do ng = 1,ngrids

         !----- Defining grid dimensions --------------------------------------------------!
         mxp   = nxend(nm,ng) - nxbeg(nm,ng) + 1
         myp   = nyend(nm,ng) - nybeg(nm,ng) + 1
         mxyp  = mxp * myp
         mxyzp = mxyp * nnzp(ng)

         !----- Send variables that the nodes will only have the subdomain portion of. ----!
         varloop: do nv = 1,num_var(ng)
            if ( vtab_r(nv,ng)%impti == 1) then
               !----- This will create a unique ID for this package -----------------------!
               mpiid = mpiid + 1

               !---------------------------------------------------------------------------!
               !     We will always send scratch%scr2, which is a scratch variable enough  !
               ! large to receive any subdomain. Here all we do is to copy the right       !
               ! amount of data to scr2 depending on the variable type, and then define    !
               ! npts accordingly.                                                         !
               !---------------------------------------------------------------------------!
               select case(vtab_r(nv,ng)%idim_type)
               case (2) !----- 2D variables (nxp,nyp) -------------------------------------!
                  call mk_2_buff(vtab_r(nv,ng)%var_p ,scratch%scr2,nnxp(ng),nnyp(ng)       &
                                ,mxp,myp,nxbeg(nm,ng),nxend(nm,ng),nybeg(nm,ng)            &
                                ,nyend(nm,ng))
                  npts=mxyp

               case (3) !----- 3D variables (nzp,nxp,nyp) ---------------------------------!
                  call mk_3_buff(vtab_r(nv,ng)%var_p,scratch%scr2,nnzp(ng),nnxp(ng)        &
                                ,nnyp(ng),mxp,myp,nxbeg(nm,ng),nxend(nm,ng),nybeg(nm,ng)   &
                                ,nyend(nm,ng))
                  npts=mxyzp

               case (4) !----- 4D variables (nzg,nxp,nyp,npatch) --------------------------!
                  call mk_4_buff(vtab_r(nv,ng)%var_p,scratch%scr2,nzg,nnxp(ng),nnyp(ng)    &
                                ,npatch,mxp,myp,nxbeg(nm,ng),nxend(nm,ng),nybeg(nm,ng)     &
                                ,nyend(nm,ng))
                  npts=mxyp*nzg*npatch

               case (5) !----- 4D variables (nzs,nxp,nyp,npatch) --------------------------!
                  call mk_4_buff(vtab_r(nv,ng)%var_p,scratch%scr2,nzs,nnxp(ng),nnyp(ng)    &
                                ,npatch,mxp,myp,nxbeg(nm,ng),nxend(nm,ng),nybeg(nm,ng)     &
                                ,nyend(nm,ng))
                  npts=mxyp*nzs*npatch

               case (6) !----- 3D variables (nxp,nyp,npatch) ------------------------------!
                  call mk_2p_buff(vtab_r(nv,ng)%var_p,scratch%scr2,nnxp(ng),nnyp(ng)       &
                                 ,npatch,mxp,myp,nxbeg(nm,ng),nxend(nm,ng),nybeg(nm,ng)    &
                                 ,nyend(nm,ng))
                  npts=mxyp*npatch

               case (7) !----- 3D variables (nxp,nyp,nwave) -------------------------------!
                  call mk_2p_buff(vtab_r(nv,ng)%var_p,scratch%scr2,nnxp(ng),nnyp(ng)       &
                                    ,nwave,mxp,myp,nxbeg(nm,ng),nxend(nm,ng),nybeg(nm,ng)  &
                                    ,nyend(nm,ng))
                  npts=mxyp*nwave

               case (8) !----- 4D variables (nzp,nxp,nyp,nclouds) -------------------------!
                  call mk_4_buff(vtab_r(nv,ng)%var_p,scratch%scr2,nnzp(ng),nnxp(ng)        &
                                ,nnyp(ng),nclouds,mxp,myp,nxbeg(nm,ng),nxend(nm,ng)        &
                                ,nybeg(nm,ng),nyend(nm,ng))
                  npts=mxyzp*nclouds

               case (9) !----- 3D variables (nxp,nyp,nclouds) -----------------------------!
                  call mk_2p_buff(vtab_r(nv,ng)%var_p,scratch%scr2,nnxp(ng),nnyp(ng)       &
                                 ,nclouds,mxp,myp,nxbeg(nm,ng),nxend(nm,ng),nybeg(nm,ng)   &
                                 ,nyend(nm,ng))
                  npts=mxyp*nclouds

               end select
               !----- Sending the buffer to the node --------------------------------------!
               call MPI_Send(scratch%scr2,npts,MPI_REAL,machnum(nm),mpiid                  &
                            ,MPI_COMM_WORLD,ierr)
            end if
         end do varloop
      end do gridloop
   end do machloop

   return
end subroutine master_sendinit
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine simply copies the 2D subdomain data to the buffer matrix. Here we    !
! use the following notation: n? is the full domain dimension, whereas m? is the machine   !
! subdomain dimension.                                                                     !
!------------------------------------------------------------------------------------------!
subroutine mk_2_buff(mydata,buff,nx,ny,mx,my,ibeg,iend,jbeg,jend)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer               , intent(in)  :: nx,ny,mx,my,ibeg,iend,jbeg,jend
   real, dimension(nx,ny), intent(in)  :: mydata
   real, dimension(mx,my), intent(out) :: buff
   !----- Local variables. ----------------------------------------------------------------!
   integer                             :: mi
   integer                             :: mj
   integer                             :: ni
   integer                             :: nj
   !---------------------------------------------------------------------------------------!
   nj = jbeg
   do mj=1,my
      ni=ibeg

      do mi=1,mx
         buff(mi,mj)=mydata(ni,nj)
         ni = ni + 1
      end do

      nj = nj +1
   end do

   return
end subroutine mk_2_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine simply copies a 2D subdomain data that has an extra dimension (such  !
! as patch, wavelength, or cloud size) to the buffer matrix. X and Y are the only dimen-   !
! sions allowed to have sub-domains. Here we use the following notation: n? is the full    !
! domain dimension, whereas m? is the machine subdomain dimension.                         !
!------------------------------------------------------------------------------------------!
subroutine mk_2p_buff(mydata,buff,nx,ny,ne,mx,my,ibeg,iend,jbeg,jend)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                  , intent(in)  :: nx,ny,ne,mx,my,ibeg,iend,jbeg,jend
   real, dimension(nx,ny,ne), intent(in)  :: mydata
   real, dimension(mx,my,ne), intent(out) :: buff
   !----- Local variables. ----------------------------------------------------------------!
   integer                             :: mi
   integer                             :: mj
   integer                             :: ni
   integer                             :: nj
   integer                             :: nl
   !---------------------------------------------------------------------------------------!

   do nl=1,ne
      nj = jbeg
      do mj=1,my
         ni=ibeg

         do mi=1,mx

            buff(mi,mj,nl)=mydata(ni,nj,nl)

            ni = ni + 1
         end do

         nj = nj +1
      end do
   end do

   return
end subroutine mk_2p_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine simply copies a 3D subdomain data to the buffer matrix. Here we use  !
! the following notation: n? is the full domain dimension, whereas m? is the machine sub-  !
! domain dimension. Note that X and Y are the only dimensions allowed to have sub-domains. !
!------------------------------------------------------------------------------------------!
subroutine mk_3_buff(mydata,buff,nz,nx,ny,mx,my,ibeg,iend,jbeg,jend)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                  , intent(in)  :: nz,nx,ny,mx,my,ibeg,iend,jbeg,jend
   real, dimension(nz,nx,ny), intent(in)  :: mydata
   real, dimension(nz,mx,my), intent(out) :: buff
   !----- Local variables. ----------------------------------------------------------------!
   integer                             :: mi
   integer                             :: mj
   integer                             :: ni
   integer                             :: nj
   integer                             :: nk
   !---------------------------------------------------------------------------------------!

   nj = jbeg
   do mj=1,my
      ni=ibeg

      do mi=1,mx

         do nk=1,nz
            buff(nk,mi,mj)=mydata(nk,ni,nj)
         end do

         ni = ni + 1
      end do

      nj = nj +1
   end do

   return
end subroutine mk_3_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine simply copies a 3D subdomain data that has an extra dimension (such  !
! as patch, wavelength, or cloud size) to the buffer matrix. X and Y are the only dimen-   !
! sions allowed to have sub-domains.                                                       !
!------------------------------------------------------------------------------------------!
subroutine mk_4_buff(mydata,buff,nz,nx,ny,ne,mx,my,ibeg,iend,jbeg,jend)
  implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                     , intent(in)  :: nz,nx,ny,ne,mx,my,ibeg,iend,jbeg,jend
   real, dimension(nz,nx,ny,ne), intent(in)  :: mydata
   real, dimension(nz,mx,my,ne), intent(out) :: buff
   !----- Local variables. ----------------------------------------------------------------!
   integer                             :: mi
   integer                             :: mj
   integer                             :: ni
   integer                             :: nj
   integer                             :: nk
   integer                             :: nl
   !---------------------------------------------------------------------------------------!

   do nl=1,ne
      nj = jbeg
      do mj=1,my
         ni=ibeg

         do mi=1,mx
            do nk=1,nz
               buff(nk,mi,mj,nl)=mydata(nk,ni,nj,nl)
            end do
            ni = ni + 1
         end do

         nj = nj +1
      end do
   end do


   return
end subroutine mk_4_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine is called by the slave nodes to retrieve the full sub-domain inform-  !
! ation sent by the head node.                                                             !
!------------------------------------------------------------------------------------------!
subroutine node_getinit()

   use var_tables
   use mem_scratch
   use node_mod
   use mem_cuparm
   use mem_varinit
   use mem_grid
   use io_params

   implicit none

   !----- External variable declaration ---------------------------------------------------!
   include 'mpif.h'
   include 'interface.h'
   !----- Local variables -----------------------------------------------------------------!
   integer :: ierr
   integer :: mxyp,mxyzp,nv,nmp
   integer :: ng, mpiid
   !---------------------------------------------------------------------------------------!

   !----- Getting the time-related information --------------------------------------------!
   call MPI_Recv(vtime1,1,MPI_DOUBLE_PRECISION,master_num,200,MPI_COMM_WORLD               &
                ,MPI_STATUS_IGNORE,ierr)
   call MPI_Recv(vtime2,1,MPI_DOUBLE_PRECISION,master_num,201,MPI_COMM_WORLD               &
                ,MPI_STATUS_IGNORE,ierr)
   call MPI_Recv(htime1,1,MPI_DOUBLE_PRECISION,master_num,202,MPI_COMM_WORLD               &
                ,MPI_STATUS_IGNORE,ierr)
   call MPI_Recv(htime2,1,MPI_DOUBLE_PRECISION,master_num,203,MPI_COMM_WORLD               &
                ,MPI_STATUS_IGNORE,ierr)
   call MPI_Recv(condtime1,1,MPI_DOUBLE_PRECISION,master_num,204,MPI_COMM_WORLD            &
                ,MPI_STATUS_IGNORE,ierr)
   call MPI_Recv(condtime2,1,MPI_DOUBLE_PRECISION,master_num,205,MPI_COMM_WORLD            &
                ,MPI_STATUS_IGNORE,ierr)
   call MPI_Recv(cutime1,1,MPI_DOUBLE_PRECISION,master_num,206,MPI_COMM_WORLD              &
                ,MPI_STATUS_IGNORE,ierr)
   call MPI_Recv(cutime2,1,MPI_DOUBLE_PRECISION,master_num,207,MPI_COMM_WORLD              &
                ,MPI_STATUS_IGNORE,ierr)
   call MPI_Recv(ssttime1,maxgrds,MPI_DOUBLE_PRECISION,master_num,208,MPI_COMM_WORLD       &
                ,MPI_STATUS_IGNORE,ierr)
   call MPI_Recv(ssttime2,maxgrds,MPI_DOUBLE_PRECISION,master_num,209,MPI_COMM_WORLD       &
                ,MPI_STATUS_IGNORE,ierr)
   call MPI_Recv(ndvitime1,maxgrds,MPI_DOUBLE_PRECISION,master_num,210,MPI_COMM_WORLD      &
                ,MPI_STATUS_IGNORE,ierr)
   call MPI_Recv(ndvitime2,maxgrds,MPI_DOUBLE_PRECISION,master_num,211,MPI_COMM_WORLD      &
                ,MPI_STATUS_IGNORE,ierr)
   !---------------------------------------------------------------------------------------!

   mpiid = 100000
   gridloop: do ng = 1,ngrids
      varloop: do nv=1,num_var(ng)

         if (vtab_r(nv,ng)%impti == 1) then
            mpiid = mpiid + 1
            !------------------------------------------------------------------------------!
            !    Receiving the variable sub-domain. Using the variable table pointer to    !
            ! receive the information. It is already set up to the actual place so this    !
            ! will direct it to the right place.                                           !
            !------------------------------------------------------------------------------!
            call MPI_Recv(vtab_r(nv,ng)%var_p,vtab_r(nv,ng)%npts,MPI_REAL,master_num       &
                         ,mpiid,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
         end if

      end do varloop
   end do gridloop

   return
end subroutine node_getinit
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will be called by the nodes when it is time to send all information   !
! in their subdomains that is expected by the head node.                                   !
!------------------------------------------------------------------------------------------!
subroutine node_sendall()

   use node_mod
   use mem_grid
   use var_tables
   use mem_scratch
   use grid_dims

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'mpif.h'
   include 'interface.h'
   !----- Local constants -----------------------------------------------------------------!
   integer                  , parameter :: ntags=5
   !----- Local variables -----------------------------------------------------------------!
   integer                              :: ierr,mpiidtags,mpiidvtab
   integer                              :: mxyp,mxyzp,nmp,nv,ind,ng
   integer, dimension(ntags)            :: msgtags
   !---------------------------------------------------------------------------------------!

   mpiidtags=200000000 + maxvars*maxgrds*(mynum-1)
   mpiidvtab=300000000 + maxvars*maxgrds*(mynum-1)
   gridloop: do ng=1,ngrids
      !----- Preparing the tags with information pertinent to this grid and this machine. -!
      msgtags(1)=mynum
      msgtags(2)=ng

      varloop: do nv=1,num_var(ng)

         if (vtab_r(nv,ng)%impt3 == 1) then
            !----- Filling the remainder tags with variable specific information. ---------!
            msgtags(3)= vtab_r(nv,ng)%npts
            msgtags(4)= vtab_r(nv,ng)%idim_type
            msgtags(5)= nv
            mpiidtags = mpiidtags + 1
            mpiidvtab = mpiidvtab + 1

            !----- Sending the tags then the variable to the head node --------------------!
            call MPI_Send(msgtags,ntags,MPI_INTEGER,master_num,mpiidtags                   &
                         ,MPI_COMM_WORLD,ierr)
            call MPI_Send(vtab_r(nv,ng)%var_p,vtab_r(nv,ng)%npts,MPI_REAL,master_num       &
                         ,mpiidvtab,MPI_COMM_WORLD,ierr)
         end if
      end do varloop
   end do gridloop


   return
end subroutine node_sendall
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will be called by the head node when it is time to collect the        !
! sub-domain information sent by each slave node.                                          !
!------------------------------------------------------------------------------------------!
subroutine master_getall()
   use mem_grid
   use rpara
   use var_tables
   use mem_scratch
   use mem_cuparm, only: nclouds
   use mem_aerad, only: nwave
   use grid_dims
  
   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'mpif.h'
   include 'interface.h'
   !----- Local constants -----------------------------------------------------------------!
   integer                  , parameter :: ntags=5
   !----- Local variables -----------------------------------------------------------------!
   integer                              :: ierr,nmi,numvars,nvvv,mxp,myp,mxyp,mxyzp
   integer                              :: nm,il1,ir2,jb1,jt2,ind,nv,idim_type,npts,ng
   integer                              :: mpiidtags,mpiidvtab
   integer, dimension(ntags)            :: msgtags
   !---------------------------------------------------------------------------------------!

   !----- Defining how many variables will be collected -----------------------------------!
   numvars=0
   do ng=1,ngrids
      do nv=1,num_var(ng)
         if (vtab_r(nv,ng)%impt3 == 1) numvars=numvars+1
      end do
   end do

   machloop: do nmi = 1,nmachs
      mpiidtags=200000000 + maxvars*maxgrds*(machnum(nmi)-1)
      mpiidvtab=300000000 + maxvars*maxgrds*(machnum(nmi)-1)
      varloop: do nvvv=1,numvars
         mpiidtags = mpiidtags + 1
         mpiidvtab = mpiidvtab + 1
         !----- Receiving the tag information and assigning the variables. ----------------!
         call MPI_Recv(msgtags,ntags,MPI_INTEGER,machnum(nmi),mpiidtags,MPI_COMM_WORLD     &
                      ,MPI_STATUS_IGNORE,ierr)
         nm        = msgtags(1)
         ng        = msgtags(2)
         npts      = msgtags(3)
         idim_type = msgtags(4)
         nv        = msgtags(5)

         !----- Copying the domain edges to scratch variables -----------------------------!
         il1       = nxbegc(nm,ng)
         ir2       = nxendc(nm,ng)
         jb1       = nybegc(nm,ng)
         jt2       = nyendc(nm,ng)
         
         !----- Receiving the data in the scratch array. ----------------------------------!
         call MPI_Recv(scratch%scr1,npts,MPI_REAL,nm,mpiidvtab,MPI_COMM_WORLD              &
                      ,MPI_STATUS_IGNORE,ierr)
         if(iand(ibcflg(nm,ng),1) /= 0) il1 = il1 - 1
         if(iand(ibcflg(nm,ng),2) /= 0) ir2 = ir2 + 1
         if(iand(ibcflg(nm,ng),4) /= 0) jb1 = jb1 - 1
         if(iand(ibcflg(nm,ng),8) /= 0) jt2 = jt2 + 1
         
         !----- Copying the subdomain dimensions to scratch variables ---------------------! 
         mxp  = nxend(nm,ng) - nxbeg(nm,ng) + 1
         myp  = nyend(nm,ng) - nybeg(nm,ng) + 1
         mxyp = mxp * myp

         !---------------------------------------------------------------------------------!
         !     Copying the scratch array to the proper place, using the var table pointer. !
         ! This is already set up so it will copy to the appropriate place.                !
         !---------------------------------------------------------------------------------!
         select case (idim_type)
         case (2) !----- 2D variables (nxp,nyp) -------------------------------------------!
            call ex_2_buff(vtab_r(nv,ng)%var_p,scratch%scr1,nnxp(ng),nnyp(ng),mxp,myp      &
                          ,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

         case (3) !----- 3D variables (nzp,nxp,nyp) ---------------------------------------!
            call ex_3_buff(vtab_r(nv,ng)%var_p,scratch%scr1,nnzp(ng),nnxp(ng),nnyp(ng)     &
                          ,mxp,myp,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

         case (4) !----- 4D variables (nzg,nxp,nyp,npatch) --------------------------------!
            call ex_4_buff(vtab_r(nv,ng)%var_p,scratch%scr1,nzg,nnxp(ng),nnyp(ng)          &
                          ,npatch,mxp,myp,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

         case (5) !----- 4D variables (nzs,nxp,nyp,npatch) --------------------------------!
            call ex_4_buff(vtab_r(nv,ng)%var_p,scratch%scr1,nzs,nnxp(ng),nnyp(ng)          &
                          ,npatch,mxp,myp,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

         case (6) !----- 3D variables (nxp,nyp,npatch) ------------------------------------!
            call ex_2p_buff(vtab_r(nv,ng)%var_p,scratch%scr1,nnxp(ng),nnyp(ng),npatch      &
                           ,mxp,myp,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

         case (7) !----- 3D variables (nxp,nyp,nwave) -------------------------------------!
            call ex_2p_buff(vtab_r(nv,ng)%var_p,scratch%scr1,nnxp(ng),nnyp(ng),nwave       &
                           ,mxp,myp,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

         case (8) !----- 4D variables (nzp,nxp,nyp,nclouds) -------------------------------!
            call ex_4_buff(vtab_r(nv,ng)%var_p,scratch%scr1,nnzp(ng),nnxp(ng),nnyp(ng)     &
                          ,nclouds,mxp,myp,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

         case (9) !----- 3D variables (nxp,nyp,nclouds) -----------------------------------!
            call ex_2p_buff(vtab_r(nv,ng)%var_p,scratch%scr1,nnxp(ng),nnyp(ng),nclouds     &
                           ,mxp,myp,ixoff(nm,ng),iyoff(nm,ng),il1,ir2,jb1,jt2)

         end select

      end do varloop
   end do machloop
   return
end subroutine master_getall
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine copies the 2D variable sub-domain data on the buffer to their actual !
! place in the full domain. Here we use the following notation: n? is the full domain      !
! dimension, whereas m? is the machine subdomain dimension.                                !
!------------------------------------------------------------------------------------------!
subroutine ex_2_buff(mydata,buff,nx,ny,mx,my,ioff,joff,ibeg,iend,jbeg,jend)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer               , intent(in)    :: nx,ny,mx,my,ioff,joff,ibeg,iend,jbeg,jend
   real, dimension(mx,my), intent(in)    :: buff
   real, dimension(nx,ny), intent(inout) :: mydata
   !----- Local variables. ----------------------------------------------------------------!
   integer                               :: ni
   integer                               :: nj
   integer                               :: mi
   integer                               :: mj
   !---------------------------------------------------------------------------------------!

   do mj=jbeg,jend
      nj = mj + joff
      do mi=ibeg,iend
         ni = mi + ioff
         mydata(ni,nj) = buff(mi,mj)
      end do
   end do

   return
end subroutine ex_2_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine copies the 2D variable (with an extra dimension that is usually      !
! patch, wavelength or cloud size) sub-domain data on the buffer to their actual           !
! place in the full domain. Here we use the following notation: n? is the full domain      !
! dimension, whereas m? is the machine subdomain dimension.                                !
!------------------------------------------------------------------------------------------!
subroutine ex_2p_buff(mydata,buff,nx,ny,ne,mx,my,ioff,joff,ibeg,iend,jbeg,jend)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                  , intent(in)    :: nx,ny,ne,mx,my,ioff,joff,ibeg,iend,jbeg,jend
   real, dimension(mx,my,ne), intent(in)    :: buff
   real, dimension(nx,ny,ne), intent(inout) :: mydata
   !----- Local variables. ----------------------------------------------------------------!
   integer                               :: ni
   integer                               :: nj
   integer                               :: nl
   integer                               :: mi
   integer                               :: mj
   !---------------------------------------------------------------------------------------!

   do nl = 1, ne
      do mj=jbeg,jend
         nj = mj + joff
         do mi=ibeg,iend
            ni = mi + ioff
            mydata(ni,nj,nl) = buff(mi,mj,nl)
         end do
      end do
   end do

   return
end subroutine ex_2p_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine copies the 3D variable sub-domain data on the buffer to their actual !
! place in the full domain. X and Y are the only dimensions allowed to have sub-domains.   !
! Here we use the following notation: n? is the full domain dimension, whereas m? is the   !
! machine subdomain dimension.                                                             !
!------------------------------------------------------------------------------------------!
subroutine ex_3_buff(mydata,buff,nz,nx,ny,mx,my,ioff,joff,ibeg,iend,jbeg,jend)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                  , intent(in)    :: nz,nx,ny,mx,my,ioff,joff,ibeg,iend,jbeg,jend
   real, dimension(nz,mx,my), intent(in)    :: buff
   real, dimension(nz,nx,ny), intent(inout) :: mydata
   !----- Local variables. ----------------------------------------------------------------!
   integer                               :: ni
   integer                               :: nj
   integer                               :: nk
   integer                               :: mi
   integer                               :: mj
   !---------------------------------------------------------------------------------------!

   do mj=jbeg,jend
      nj = mj + joff
      do mi=ibeg,iend
         ni = mi + ioff
         do nk=1,nz
            mydata(nk,ni,nj) = buff(nk,mi,mj)
         end do
      end do
   end do

   return
end subroutine ex_3_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine copies the 3D variable (with an extra dimension that is usually      !
! patch, wavelength or cloud size) sub-domain data on the buffer to their actual           !
! place in the full domain. Here we use the following notation: n? is the full domain      !
! dimension, whereas m? is the machine subdomain dimension.                                !
!------------------------------------------------------------------------------------------!
subroutine ex_4_buff(mydata,buff,nz,nx,ny,ne,mx,my,ioff,joff,ibeg,iend,jbeg,jend)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                     , intent(in)    :: nz,nx,ny,ne,mx,my,ioff,joff
   integer                     , intent(in)    :: ibeg,iend,jbeg,jend
   real, dimension(nz,mx,my,ne), intent(in)    :: buff
   real, dimension(nz,nx,ny,ne), intent(inout) :: mydata
   !----- Local variables. ----------------------------------------------------------------!
   integer                               :: ni
   integer                               :: nj
   integer                               :: nk
   integer                               :: nl
   integer                               :: mi
   integer                               :: mj
   !---------------------------------------------------------------------------------------!

   do nl=1,ne
      do mj=jbeg,jend
         nj = mj + joff
         do mi=ibeg,iend
            ni = mi + ioff
            do nk=1,nz
               mydata(nk,ni,nj,nl) = buff(nk,mi,mj,nl)
            end do
         end do
      end do
   end do

   return
end subroutine ex_4_buff
!==========================================================================================!
!==========================================================================================!
