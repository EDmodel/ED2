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
#if defined(RAMS_MPI)
   use mpi
#endif

   implicit none

   !----- External variable declaration ---------------------------------------------------!
#if defined(RAMS_MPI)
   include 'interface.h'
#endif
   !----- Local constants. ----------------------------------------------------------------!
   integer, parameter        :: ntags=9
   !----- Local variables. ----------------------------------------------------------------!
   integer, dimension(ntags) :: msgtags  ! List of variable dims to be sent to the node
   integer                   :: ierr     !
   integer                   :: nm       ! Machine number counter for machine loop
   integer                   :: ng       ! Grid number counter for grid loop
   integer                   :: nv       ! Variable number counter for variable loop
   integer                   :: sdxp     ! Number of X points for this sub-domain
   integer                   :: sdyp     ! Number of Y points for this sub-domain
   integer                   :: fdzp     ! Number of Z points for the full domain
   integer                   :: fdep     ! Number of E points for the full domain
                                         !    (clouds, patches, waves...)
   integer                   :: sdxyp    ! Number of horizontal points for this sub-domain
   integer                   :: npts     ! Number of points that this package contains
   integer                   :: mpivarid ! Flag for MPI so it won't be messed up.
   integer                   :: mpitagid ! Flag for MPI so it won't be messed up.
   integer                   :: mpivnmid ! Flag for MPI so it won't be messed up.
   !---------------------------------------------------------------------------------------!

#if defined(RAMS_MPI)
   machloop: do nm = 1,nmachs
      !----- Send the time information. ---------------------------------------------------!
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

      mpivarid = 100000
      mpitagid = 600000
      mpivnmid = 900000
      gridloop: do ng = 1,ngrids

         !----- Define sub-domain dimensions. ---------------------------------------------!
         sdxp   = nxend(nm,ng) - nxbeg(nm,ng) + 1
         sdyp   = nyend(nm,ng) - nybeg(nm,ng) + 1
         sdxyp  = sdxp * sdyp

         !----- Send variables that the nodes will only have the subdomain portion of. ----!
         varloop: do nv = 1,num_var(ng)
            if ( vtab_r(nv,ng)%impti == 1) then
               !----- This will create a unique flag for this package ---------------------!
               mpivarid = mpivarid + 1
               mpitagid = mpitagid + 1
               mpivnmid = mpivnmid + 1
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Define the size of the vertical and "environmental" dimensions depend- !
               ! ing on the variable type.                                                 !
               !---------------------------------------------------------------------------!
               call ze_dims(ng,vtab_r(nv,ng)%idim_type,.true.,fdzp,fdep)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !       The total number of points in the buffer is the product of the four !
               ! dimensions.                                                               !
               !---------------------------------------------------------------------------!
               npts = fdzp * sdxyp * fdep
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Build the vector with the tag information.                            !
               !---------------------------------------------------------------------------!
               msgtags(1) = nm
               msgtags(2) = ng
               msgtags(3) = npts
               msgtags(4) = vtab_r(nv,ng)%idim_type
               msgtags(5) = nv
               msgtags(6) = sdxp
               msgtags(7) = sdyp
               msgtags(8) = fdzp
               msgtags(9) = fdep
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Copy the information that will be sent to scratch%scr2.  This is a    !
               ! scratch variable that is large enough to receive any subdomain.           !
               !---------------------------------------------------------------------------!
               call azero(npts,scratch%scr2)
               call mk_full_buff(vtab_r(nv,ng)%var_p,scratch%scr2,npts                     &
                                ,fdzp,nnxp(ng),nnyp(ng),fdep                               &
                                ,nxbeg(nm,ng),nxend(nm,ng),nybeg(nm,ng),nyend(nm,ng))
               !---------------------------------------------------------------------------!




               !----- Send the variable name to the node. ---------------------------------!
               call MPI_Send(vtab_r(nv,ng)%name,16,MPI_CHARACTER,machnum(nm),mpivnmid      &
                            ,MPI_COMM_WORLD,ierr)
               !---------------------------------------------------------------------------!



               !----- Send the tag info to the node. --------------------------------------!
               call MPI_Send(msgtags,ntags,MPI_INTEGER,machnum(nm),mpitagid                &
                            ,MPI_COMM_WORLD,ierr)
               !---------------------------------------------------------------------------!



               !----- Send the buffer to the node -----------------------------------------!
               call MPI_Send(scratch%scr2,npts,MPI_REAL,machnum(nm),mpivarid               &
                            ,MPI_COMM_WORLD,ierr)
               !---------------------------------------------------------------------------!
            end if
         end do varloop
      end do gridloop
   end do machloop
#endif
   return
end subroutine master_sendinit
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
   use mem_aerad, only: nwave
#if defined(RAMS_MPI)
   use mpi
#endif

   implicit none

   !----- External variable declaration ---------------------------------------------------!
#if defined(RAMS_MPI)
   include 'interface.h'
#endif
   !----- Local constants. ----------------------------------------------------------------!
   integer, parameter        :: ntags=9
   !----- Local variables -----------------------------------------------------------------!
   character(len=16)         :: name_exp ! Expected variable name
   integer, dimension(ntags) :: msgtags  ! List of variable dims to be sent to the node
   integer :: ierr                       ! Error flag
   integer :: sdxp                       ! Number of points in X direction
   integer :: sdyp                       ! Number of points in Y direction
   integer :: fdzp                       ! Number of points in Z direction
   integer :: fdep                       ! Number of points in "E" direction
   integer :: nv                         ! variable number
   integer :: ng                         ! grid number
   integer :: nm                         ! machine number
   integer :: nv_exp                     ! Expected variable number
   integer :: ng_exp                     ! Expected grid number
   integer :: nm_exp                     ! Expected machine number
   integer :: npts_exp                   ! Expected number of points
   integer :: sdxp_exp                   ! Expected number of points in X direction
   integer :: sdyp_exp                   ! Expected number of points in Y direction
   integer :: fdzp_exp                   ! Expected number of points in Z direction
   integer :: fdep_exp                   ! Expected number of points in "E" direction
   integer :: idim_exp                   ! Expected number of dimensions.
   integer :: mpivarid                   ! Unique identifier for the actual variable.
   integer :: mpitagid                   ! Unique identifier for the tag.
   integer :: mpivnmid                   ! Unique identifier for the variable name.
   !---------------------------------------------------------------------------------------!

#if defined(RAMS_MPI)
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

   mpivarid = 100000
   mpitagid = 600000
   mpivnmid = 900000
   gridloop: do ng = 1,ngrids
      varloop: do nv = 1,num_var(ng)

         if (vtab_r(nv,ng)%impti == 1) then

            !----- Update the unique flag counters. ---------------------------------------!
            mpivarid = mpivarid + 1
            mpitagid = mpitagid + 1
            mpivnmid = mpivnmid + 1
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Receive the variable name.                                               !
            !------------------------------------------------------------------------------!
            call MPI_Recv(name_exp,16,MPI_CHARACTER,master_num,mpivnmid,MPI_COMM_WORLD     &
                         ,MPI_STATUS_IGNORE,ierr)
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !     Receive the tag.                                                         !
            !------------------------------------------------------------------------------!
            call MPI_Recv(msgtags,ntags,MPI_INTEGER,master_num,mpitagid,MPI_COMM_WORLD     &
                         ,MPI_STATUS_IGNORE,ierr)
            !------------------------------------------------------------------------------!




            !----- Define the dimensions of this receive. ---------------------------------!
            sdxp = mmxp(ng)
            sdyp = mmyp(ng)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Define the size of the vertical and "environmental" dimensions depend-    !
            ! ing on the variable type.                                                    !
            !------------------------------------------------------------------------------!
            call ze_dims(ng,vtab_r(nv,ng)%idim_type,.false.,fdzp,fdep)
            !------------------------------------------------------------------------------!


            !----- Retrieve the expected values. ------------------------------------------!
            nm_exp   = msgtags(1)
            ng_exp   = msgtags(2)
            npts_exp = msgtags(3)
            idim_exp = msgtags(4)
            nv_exp   = msgtags(5)
            sdxp_exp = msgtags(6)
            sdyp_exp = msgtags(7)
            fdzp_exp = msgtags(8)
            fdep_exp = msgtags(9)
            !------------------------------------------------------------------------------!


            !------ Sanity check. ---------------------------------------------------------!
            if ( vtab_r(nv,ng)%name      /= name_exp .or.                                  &
                 mynum                   /= nm_exp   .or.                                  &
                 ng                      /= ng_exp   .or.                                  &
                 nv                      /= nv_exp   .or.                                  &
                 vtab_r(nv,ng)%npts      /= npts_exp .or.                                  &
                 vtab_r(nv,ng)%idim_type /= idim_exp .or.                                  &
                 sdxp                    /= sdxp_exp .or.                                  &
                 sdyp                    /= sdyp_exp .or.                                  &
                 fdzp                    /= fdzp_exp .or.                                  &
                 fdep                    /= fdep_exp )                                     &
            then
               write (unit=*,fmt='(a)')       '---- Actual point ------------------------'
               write (unit=*,fmt='(a,1x,a)')  '  + NAME      =',trim(vtab_r(nv,ng)%name)
               write (unit=*,fmt='(a,1x,i9)') '  + NM        =',mynum
               write (unit=*,fmt='(a,1x,i9)') '  + NG        =',ng
               write (unit=*,fmt='(a,1x,i9)') '  + NV        =',nv
               write (unit=*,fmt='(a,1x,i9)') '  + NPTS      =',vtab_r(nv,nv)%npts
               write (unit=*,fmt='(a,1x,i9)') '  + IDIM_TYPE =',vtab_r(nv,nv)%idim_type
               write (unit=*,fmt='(a,1x,i9)') '  + SDXP      =',sdxp
               write (unit=*,fmt='(a,1x,i9)') '  + SDYP      =',sdyp
               write (unit=*,fmt='(a,1x,i9)') '  + FDZP      =',fdzp
               write (unit=*,fmt='(a,1x,i9)') '  + FDEP      =',fdep
               write (unit=*,fmt='(a)')       '---- Message received --------------------'
               write (unit=*,fmt='(a,1x,a)')  '  + NAME      =',trim(name_exp)
               write (unit=*,fmt='(a,1x,i9)') '  + NM        =',nm_exp
               write (unit=*,fmt='(a,1x,i9)') '  + NG        =',ng_exp
               write (unit=*,fmt='(a,1x,i9)') '  + NV        =',nv_exp
               write (unit=*,fmt='(a,1x,i9)') '  + NPTS      =',npts_exp
               write (unit=*,fmt='(a,1x,i9)') '  + IDIM_TYPE =',idim_exp
               write (unit=*,fmt='(a,1x,i9)') '  + SDXP      =',sdxp_exp
               write (unit=*,fmt='(a,1x,i9)') '  + SDYP      =',sdyp_exp
               write (unit=*,fmt='(a,1x,i9)') '  + FDZP      =',fdzp_exp
               write (unit=*,fmt='(a,1x,i9)') '  + FDEP      =',fdep_exp
               call abort_run('Conflict in the message received...','node_getinit'         &
                             ,'mpass_full.f90')
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Receive the variable sub-domain, using the variable table pointer to      !
            ! receive the information.  It is already set up to the actual place so this   !
            ! will direct it to the right place.                                           !
            !------------------------------------------------------------------------------!
            call azero(npts_exp,scratch%scr1)
            call MPI_Recv(scratch%scr1,npts_exp,MPI_REAL,master_num,mpivarid               &
                         ,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Extract the received information and place it into the right place in    !
            ! the actual variable to which the variable table points.                      !
            !------------------------------------------------------------------------------!
            call ex_full_buff(vtab_r(nv,ng)%var_p,scratch%scr1,vtab_r(nv,ng)%npts,fdzp     &
                             ,sdxp,sdyp,fdep,sdxp,sdyp,0,0,1,sdxp,1,sdyp)
            !------------------------------------------------------------------------------!
         end if

      end do varloop
   end do gridloop
#endif
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
   use mem_cuparm , only : nclouds ! ! intent(in)
   use mem_aerad  , only : nwave   ! ! intent(in)
#if defined(RAMS_MPI)
   use mpi
#endif

   implicit none
   !----- External variable declaration ---------------------------------------------------!
#if defined(RAMS_MPI)
   include 'interface.h'
#endif
   !----- Local constants -----------------------------------------------------------------!
   integer                  , parameter :: ntags = 9
   !----- Local variables -----------------------------------------------------------------!
   integer, dimension(ntags)            :: msgtags
   integer                              :: ierr
   integer                              :: mpitagsid
   integer                              :: mpivtabid
   integer                              :: mpinameid
   integer                              :: sdxp
   integer                              :: sdyp
   integer                              :: fdzp
   integer                              :: fdep
   integer                              :: nmp
   integer                              :: nv
   integer                              :: ind
   integer                              :: ng
   !---------------------------------------------------------------------------------------!

#if defined(RAMS_MPI)
   mpitagsid=200000000 + maxvars*maxgrds*(mynum-1)
   mpivtabid=300000000 + maxvars*maxgrds*(mynum-1)
   mpinameid=400000000 + maxvars*maxgrds*(mynum-1)
   gridloop: do ng=1,ngrids

      varloop: do nv=1,num_var(ng)

         if (vtab_r(nv,ng)%impt3 == 1) then
            !----- Update the tag counters. -----------------------------------------------!
            mpitagsid = mpitagsid + 1
            mpivtabid = mpivtabid + 1
            mpinameid = mpinameid + 1
            !------------------------------------------------------------------------------!



            !----- Define the horizontal dimensions of this sub-domain. -------------------!
            sdxp = mmxp(ng)
            sdyp = mmyp(ng)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Find the size of the vertical and "environmental" dimensions depending on !
            ! the variable type.                                                           !
            !------------------------------------------------------------------------------!
            call ze_dims(ng,vtab_r(nv,ng)%idim_type,.false.,fdzp,fdep)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Prepare the tags with information about the message that is about to be !
            ! sent to the head node.                                                       !
            !------------------------------------------------------------------------------!
            msgtags(1) = mynum
            msgtags(2) = ng
            msgtags(3) = vtab_r(nv,ng)%npts
            msgtags(4) = vtab_r(nv,ng)%idim_type
            msgtags(5) = nv
            msgtags(6) = sdxp
            msgtags(7) = sdyp
            msgtags(8) = fdzp
            msgtags(9) = fdep
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Copy the information that will be sent to scratch%scr6.  This is a       !
            ! scratch variable that is large enough to receive any subdomain.              !
            !------------------------------------------------------------------------------!
            call azero(vtab_r(nv,ng)%npts,scratch%scr6)
            call mk_full_buff(vtab_r(nv,ng)%var_p,scratch%scr6,vtab_r(nv,ng)%npts          &
                             ,fdzp,sdxp,sdyp,fdep,1,sdxp,1,sdyp)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------! 
            !    Send the messages to the head node.  They go in three steps, the first    !
            ! with the variable name, the second with the dimensions, and the third with   !
            ! the field.                                                                   !
            !------------------------------------------------------------------------------! 
            call MPI_Send(vtab_r(nv,ng)%name,16,MPI_CHARACTER,master_num,mpinameid         &
                         ,MPI_COMM_WORLD,ierr)
            call MPI_Send(msgtags,ntags,MPI_INTEGER,master_num,mpitagsid                   &
                         ,MPI_COMM_WORLD,ierr)
            call MPI_Send(scratch%scr6,vtab_r(nv,ng)%npts,MPI_REAL,master_num,mpivtabid    &
                         ,MPI_COMM_WORLD,ierr)
            !------------------------------------------------------------------------------! 

         end if
      end do varloop
   end do gridloop
#endif

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
   use mem_cuparm , only: nclouds ! ! intent(in)
   use mem_aerad  , only: nwave   ! ! intent(in)
   use grid_dims
#if defined(RAMS_MPI)
   use mpi
#endif
  
   implicit none
   !----- External variable declaration ---------------------------------------------------!
#if defined(RAMS_MPI)
   include 'interface.h'
#endif
   !----- Local constants -----------------------------------------------------------------!
   integer                  , parameter :: ntags=9
   !----- Local variables -----------------------------------------------------------------!
   character(len=16)                    :: name_exp
   integer, dimension(ntags)            :: msgtags
   integer                              :: ierr
   integer                              :: nmi
   integer                              :: numvars
   integer                              :: nvvv
   integer                              :: sdxp
   integer                              :: sdyp
   integer                              :: sdxyp
   integer                              :: fdzp
   integer                              :: fdep
   integer                              :: sdxp_exp
   integer                              :: sdyp_exp
   integer                              :: fdzp_exp
   integer                              :: fdep_exp
   integer                              :: npts_exp
   integer                              :: nm
   integer                              :: iwest
   integer                              :: ieast
   integer                              :: jsouth
   integer                              :: jnorth
   integer                              :: ioff
   integer                              :: joff
   integer                              :: nv
   integer                              :: idim_type
   integer                              :: npts
   integer                              :: ng
   integer                              :: mpitagsid
   integer                              :: mpivtabid
   integer                              :: mpinameid
   !---------------------------------------------------------------------------------------!

   !----- Determine how many variables will be collected ---------------------------------!
   numvars=0
   do ng=1,ngrids
      do nv=1,num_var(ng)
         if (vtab_r(nv,ng)%impt3 == 1) numvars=numvars+1
      end do
   end do
   !---------------------------------------------------------------------------------------!


#if defined(RAMS_MPI)

   !---------------------------------------------------------------------------------------!
   !     Loop over all machines to receive the nodes information.                          !
   !---------------------------------------------------------------------------------------!
   machloop: do nmi = 1,nmachs
      mpitagsid = 200000000 + maxvars*maxgrds*(machnum(nmi)-1)
      mpivtabid = 300000000 + maxvars*maxgrds*(machnum(nmi)-1)
      mpinameid = 400000000 + maxvars*maxgrds*(machnum(nmi)-1)

      varloop: do nvvv=1,numvars
         !----- Update the message flag. --------------------------------------------------!
         mpitagsid = mpitagsid + 1
         mpivtabid = mpivtabid + 1
         mpinameid = mpinameid + 1
         !---------------------------------------------------------------------------------!



         !----- Receive the expected variable name. ---------------------------------------!
         call MPI_Recv(name_exp,16,MPI_CHARACTER,machnum(nmi),mpinameid,MPI_COMM_WORLD     &
                      ,MPI_STATUS_IGNORE,ierr)
         !---------------------------------------------------------------------------------!


         !----- Receive the tag information and assign some variables. --------------------!
         call MPI_Recv(msgtags,ntags,MPI_INTEGER,machnum(nmi),mpitagsid,MPI_COMM_WORLD     &
                      ,MPI_STATUS_IGNORE,ierr)
         nm        = msgtags(1)
         ng        = msgtags(2)
         npts_exp  = msgtags(3)
         idim_type = msgtags(4)
         nv        = msgtags(5)
         sdxp_exp  = msgtags(6)
         sdyp_exp  = msgtags(7)
         fdzp_exp  = msgtags(8)
         fdep_exp  = msgtags(9)
         !---------------------------------------------------------------------------------!



         !----- Copy the domain edges to scratch variables. -------------------------------!
         sdxp   = nxend  (nm,ng) - nxbeg  (nm,ng) + 1
         sdyp   = nyend  (nm,ng) - nybeg  (nm,ng) + 1
         iwest  = nxbegc (nm,ng)
         ieast  = nxendc (nm,ng)
         jsouth = nybegc (nm,ng)
         jnorth = nyendc (nm,ng)
         ioff   = ixoff  (nm,ng)
         joff   = iyoff  (nm,ng)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the vertical and "environmental" dimensions depending on the variable  !
         ! type.                                                                           !
         !---------------------------------------------------------------------------------!
         call ze_dims(ng,idim_type,.true.,fdzp,fdep)
         !---------------------------------------------------------------------------------!



         !----- Determine the number of points. -------------------------------------------!
         npts = sdxp * sdyp * fdzp * fdep
         !---------------------------------------------------------------------------------!



         !------ Sanity check. ------------------------------------------------------------!
         if ( vtab_r(nv,ng)%name      /= name_exp     .or.                                 &
              nm                      /= machnum(nmi) .or.                                 &
              npts                    /= npts_exp     .or.                                 &
              vtab_r(nv,ng)%idim_type /= idim_type    .or.                                 &
              sdxp                    /= sdxp_exp     .or.                                 &
              sdyp                    /= sdyp_exp     .or.                                 &
              fdzp                    /= fdzp_exp     .or.                                 &
              fdep                    /= fdep_exp     )                                    &
         then
            write (unit=*,fmt='(a)')       '---- Actual point ------------------------'
            write (unit=*,fmt='(a,1x,a)')  '  + NAME      =',trim(vtab_r(nv,ng)%name)
            write (unit=*,fmt='(a,1x,i9)') '  + NM        =',machnum(nmi)
            write (unit=*,fmt='(a,1x,i9)') '  + NG        =',ng
            write (unit=*,fmt='(a,1x,i9)') '  + NV        =',nv
            write (unit=*,fmt='(a,1x,i9)') '  + NPTS      =',npts
            write (unit=*,fmt='(a,1x,i9)') '  + IDIM_TYPE =',vtab_r(nv,nv)%idim_type
            write (unit=*,fmt='(a,1x,i9)') '  + SDXP      =',sdxp
            write (unit=*,fmt='(a,1x,i9)') '  + SDYP      =',sdyp
            write (unit=*,fmt='(a,1x,i9)') '  + FDZP      =',fdzp
            write (unit=*,fmt='(a,1x,i9)') '  + FDEP      =',fdep
            write (unit=*,fmt='(a)')       '---- Message received --------------------'
            write (unit=*,fmt='(a,1x,a)')  '  + NAME      =',trim(name_exp)
            write (unit=*,fmt='(a,1x,i9)') '  + NM        =',nm
            write (unit=*,fmt='(a,1x,i9)') '  + NG        =',ng
            write (unit=*,fmt='(a,1x,i9)') '  + NV        =',nv
            write (unit=*,fmt='(a,1x,i9)') '  + NPTS      =',npts_exp
            write (unit=*,fmt='(a,1x,i9)') '  + IDIM_TYPE =',idim_type
            write (unit=*,fmt='(a,1x,i9)') '  + SDXP      =',sdxp_exp
            write (unit=*,fmt='(a,1x,i9)') '  + SDYP      =',sdyp_exp
            write (unit=*,fmt='(a,1x,i9)') '  + FDZP      =',fdzp_exp
            write (unit=*,fmt='(a,1x,i9)') '  + FDEP      =',fdep_exp
            call abort_run('Conflict in the message received...','master_getall'           &
                          ,'mpass_full.f90')
         end if
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Shift the usable domain to the actual edge for edge and corner sub-domains. !
         !---------------------------------------------------------------------------------!
         if(iand(ibcflg(nm,ng),1) /= 0) iwest  = iwest  - 1
         if(iand(ibcflg(nm,ng),2) /= 0) ieast  = ieast  + 1
         if(iand(ibcflg(nm,ng),4) /= 0) jsouth = jsouth - 1
         if(iand(ibcflg(nm,ng),8) /= 0) jnorth = jnorth + 1
         !---------------------------------------------------------------------------------!




         !----- Receive the data and put it temporarily into the scratch array. -----------!
         call azero(npts,scratch%scr1)
         call MPI_Recv(scratch%scr1,npts,MPI_REAL,nm,mpivtabid,MPI_COMM_WORLD              &
                      ,MPI_STATUS_IGNORE,ierr)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Extract the received information and place it into the right place in the   !
         ! actual variable to which the variable table points.                             !
         !---------------------------------------------------------------------------------!
         call ex_full_buff(vtab_r(nv,ng)%var_p,scratch%scr1,npts,fdzp,nnxp(ng),nnyp(ng)    &
                          ,fdep,sdxp,sdyp,ioff,joff,iwest,ieast,jsouth,jnorth)
         !---------------------------------------------------------------------------------!
      end do varloop
   end do machloop
#endif
   return
end subroutine master_getall
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine copies part of the full domain to the buffer matrix, which will be   !
! sent to the nodes.                                                                       !
!------------------------------------------------------------------------------------------!
subroutine mk_full_buff(mydata,buff,npts,nz,nx,ny,ne,iwest,ieast,jsouth,jnorth)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: ne
   integer                     , intent(in)    :: npts
   integer                     , intent(in)    :: iwest
   integer                     , intent(in)    :: ieast
   integer                     , intent(in)    :: jsouth
   integer                     , intent(in)    :: jnorth
   real, dimension(nz,nx,ny,ne), intent(in)    :: mydata
   real, dimension(npts)       , intent(inout) :: buff
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: ind
   integer                                     :: i
   integer                                     :: j
   integer                                     :: k
   integer                                     :: l
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Copy the information to the buffer.                                               !
   !---------------------------------------------------------------------------------------!
   ind = 0
   do l=1,ne
      do j=jsouth,jnorth
         do i=iwest,ieast
            do k=1,nz
               ind       = ind + 1
               buff(ind) = mydata(k,i,j,l)
            end do
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Crash in case there is a mismatch between the number of points in copied and the !
   ! dimension of the buffer.                                                              !
   !---------------------------------------------------------------------------------------!
   if (ind /= npts) then
      write (unit=*,fmt='(a)') '----------------------------------------------------------'
      write (unit=*,fmt='(a)') ' Mismatch between expected buffer size and amount copied! '
      write (unit=*,fmt='(a)') '----------------------------------------------------------'
      write (unit=*,fmt='(a,1x,i6)') ' NZ     = ',nz
      write (unit=*,fmt='(a,1x,i6)') ' NX     = ',nx
      write (unit=*,fmt='(a,1x,i6)') ' NY     = ',ny
      write (unit=*,fmt='(a,1x,i6)') ' NE     = ',ne
      write (unit=*,fmt='(a,1x,i6)') ' IWEST  = ',iwest
      write (unit=*,fmt='(a,1x,i6)') ' IEAST  = ',ieast
      write (unit=*,fmt='(a,1x,i6)') ' JSOUTH = ',jsouth
      write (unit=*,fmt='(a,1x,i6)') ' JNORTH = ',jnorth
      write (unit=*,fmt='(a,1x,i6)') ' IND    = ',ind
      write (unit=*,fmt='(a,1x,i6)') ' NPTS   = ',npts
      write (unit=*,fmt='(a)') '----------------------------------------------------------'
      call abort_run('Incorrect buffer size','mk_full_buff','mpass_full.f90')
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine mk_full_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine extracts the information stored in the buffer, and copies to the     !
! actual place that will be used by the head node.                                         !
!------------------------------------------------------------------------------------------!
subroutine ex_full_buff(mydata,buff,npts,nz,nx,ny,ne,mx,my,ioff,joff,iwest,ieast           &
                       ,jsouth,jnorth)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                     , intent(in)    :: npts
   integer                     , intent(in)    :: ne
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: mx
   integer                     , intent(in)    :: my
   integer                     , intent(in)    :: ioff
   integer                     , intent(in)    :: joff
   integer                     , intent(in)    :: iwest
   integer                     , intent(in)    :: ieast
   integer                     , intent(in)    :: jsouth
   integer                     , intent(in)    :: jnorth
   real, dimension(npts)       , intent(in)    :: buff
   real, dimension(nz,nx,ny,ne), intent(inout) :: mydata
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: iabs
   integer                                     :: jabs
   integer                                     :: i
   integer                                     :: j
   integer                                     :: k
   integer                                     :: l
   integer                                     :: ind
   logical                                     :: icopy
   logical                                     :: jcopy
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Copy the buffer information.                                                       !
   !---------------------------------------------------------------------------------------!
   !----- ind is the buffer index, notice that not everything will be copied. -------------!
   ind = 0
   eloop: do l=1,ne
      yloop: do j=1,my
         !----- Jabs is the absolute Y position. ------------------------------------------!
         jabs = j + joff
         !----- Jcopy is a flag telling whether this Y column is a resolved grid point. ---!
         jcopy = j >= jsouth .and. j <= jnorth
         xloop: do i=1,mx
            !----- Iabs is the absolute X position. ---------------------------------------!
            iabs  = i + ioff
            !------------------------------------------------------------------------------!
            !     Sanity check.                                                            !
            !------------------------------------------------------------------------------!
            if (iabs < 1 .or. iabs > nx .or. jabs < 1 .or. jabs > ny) then
               write (unit=*,fmt='(a)') '-------------------------------------------------'
               write (unit=*,fmt='(a)') ' Point overboard!!!'
               write (unit=*,fmt='(a)') '-------------------------------------------------'
               write (unit=*,fmt='(a,1x,i6)') ' I      = ',i
               write (unit=*,fmt='(a,1x,i6)') ' IABS   = ',iabs
               write (unit=*,fmt='(a,1x,i6)') ' J      = ',j
               write (unit=*,fmt='(a,1x,i6)') ' JABS   = ',jabs
               write (unit=*,fmt='(a,1x,i6)') ' NZ     = ',nz
               write (unit=*,fmt='(a,1x,i6)') ' NX     = ',nx
               write (unit=*,fmt='(a,1x,i6)') ' NY     = ',ny
               write (unit=*,fmt='(a,1x,i6)') ' NE     = ',ne
               write (unit=*,fmt='(a,1x,i6)') ' MX     = ',mx
               write (unit=*,fmt='(a,1x,i6)') ' MY     = ',my
               write (unit=*,fmt='(a,1x,i6)') ' IOFF   = ',ioff
               write (unit=*,fmt='(a,1x,i6)') ' IWEST  = ',iwest
               write (unit=*,fmt='(a,1x,i6)') ' IEAST  = ',ieast
               write (unit=*,fmt='(a,1x,i6)') ' JOFF   = ',joff
               write (unit=*,fmt='(a,1x,i6)') ' JSOUTH = ',jsouth
               write (unit=*,fmt='(a,1x,i6)') ' JNORTH = ',jnorth
               write (unit=*,fmt='(a,1x,i6)') ' IND    = ',ind
               write (unit=*,fmt='(a,1x,i6)') ' NPTS   = ',npts
               write (unit=*,fmt='(a)') '-------------------------------------------------'
               call abort_run('Grid point offset is wrong','ex_full_buff','mpass_full.f90')
            end if
            !------------------------------------------------------------------------------!


            !----- Same as jcopy but for the X rows. --------------------------------------!
            icopy = i >= iwest .and. i <= ieast

            zloop: do k=1,nz
               !----- We always update ind, but we do not always copy to mydata. ----------!
               ind = ind + 1

               !----- Check whether this is copied or not. --------------------------------!
               if (icopy .and. jcopy) then
                  mydata(k,iabs,jabs,l) = buff(ind)
               end if
            end do zloop
         end do xloop
      end do yloop
   end do eloop
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Sanity check.                                                                     !
   !---------------------------------------------------------------------------------------!
   if (ind /= npts) then
      write (unit=*,fmt='(a)') '----------------------------------------------------------'
      write (unit=*,fmt='(a)') ' Mismatch between expected buffer size and amount copied! '
      write (unit=*,fmt='(a)') '----------------------------------------------------------'
      write (unit=*,fmt='(a,1x,i6)') ' NZ     = ',nz
      write (unit=*,fmt='(a,1x,i6)') ' NX     = ',nx
      write (unit=*,fmt='(a,1x,i6)') ' NY     = ',ny
      write (unit=*,fmt='(a,1x,i6)') ' NE     = ',ne
      write (unit=*,fmt='(a,1x,i6)') ' MX     = ',mx
      write (unit=*,fmt='(a,1x,i6)') ' MY     = ',my
      write (unit=*,fmt='(a,1x,i6)') ' IOFF   = ',ioff
      write (unit=*,fmt='(a,1x,i6)') ' IWEST  = ',iwest
      write (unit=*,fmt='(a,1x,i6)') ' IEAST  = ',ieast
      write (unit=*,fmt='(a,1x,i6)') ' JOFF   = ',joff
      write (unit=*,fmt='(a,1x,i6)') ' JSOUTH = ',jsouth
      write (unit=*,fmt='(a,1x,i6)') ' JNORTH = ',jnorth
      write (unit=*,fmt='(a,1x,i6)') ' IND    = ',ind
      write (unit=*,fmt='(a,1x,i6)') ' NPTS   = ',npts
      write (unit=*,fmt='(a)') '----------------------------------------------------------'
      call abort_run('Incorrect buffer size','ex_full_buff','mpass_full.f90')
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine ex_full_buff
!==========================================================================================!
!==========================================================================================!
