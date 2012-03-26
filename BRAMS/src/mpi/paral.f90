!===================================== Change Log =========================================!
! 5.0.0                                                                                    !
!                                                                                          !
! MLO - 09/29/08 Including the new Grell related dimensions.                               !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This routine sends the lite and/or averaged variables (structures) back to the head   !
! node if these analyses are requested.                                                    !
!------------------------------------------------------------------------------------------!
subroutine node_sendanl(vtype)
   use node_mod
   use mem_grid
   use var_tables
   use mem_scratch
   use grid_dims


   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Arguments -----------------------------------------------------------------------!
   character(len=*)         , intent(in) :: vtype
   !----- Local constants -----------------------------------------------------------------!
   integer                  , parameter  :: ntags=5
   !----- Local variables -----------------------------------------------------------------!
   integer, dimension(ntags)             :: msgtags
   integer                               :: ierr,msgnum,nv,ng,mpiidtags,mpiidvtab
   integer                               :: mxyp,mxyzp,nmp,ind,isend
   !---------------------------------------------------------------------------------------!

   mpiidtags = 400000000+maxvars*maxgrds*(mynum-1)
   mpiidvtab = 500000000+maxvars*maxgrds*(mynum-1)
   !----- First tag is always the node ID, defining it only once --------------------------!
   msgtags(1)=mynum
   do ng=1,ngrids
      !----- Grid ID, defining only when we switch grids ----------------------------------!
      msgtags(2)=ng

      do nv=1, num_var(ng)
         !---------------------------------------------------------------------------------!
         !     Lite Analysis time. Check whether the variable should be there and, if so,  !
         ! send it to the head node.                                                       !
         !---------------------------------------------------------------------------------!
         if(trim(vtype) == 'LITE' .and. vtab_r(nv,ng)%ilite == 1) then
            mpiidtags = mpiidtags + 1
            mpiidvtab = mpiidvtab + 1

            !----- Variable specific tags -------------------------------------------------!
            msgtags(3)= vtab_r(nv,ng)%npts
            msgtags(4)= vtab_r(nv,ng)%idim_type
            msgtags(5)= nv

            !----- Sending the tags and the variable --------------------------------------!
            call MPI_Send(msgtags,ntags,MPI_INTEGER,master_num,mpiidtags,MPI_COMM_WORLD    &
                         ,ierr)
            call MPI_Send(vtab_r(nv,ng)%var_p, vtab_r(nv,ng)%npts,MPI_REAL,master_num      &
                         ,mpiidvtab,MPI_COMM_WORLD,ierr)

         !---------------------------------------------------------------------------------!
         !     Averaged analysis time. Check whether the variable should be there and, if  !
         ! so, send the averaged value to the head node.                                   !
         !---------------------------------------------------------------------------------!
         elseif( (vtype == 'MEAN' .or. vtype == 'BOTH') .and. vtab_r(nv,ng)%imean == 1)    &
         then
            mpiidtags = mpiidtags + 1
            mpiidvtab = mpiidvtab + 1

            !----- Variable specific tags -------------------------------------------------!
            msgtags(3)= vtab_r(nv,ng)%npts
            msgtags(4)= vtab_r(nv,ng)%idim_type
            msgtags(5)= nv

            !----- Sending the tags and the variable --------------------------------------!
            call MPI_Send(msgtags,ntags,MPI_INTEGER,master_num,mpiidtags,MPI_COMM_WORLD    &
                         ,ierr)
            call MPI_Send(vtab_r(nv,ng)%var_m, vtab_r(nv,ng)%npts,MPI_REAL,master_num      &
                         ,mpiidvtab,MPI_COMM_WORLD,ierr)
         end if
      end do
   end do

   return
end subroutine node_sendanl
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine retrieves the information sent by the nodes, reads at a temporary      !
! place and copy it to the appropriate variable depending on its type. This happens when   !
! a "lite", "averaged", or "both" analysis are going to be output.                         !
!------------------------------------------------------------------------------------------!
subroutine master_getanl(vtype)

   use mem_grid
   use rpara
   use mem_scratch
   use var_tables
   use mem_cuparm, only : nclouds
   use mem_aerad , only : nwave
   use grid_dims

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Arguments -----------------------------------------------------------------------!
   character(len=*)         , intent(in) :: vtype
   !----- Local constants -----------------------------------------------------------------!
   integer                  , parameter  :: ntags=5
   !----- Local variables -----------------------------------------------------------------!
   integer, dimension(ntags)             :: msgtags
   integer                               :: ierr
   integer                               :: nmi
   integer                               :: numvars
   integer                               :: ng
   integer                               :: nummess
   integer                               :: nvvv
   integer                               :: ibytes
   integer                               :: msgtyp
   integer                               :: ihostnum
   integer                               :: nm
   integer                               :: mlon
   integer                               :: mlat
   integer                               :: ioff
   integer                               :: joff
   integer                               :: iwest
   integer                               :: ieast
   integer                               :: jsouth
   integer                               :: jnorth
   integer                               :: sdxp
   integer                               :: sdyp
   integer                               :: sdxyp
   integer                               :: fdzp
   integer                               :: fdep
   integer                               :: nv
   integer                               :: npts
   integer                               :: idim_type
   integer                               :: mpiidtags
   integer                               :: mpiidvtab
   !---------------------------------------------------------------------------------------!


   !----- Check the amount of information (messages) to be received. ----------------------!
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

   !----- Receive the information. First the tags, then the actual data. ------------------!
   machloop: do nmi = 1,nmachs
      mpiidtags = 400000000 + maxvars*maxgrds*(machnum(nmi)-1)
      mpiidvtab = 500000000 + maxvars*maxgrds*(machnum(nmi)-1)
      varloop: do nvvv =1,numvars
         mpiidtags = mpiidtags + 1
         mpiidvtab = mpiidvtab + 1

         !----- Get the tags, and assigning to some scratch variables. --------------------!
         call MPI_Recv(msgtags,ntags,MPI_INTEGER,machnum(nmi),mpiidtags,MPI_COMM_WORLD     &
                      ,MPI_STATUS_IGNORE,ierr)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Copy the tags to some scalar integers.                                      !
         !---------------------------------------------------------------------------------!
         nm        = msgtags(1) !----- Number of process for the message received
         ng        = msgtags(2) !----- Number of the Grid
         npts      = msgtags(3) !----- Total number of points for the real array
         idim_type = msgtags(4) !----- Dimension of the real array
         nv        = msgtags(5) !----- Number of the variable in: vtab_r(nv,ng)
         !---------------------------------------------------------------------------------!


         !----- Find the position (plane coordinates) to fit on grid. ---------------------!
         mlon   = nxend  (nm,ng) - nxbeg  (nm,ng) + 1
         mlat   = nyend  (nm,ng) - nybeg  (nm,ng) + 1
         iwest  = nxbegc (nm,ng)
         ieast  = nxendc (nm,ng)
         jsouth = nybegc (nm,ng)
         jnorth = nyendc (nm,ng)
         ioff   = ixoff  (nm,ng)
         joff   = iyoff  (nm,ng)
         !---------------------------------------------------------------------------------!


         !----- Check whether any of the node edge is an actual domain edge. --------------!
         if(iand(ibcflg(nm,ng),1) /= 0) iwest  = iwest  - 1
         if(iand(ibcflg(nm,ng),2) /= 0) ieast  = ieast  + 1
         if(iand(ibcflg(nm,ng),4) /= 0) jsouth = jsouth - 1
         if(iand(ibcflg(nm,ng),8) /= 0) jnorth = jnorth + 1
         !---------------------------------------------------------------------------------!


         !----- Determine the sub-domain sizes --------------------------------------------!
         sdxp  = nxend(nm,ng) - nxbeg(nm,ng) + 1
         sdyp  = nyend(nm,ng) - nybeg(nm,ng) + 1
         sdxyp = sdxp * sdyp
         !---------------------------------------------------------------------------------!


         !----- Receive the data; storing first at a temporary place. ---------------------!
         call MPI_Recv(scratch%scr1(1),npts,MPI_REAL,nm,mpiidvtab,MPI_COMM_WORLD           &
                      ,MPI_STATUS_IGNORE,ierr)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Define the size of the vertical and "environmental" dimensions depending on  !
         ! the variable type.                                                              !
         !---------------------------------------------------------------------------------!
         call ze_dims(ng,idim_type,.true.,fdzp,fdep)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Copying the buffer scratch in the real array. Here we check whether it is an !
         ! instantaneous or averaged data, and use the right subroutine to copy to the     !
         ! actual variable (the vtab_r%var_? is a predefined pointer to the right place).  !
         !---------------------------------------------------------------------------------!
         select case (trim(vtype))
         !----- Instantaneous variables, use var_p pointer --------------------------------!
         case ('LITE') 
              call ex_full_buff(vtab_r(nv,ng)%var_p,scratch%scr1,npts,fdzp,nnxp(ng)        &
                               ,nnyp(ng),fdep,mlon,mlat,ioff,joff,iwest,ieast              &
                               ,jsouth,jnorth)
         !----- Averaged variables, use var_m pointer -------------------------------------!
         case ('MEAN','BOTH') 
              call ex_full_buff(vtab_r(nv,ng)%var_m,scratch%scr1,npts,fdzp,nnxp(ng)        &
                               ,nnyp(ng),fdep,mlon,mlat,ioff,joff,iwest,ieast              &
                               ,jsouth,jnorth)
         end select 
      end do varloop
   end do machloop

   return
end subroutine master_getanl
!==========================================================================================!
!==========================================================================================!
