!==========================================================================================!
!==========================================================================================!
! Copyright (C) 1991-2004  ; All Rights Reserved ; Colorado State University               !
! Colorado State University Research Foundation ; ATMET, LLC                               !
!                                                                                          !
! This file is free software; you can redistribute it and/or modify it under the           !
! terms of the GNU General Public License as published by the Free Software                !
! Foundation; either version 2 of the License, or (at your option) any later version.      !
!                                                                                          !
! This software is distributed in the hope that it will be useful, but WITHOUT ANY         !
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A          !
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.                !
!                                                                                          !
! You should have received a copy of the GNU General Public License along with this        !
! program; if not, write to the Free Software Foundation, Inc.,                            !
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.                                 !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine will take the boundary conditions form one node to the other.       !
!------------------------------------------------------------------------------------------!
subroutine node_sendlbc()

   use mem_grid
   use node_mod

   use var_tables
   use mem_scratch
   use mem_cuparm , only : nclouds  ! ! intent(in)
   use grid_dims  , only : maxgrds  ! ! intent(in)
   use mem_aerad  , only : nwave    ! ! intent(in)

   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer               :: ierr
   integer               :: ipos
   integer               :: itype
   integer               :: nm
   integer               :: i1
   integer               :: i2
   integer               :: j1
   integer               :: j2
   integer               :: fdzp
   integer               :: fdep
   integer               :: nptsxy
   integer               :: nv
   integer               :: mtp
   integer               :: mpiid
   !----- Module variables. ---------------------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !---------------------------------------------------------------------------------------!


   itype=1

   !---------------------------------------------------------------------------------------!
   !     First, before we send anything, let's post the receives.  Also, make sure any     !
   ! pending sends are complete.                                                           !
   !---------------------------------------------------------------------------------------!
   do nm=1,nmachs
      if (iget_paths(itype,ngrid,nm) /= 0) then
         !----- Find the unique MPI flag to make sure the right message is got. -----------!
         mpiid = 300000 + maxgrds * (nm-1) + ngrid

         !----- Post the receive. ---------------------------------------------------------!
         call MPI_Irecv(node_buffs_lbc(nm)%pack_recv_buff(1)                               &
                       ,node_buffs_lbc(nm)%nrecv*f_ndmd_size                               &
                       ,MPI_PACKED,machs(nm),mpiid,MPI_COMM_WORLD,irecv_req(nm),ierr )
      end if
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Now we can actually go on to sending the stuff.                                   !
   !---------------------------------------------------------------------------------------!
   do nm=1,nmachs

      if (ipaths(1,itype,ngrid,nm) /= 0) then

         i1 = ipaths(1,itype,ngrid,nm)
         i2 = ipaths(2,itype,ngrid,nm)
         j1 = ipaths(3,itype,ngrid,nm)
         j2 = ipaths(4,itype,ngrid,nm)
         
         nptsxy = (i2-i1+1) * (j2-j1+1)

         ipos = 1
         call MPI_Pack(i1,1,MPI_INTEGER,node_buffs_lbc(nm)%pack_send_buff                  &
                      ,node_buffs_lbc(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(i2,1,MPI_INTEGER,node_buffs_lbc(nm)%pack_send_buff                  &
                      ,node_buffs_lbc(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(j1,1,MPI_INTEGER,node_buffs_lbc(nm)%pack_send_buff                  &
                      ,node_buffs_lbc(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(j2,1,MPI_INTEGER,node_buffs_lbc(nm)%pack_send_buff                  &
                      ,node_buffs_lbc(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(mynum,1,MPI_INTEGER,node_buffs_lbc(nm)%pack_send_buff               &
                      ,node_buffs_lbc(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)

         do nv = 1,num_var(ngrid)
            if ( vtab_r(nv,ngrid)%impt1 == 1) then

               !----- Find the size of the vertical and environmental axes. ---------------!
               call ze_dims(ngrid,vtab_r(nv,ngrid)%idim_type,.false.,fdzp,fdep)
               !---------------------------------------------------------------------------!


               !----- Find the size of the lateral boundary condition. --------------------!
               mtp = nptsxy * fdzp * fdep

               !----- Copy the variable to the scratch array. -----------------------------!
               call mk_lbc_buff(fdzp,mxp,myp,fdep,mtp,vtab_r(nv,ngrid)%var_p,scratch%scr5  &
                               ,i1-i0,i2-i0,j1-j0,j2-j0)
               !---------------------------------------------------------------------------!

               !----- Add the boundary condition to the MPI buffer. -----------------------!
               call MPI_Pack(scratch%scr5 ,mtp,MPI_REAL,node_buffs_lbc(nm)%pack_send_buff  &
                            ,node_buffs_lbc(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
               !---------------------------------------------------------------------------!
            end if
         end do





         !----- Send out the stuff to the node. -------------------------------------------!
         mpiid = 300000 + maxgrds * (mynum-1) + ngrid
         call MPI_Isend(node_buffs_lbc(nm)%pack_send_buff,ipos-1,MPI_Packed                     &
                       ,ipaths(5,itype,ngrid,nm),mpiid,MPI_COMM_WORLD,isend_req(nm),ierr)
         !---------------------------------------------------------------------------------!
      end if
   end do

   return
end subroutine node_sendlbc
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine gets the lateral boundary condition information and copy to the     !
! edges of the node sub-domain.                                                            !
!------------------------------------------------------------------------------------------!
subroutine node_getlbc()
   use mem_grid
   use node_mod

   use var_tables
   use mem_scratch
   use mem_cuparm , only : nclouds  ! ! intent(in)
   use mem_aerad  , only : nwave    ! ! intent(in)

   implicit none
   !----- Module variables. ---------------------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Local variables. ----------------------------------------------------------------!
   integer, dimension(MPI_STATUS_SIZE) :: status1
   integer                             :: ierr
   integer                             :: ipos
   integer                             :: itype
   integer                             :: nm
   integer                             :: ibytes
   integer                             :: msgid
   integer                             :: ihostnum
   integer                             :: i1
   integer                             :: i2
   integer                             :: j1
   integer                             :: j2
   integer                             :: fdzp
   integer                             :: fdep
   integer                             :: nmp
   integer                             :: nv
   integer                             :: node_src
   integer                             :: mtc
   integer                             :: mtp
   integer                             :: nptsxy
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Set itype to 1 (I think this controls the source/destination of the information,  !
   ! but I am not sure.                                                                    !
   !---------------------------------------------------------------------------------------!
   itype = 1
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     First, let's make sure our sends are all finished and de-allocated.               !
   !---------------------------------------------------------------------------------------!
   do nm=1,nmachs
      if (ipaths(1,itype,ngrid,nm) /= 0) call MPI_Wait(isend_req(nm),status1,ierr)
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Now, let's wait on our receives.                                                  !
   !---------------------------------------------------------------------------------------!
   do nm=1,nmachs
      if (iget_paths(itype,ngrid,nm) /= 0) call MPI_Wait(irecv_req(nm),status1,ierr)
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      We got all our stuff.  Now unpack it into the appropriate space.                 !
   !---------------------------------------------------------------------------------------!
   machloop: do nm=1,nmachs
      if (iget_paths(itype,ngrid,nm) /= 0) then

         ipos = 1
         call MPI_Unpack(node_buffs_lbc(nm)%pack_recv_buff                                 &
                        ,node_buffs_lbc(nm)%nrecv*f_ndmd_size                              &
                        ,ipos,i1,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_lbc(nm)%pack_recv_buff                                 &
                        ,node_buffs_lbc(nm)%nrecv*f_ndmd_size                              &
                        ,ipos,i2,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_lbc(nm)%pack_recv_buff                                 &
                        ,node_buffs_lbc(nm)%nrecv*f_ndmd_size                              &
                        ,ipos,j1,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_lbc(nm)%pack_recv_buff                                 &
                        ,node_buffs_lbc(nm)%nrecv*f_ndmd_size                              &
                        ,ipos,j2,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_lbc(nm)%pack_recv_buff                                 &
                        ,node_buffs_lbc(nm)%nrecv*f_ndmd_size                              &
                        ,ipos,node_src,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

         !----- Total number of horizontal points. ----------------------------------------!
         nptsxy=(i2-i1+1)*(j2-j1+1)

         varloop: do nv = 1,num_var(ngrid)
            if (vtab_r(nv,ngrid)%impt1 == 1) then

               !----- Find the size of the vertical and environmental axes. ---------------!
               call ze_dims(ngrid,vtab_r(nv,ngrid)%idim_type,.false.,fdzp,fdep)
               !---------------------------------------------------------------------------!


               !----- Find the total number of points. ------------------------------------!
               mtp= fdzp * nptsxy * fdep
               !---------------------------------------------------------------------------!


               !----- Unpack the buffer into the scratch array. ---------------------------!
               call MPI_Unpack(node_buffs_lbc(nm)%pack_recv_buff                           &
                              ,node_buffs_lbc(nm)%nrecv*f_ndmd_size,ipos,scratch%scr6 ,mtp &
                              ,MPI_REAL,MPI_COMM_WORLD,ierr)
               !---------------------------------------------------------------------------!


               !----- Extract the information. --------------------------------------------!
               call ex_lbc_buff(fdzp,mxp,myp,fdep,mtp,vtab_r(nv,ngrid)%var_p,scratch%scr6  &
                               ,i1-i0,i2-i0,j1-j0,j2-j0)
               !---------------------------------------------------------------------------!

            end if
         end do varloop
      end if
   end do machloop

   return
end subroutine node_getlbc
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine simply copies a stretch of the sub-domain data into a buffer vector. !
!  X and Y (n2 and n3) are the only dimensions allowed to have sub-domains, and n1 and n4  !
! should always match the original dimension.                                              !
!------------------------------------------------------------------------------------------!
subroutine mk_lbc_buff(nz,nx,ny,ne,nlbc,arr,buff,iwest,ieast,jsouth,jnorth)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: ne
   integer                     , intent(in)    :: nlbc
   integer                     , intent(in)    :: iwest
   integer                     , intent(in)    :: ieast
   integer                     , intent(in)    :: jsouth
   integer                     , intent(in)    :: jnorth
   real, dimension(nz,nx,ny,ne), intent(in)    :: arr
   real, dimension(nlbc)       , intent(out)   :: buff
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: ind
   integer                                     :: i
   integer                                     :: j
   integer                                     :: k
   integer                                     :: l
   !---------------------------------------------------------------------------------------!
   
   ind = 0
   do l=1,ne
      do j=jsouth,jnorth
         do i=iwest,ieast
            do k=1,nz
              ind=ind+1
              buff(ind) = arr(k,i,j,l)
            end do
         end do
      end do
   end do

   !----- Sanity check. -------------------------------------------------------------------!
   if (ind /= nlbc) then
      write (unit=*,fmt='(a)')       '------------------------------------------------'
      write (unit=*,fmt='(a)')       ' Mismatch betweeen buffer and LBC domain sizes.'
      write (unit=*,fmt='(a)')       ' '
      write (unit=*,fmt='(a,1x,i6)') ' - NLBC   = ',nlbc
      write (unit=*,fmt='(a,1x,i6)') ' - NBUFF  = ',ind
      write (unit=*,fmt='(a,1x,i6)') ' - IWEST  = ',iwest
      write (unit=*,fmt='(a,1x,i6)') ' - IEAST  = ',ieast
      write (unit=*,fmt='(a,1x,i6)') ' - JSOUTH = ',jsouth
      write (unit=*,fmt='(a,1x,i6)') ' - JNORTH = ',jnorth
      write (unit=*,fmt='(a,1x,i6)') ' - NZP    = ',nz
      write (unit=*,fmt='(a,1x,i6)') ' - NEP    = ',ne
      write (unit=*,fmt='(a)')       ' '
      write (unit=*,fmt='(a)')       '------------------------------------------------'
      write (unit=*,fmt='(a)')       ' '
      call abort_run('Mismatch in the buffer size!','mk_lbc_buff','mpass_lbc.f90')
   end if

   return
end subroutine mk_lbc_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine simply copies the buffer vector back to a 3D subdomain array.  X and !
! Y (n2 and n3) are the only dimensions allowed to have sub-domains, and n1 and n4         !
! should always match the original dimension.                                              !
!------------------------------------------------------------------------------------------!
subroutine ex_lbc_buff(nz,nx,ny,ne,nlbc,arr,buff,iwest,ieast,jsouth,jnorth)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: ne
   integer                     , intent(in)    :: nlbc
   integer                     , intent(in)    :: iwest
   integer                     , intent(in)    :: ieast
   integer                     , intent(in)    :: jsouth
   integer                     , intent(in)    :: jnorth
   real, dimension(nz,nx,ny,ne), intent(inout) :: arr
   real, dimension(nlbc)       , intent(in)    :: buff
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: ind
   integer                                     :: i
   integer                                     :: j
   integer                                     :: k
   integer                                     :: l
   !---------------------------------------------------------------------------------------!

   ind=0

   do l=1,ne
      do j=jsouth,jnorth
         do i=iwest,ieast
            do k=1,nz
               ind = ind + 1
               arr(k,i,j,l) = buff(ind)
            end do
         end do
      end do
   end do

   !----- Sanity check. -------------------------------------------------------------------!
   if (ind /= nlbc) then
      write (unit=*,fmt='(a)')       '------------------------------------------------'
      write (unit=*,fmt='(a)')       ' Mismatch betweeen buffer and LBC domain sizes.'
      write (unit=*,fmt='(a)')       ' '
      write (unit=*,fmt='(a,1x,i6)') ' - NLBC   = ',nlbc
      write (unit=*,fmt='(a,1x,i6)') ' - NBUFF  = ',ind
      write (unit=*,fmt='(a,1x,i6)') ' - IWEST  = ',iwest
      write (unit=*,fmt='(a,1x,i6)') ' - IEAST  = ',ieast
      write (unit=*,fmt='(a,1x,i6)') ' - JSOUTH = ',jsouth
      write (unit=*,fmt='(a,1x,i6)') ' - JNORTH = ',jnorth
      write (unit=*,fmt='(a,1x,i6)') ' - NZP    = ',nz
      write (unit=*,fmt='(a,1x,i6)') ' - NEP    = ',ne
      write (unit=*,fmt='(a)')       ' '
      write (unit=*,fmt='(a)')       '------------------------------------------------'
      write (unit=*,fmt='(a)')       ' '
      call abort_run('Mismatch in the buffer size!!!','ex_lbc_buff','mpass_lbc.f90')
   end if

   return
end subroutine ex_lbc_buff
!==========================================================================================!
!==========================================================================================!
