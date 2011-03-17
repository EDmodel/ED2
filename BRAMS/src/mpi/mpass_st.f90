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
!     This sub-routine sends the boundary conditions of the staggered variables to the     !
! nodes that are solving the neighbour grid points.                                        !
!------------------------------------------------------------------------------------------!
subroutine node_sendst(isflag)

   use mem_grid
   use grid_dims
   use node_mod

   use mem_scratch
   use mem_basic

   implicit none

   !----- Module variables. ---------------------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: isflag
   !----- Local variables -----------------------------------------------------------------!
   integer             :: ierr
   integer             :: ipos
   integer             :: itype
   integer             :: nm
   integer             :: i1
   integer             :: i2
   integer             :: j1
   integer             :: j2
   integer             :: mtp
   integer             :: mpiid
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Define the type of variable based on the input flag.                              !
   !---------------------------------------------------------------------------------------!
   select case (isflag)
   case (2:4)
      itype = isflag
   case (5:6)
      itype = 1
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      First, before we send anything, let's post the receives.                         !
   !---------------------------------------------------------------------------------------!
   do nm=1,nmachs
      if (iget_paths(itype,ngrid,nm) /= 0) then
         mpiid= 2000000 + 10*maxgrds*(machs(nm)-1) + 10*(ngrid-1) + isflag
         call MPI_Irecv(node_buffs_st(nm)%pack_recv_buff                                   &
                       ,node_buffs_st(nm)%nrecv*f_ndmd_size,MPI_PACKED                     &
                       ,machs(nm),mpiid,MPI_COMM_WORLD,irecv_req(nm),ierr )
      end if
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Now we can actually go on to sending the stuff.                                  !
   !---------------------------------------------------------------------------------------!
   do nm=1,nmachs

      if (ipaths(1,itype,ngrid,nm) /= 0) then

         i1=ipaths(1,itype,ngrid,nm)
         i2=ipaths(2,itype,ngrid,nm)
         j1=ipaths(3,itype,ngrid,nm)
         j2=ipaths(4,itype,ngrid,nm)

         !----- In this case all variables are 3-D, so compute the number of points here. -!
         mtp = (i2-i1+1) * (j2-j1+1) * mzp

         !----- Include the size information in the buffer. -------------------------------!
         ipos = 1
         call MPI_Pack(i1,1,MPI_INTEGER,node_buffs_st(nm)%pack_send_buff                   &
                      ,node_buffs_st(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(i2,1,MPI_INTEGER,node_buffs_st(nm)%pack_send_buff                   &
                      ,node_buffs_st(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(j1,1,MPI_INTEGER,node_buffs_st(nm)%pack_send_buff                   &
                      ,node_buffs_st(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(j2,1,MPI_INTEGER,node_buffs_st(nm)%pack_send_buff                   &
                      ,node_buffs_st(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(mynum,1,MPI_INTEGER,node_buffs_st(nm)%pack_send_buff &
                      ,node_buffs_st(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
         !---------------------------------------------------------------------------------!


         !------ Decide which variable to include in the package based on the input flag. -!
         select case (isflag)
         case (2)
            !----- Zonal wind. ------------------------------------------------------------!
            call mk_st_buff(mzp,mxp,myp,1,mtp,basic_g(ngrid)%up,scratch%scr3               &
                           ,i1-i0,i2-i0,j1-j0,j2-j0)
            call MPI_Pack(scratch%scr3,mtp,MPI_REAL,node_buffs_st(nm)%pack_send_buff       &
                         ,node_buffs_st(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)

         case (3)
            !----- Meridional wind. -------------------------------------------------------!
            call mk_st_buff(mzp,mxp,myp,1,mtp,basic_g(ngrid)%vp,scratch%scr3               &
                           ,i1-i0,i2-i0,j1-j0,j2-j0)
            call MPI_Pack(scratch%scr3,mtp,MPI_REAL,node_buffs_st(nm)%pack_send_buff       &
                         ,node_buffs_st(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)

         case (4)
            !----- Perturbation of the Exner function. ------------------------------------!
            call mk_st_buff(mzp,mxp,myp,1,mtp,basic_g(ngrid)%pp,scratch%scr3               &
                           ,i1-i0,i2-i0,j1-j0,j2-j0)
            call MPI_Pack(scratch%scr3,mtp,MPI_REAL,node_buffs_st(nm)%pack_send_buff       &
                         ,node_buffs_st(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)

         case (5)
            !----- Zonal + Meridional wind speeds. ----------------------------------------!
            call mk_st_buff(mzp,mxp,myp,1,mtp,basic_g(ngrid)%up,scratch%scr3               &
                           ,i1-i0,i2-i0,j1-j0,j2-j0)
            call MPI_Pack(scratch%scr3,mtp,MPI_REAL,node_buffs_st(nm)%pack_send_buff       &
                         ,node_buffs_st(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
            call mk_st_buff(mzp,mxp,myp,1,mtp,basic_g(ngrid)%vp,scratch%scr3               &
                           ,i1-i0,i2-i0,j1-j0,j2-j0)
            call MPI_Pack(scratch%scr3,mtp,MPI_REAL,node_buffs_st(nm)%pack_send_buff       &
                         ,node_buffs_st(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
         case (6)
            !----- Vertical velocity and perturbation of the Exner function. --------------!
            call mk_st_buff(mzp,mxp,myp,1,mtp,basic_g(ngrid)%wp,scratch%scr3               &
                           ,i1-i0,i2-i0,j1-j0,j2-j0)
            call MPI_Pack(scratch%scr3,mtp,MPI_REAL,node_buffs_st(nm)%pack_send_buff       &
                         ,node_buffs_st(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)
            call mk_st_buff(mzp,mxp,myp,1,mtp,basic_g(ngrid)%pp,scratch%scr3               &
                           ,i1-i0,i2-i0,j1-j0,j2-j0)
            call MPI_Pack(scratch%scr3,mtp,MPI_REAL,node_buffs_st(nm)%pack_send_buff       &
                         ,node_buffs_st(nm)%nsend*f_ndmd_size,ipos,MPI_COMM_WORLD,ierr)

         end select
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Send the boundary condition.                                                 !
         !---------------------------------------------------------------------------------!
         mpiid = 2000000 + 10*maxgrds*(mchnum-1) + 10*(ngrid-1) + isflag
         call MPI_Isend(node_buffs_st(nm)%pack_send_buff,ipos-1,MPI_PACKED                 &
                       ,ipaths(5,itype,ngrid,nm),mpiid,MPI_COMM_WORLD,isend_req(nm),ierr)
         !---------------------------------------------------------------------------------!
      end if
   end do

   return
end subroutine node_sendst
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine receives the boundary conditions of the staggered variables from    !
! the nodes that are solving the neighbour grid points.                                    !
!------------------------------------------------------------------------------------------!
subroutine node_getst(isflag)
   use mem_grid
   use node_mod
   use mem_scratch
   use mem_basic
   implicit none

   !----- Module variables. ---------------------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Arguments -----------------------------------------------------------------------!
   integer                            , intent(in) :: isflag
   !----- Local variables -----------------------------------------------------------------!
   integer, dimension(MPI_STATUS_SIZE)             :: status
   integer                                         :: ierr
   integer                                         :: ipos
   integer                                         :: itype
   integer                                         :: nm
   integer                                         :: i1
   integer                                         :: i2
   integer                                         :: j1
   integer                                         :: j2
   integer                                         :: mtp
   integer                                         :: node_src
   integer                                         :: nptsxy
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Define the type of variable based on the input flag.                              !
   !---------------------------------------------------------------------------------------!
   select case (isflag)
   case (2:4)
      itype = isflag
   case (5:6)
      itype = 1
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     First, let's make sure our sends are all finished and de-allocated.               !
   !---------------------------------------------------------------------------------------!
   do nm=1,nmachs
      if (ipaths(1,itype,ngrid,nm) /= 0) call MPI_Wait(isend_req(nm),status,ierr)
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Now, let's wait on our receives.                                                  !
   !---------------------------------------------------------------------------------------!
   do nm=1,nmachs
      if (iget_paths(itype,ngrid,nm) /= 0) call MPI_Wait(irecv_req(nm),status,ierr)
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     We got all our stuff. Now unpack it into appropriate space.                       !
   !---------------------------------------------------------------------------------------!
   do nm=1,nmachs

      if (iget_paths(itype,ngrid,nm) /= 0) then
         ipos = 1

         !----- Obtain the size information from the buffer. ------------------------------!
         call MPI_Unpack(node_buffs_st(nm)%pack_recv_buff                                  &
                        ,node_buffs_st(nm)%nrecv*f_ndmd_size,ipos,i1,1,MPI_INTEGER         &
                        ,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_st(nm)%pack_recv_buff                                  &
                        ,node_buffs_st(nm)%nrecv*f_ndmd_size,ipos,i2,1,MPI_INTEGER         &
                        ,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_st(nm)%pack_recv_buff                                  &
                        ,node_buffs_st(nm)%nrecv*f_ndmd_size,ipos,j1,1,MPI_INTEGER         &
                        ,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_st(nm)%pack_recv_buff                                  &
                        ,node_buffs_st(nm)%nrecv*f_ndmd_size,ipos,j2,1,MPI_INTEGER         &
                        ,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_st(nm)%pack_recv_buff                                  &
                        ,node_buffs_st(nm)%nrecv*f_ndmd_size,ipos,node_src,1,MPI_INTEGER   &
                        ,MPI_COMM_WORLD,ierr)
         !---------------------------------------------------------------------------------!

         !----- Find the number of points in the horizontal and the number of points. -----!
         nptsxy=(i2-i1+1)*(j2-j1+1)
         mtp = nzp * nptsxy



         !---------------------------------------------------------------------------------!
         !     Decide which variable to extract from the package based on the input flag.  !
         !---------------------------------------------------------------------------------!
         select case (isflag)
         case (2)
            !----- Zonal wind. ------------------------------------------------------------!
            call MPI_Unpack(node_buffs_st(nm)%pack_recv_buff                               &
                           ,node_buffs_st(nm)%nrecv*f_ndmd_size,ipos,scratch%scr4,mtp      &
                           ,MPI_REAL,MPI_COMM_WORLD,ierr)
            call ex_st_buff(mzp,mxp,myp,1,mtp,basic_g(ngrid)%up,scratch%scr4               &
                           ,i1-i0,i2-i0,j1-j0,j2-j0)

         case (3)
            call MPI_Unpack(node_buffs_st(nm)%pack_recv_buff                               &
                           ,node_buffs_st(nm)%nrecv*f_ndmd_size,ipos,scratch%scr4,mtp      &
                           ,MPI_REAL,MPI_COMM_WORLD,ierr)
            call ex_st_buff(mzp,mxp,myp,1,mtp,basic_g(ngrid)%vp,scratch%scr4               &
                           ,i1-i0,i2-i0,j1-j0,j2-j0)

         case (4)
            call MPI_Unpack(node_buffs_st(nm)%pack_recv_buff                               &
                           ,node_buffs_st(nm)%nrecv*f_ndmd_size,ipos,scratch%scr4,mtp      &
                           ,MPI_REAL,MPI_COMM_WORLD,ierr)
            call ex_st_buff(mzp,mxp,myp,1,mtp,basic_g(ngrid)%pp,scratch%scr4               &
                           ,i1-i0,i2-i0,j1-j0,j2-j0)

         case (5)
            call MPI_Unpack(node_buffs_st(nm)%pack_recv_buff                               &
                           ,node_buffs_st(nm)%nrecv*f_ndmd_size,ipos,scratch%scr4,mtp      &
                           ,MPI_REAL,MPI_COMM_WORLD,ierr)
            call ex_st_buff(mzp,mxp,myp,1,mtp,basic_g(ngrid)%up,scratch%scr4               &
                           ,i1-i0,i2-i0,j1-j0,j2-j0)
            call MPI_Unpack(node_buffs_st(nm)%pack_recv_buff                               &
                           ,node_buffs_st(nm)%nrecv*f_ndmd_size,ipos,scratch%scr4,mtp      &
                           ,MPI_REAL,MPI_COMM_WORLD,ierr)
            call ex_st_buff(mzp,mxp,myp,1,mtp,basic_g(ngrid)%vp,scratch%scr4               &
                           ,i1-i0,i2-i0,j1-j0,j2-j0)

         case (6)
            call MPI_Unpack(node_buffs_st(nm)%pack_recv_buff                               &
                           ,node_buffs_st(nm)%nrecv*f_ndmd_size,ipos,scratch%scr4,mtp      &
                           ,MPI_REAL,MPI_COMM_WORLD,ierr)
            call ex_st_buff(mzp,mxp,myp,1,mtp,basic_g(ngrid)%wp,scratch%scr4               &
                           ,i1-i0,i2-i0,j1-j0,j2-j0)
            call MPI_Unpack(node_buffs_st(nm)%pack_recv_buff                               &
                           ,node_buffs_st(nm)%nrecv*f_ndmd_size,ipos,scratch%scr4,mtp      &
                           ,MPI_REAL,MPI_COMM_WORLD,ierr)
            call ex_st_buff(mzp,mxp,myp,1,mtp,basic_g(ngrid)%pp,scratch%scr4               &
                           ,i1-i0,i2-i0,j1-j0,j2-j0)

         end select
      end if
   end do

   return
end subroutine node_getst
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine simply copies a stretch of the sub-domain data into a buffer vector. !
!  X and Y (n2 and n3) are the only dimensions allowed to have sub-domains, and n1 and n4  !
! should always match the original dimension.                                              !
!------------------------------------------------------------------------------------------!
subroutine mk_st_buff(nz,nx,ny,ne,nlbc,arr,buff,iwest,ieast,jsouth,jnorth)
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
   real, dimension(nlbc)       , intent(inout) :: buff
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
      write (unit=*,fmt='(a)')       ' Mismatch betweeen buffer and ST domain sizes.'
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
      call abort_run('Mismatch in the buffer size!','mk_st_buff','mpass_st.f90')
   end if

   return
end subroutine mk_st_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine simply copies the buffer vector back to a 3D subdomain array.  X and !
! Y (n2 and n3) are the only dimensions allowed to have sub-domains, and n1 and n4 can be  !
! the third dimension, but they cannot be greater than 1 at the same time.                 !
!------------------------------------------------------------------------------------------!
subroutine ex_st_buff(nz,nx,ny,ne,nlbc,arr,buff,iwest,ieast,jsouth,jnorth)
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
      write (unit=*,fmt='(a)')       ' Mismatch betweeen buffer and ST domain sizes.'
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
      call abort_run('Mismatch in the buffer size!!!','ex_st_buff','mpass_st.f90')
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine ex_st_buff
!==========================================================================================!
!==========================================================================================!
