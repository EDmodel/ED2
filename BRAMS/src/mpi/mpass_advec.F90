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
subroutine node_sendadv(iaflag)

   use mem_grid
   use node_mod

   use var_tables
   use mem_scratch
   use mem_cuparm , only : nclouds  ! ! intent(in)
   use grid_dims  , only : maxgrds  ! ! intent(in)
   use mem_aerad  , only : nwave    ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)   :: iaflag
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
   logical               :: sendnow
   !----- Module variables. ---------------------------------------------------------------!
#if defined(RAMS_MPI)
   include 'interface.h'
   include 'mpif.h'
#endif
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Select case depending on the type of variable to be exchanged.                    !
   !---------------------------------------------------------------------------------------!
   select case (iaflag)
   case (2)
      !----- Zonal wind variables. --------------------------------------------------------!
      itype = 2

   case (3)
      !----- Meridional wind variables. ---------------------------------------------------!
      itype = 3

   case default
      !----- Thermodynamic variables or vertical wind. ------------------------------------!
      itype = 1
   end select
   !---------------------------------------------------------------------------------------!

#if defined(RAMS_MPI)

   !---------------------------------------------------------------------------------------!
   !     First, before we send anything, let's post the receives.  Also, make sure any     !
   ! pending sends are complete.                                                           !
   !---------------------------------------------------------------------------------------!
   do nm=1,nmachs
      if (iget_paths(itype,ngrid,nm) /= 0) then
         !----- Find the unique MPI flag to make sure the right message is got. -----------!
         mpiid= 6000000 + 10*maxgrds*(machs(nm)-1) + 10*(ngrid-1) + iaflag

         !----- Post the receive. ---------------------------------------------------------!
         call MPI_Irecv(node_buffs_adv(iaflag,nm)%pack_recv_buff                           &
                       ,node_buffs_adv(iaflag,nm)%nrecv                                    &
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
         call MPI_Pack(i1,1,MPI_INTEGER,node_buffs_adv(iaflag,nm)%pack_send_buff           &
                      ,node_buffs_adv(iaflag,nm)%nsend,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(i2,1,MPI_INTEGER,node_buffs_adv(iaflag,nm)%pack_send_buff           &
                      ,node_buffs_adv(iaflag,nm)%nsend,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(j1,1,MPI_INTEGER,node_buffs_adv(iaflag,nm)%pack_send_buff           &
                      ,node_buffs_adv(iaflag,nm)%nsend,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(j2,1,MPI_INTEGER,node_buffs_adv(iaflag,nm)%pack_send_buff           &
                      ,node_buffs_adv(iaflag,nm)%nsend,ipos,MPI_COMM_WORLD,ierr)
         call MPI_Pack(mynum,1,MPI_INTEGER,node_buffs_adv(iaflag,nm)%pack_send_buff        &
                      ,node_buffs_adv(iaflag,nm)%nsend,ipos,MPI_COMM_WORLD,ierr)



         !----- Send boundaries from the sought advection variables. ----------------------!

         do nv = 1,num_var(ngrid)

            !------------------------------------------------------------------------------!
            !     Decide whether to send the variable or not depending on the value of     !
            ! IAFLAG.                                                                      !
            !------------------------------------------------------------------------------!
            select case (iaflag)
            case (1)
               sendnow = vtab_r(nv,ngrid)%iadvt == 1
            case (2)
               sendnow = vtab_r(nv,ngrid)%iadvu == 1
            case (3)
               sendnow = vtab_r(nv,ngrid)%iadvv == 1
            case (4)
               sendnow = vtab_r(nv,ngrid)%iadvw == 1
            case (5)
               sendnow = trim(vtab_r(nv,ngrid)%name) == 'SCAL_IN'
            end select
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Add the sought variables to the buffer.                                  !
            !------------------------------------------------------------------------------!
            if (sendnow) then

               !----- Find the size of the vertical and environmental axes. ---------------!
               call ze_dims(ngrid,vtab_r(nv,ngrid)%idim_type,.false.,fdzp,fdep)
               !---------------------------------------------------------------------------!


               !----- Find the size of the lateral boundary condition. --------------------!
               mtp = nptsxy * fdzp * fdep

               !----- Copy the variable to the scratch array. -----------------------------!
               call azero(mtp,scratch%scr3)
               call mk_adv_buff(fdzp,mxp,myp,fdep,mtp,vtab_r(nv,ngrid)%var_p               &
                               ,scratch%scr3,i1-i0,i2-i0,j1-j0,j2-j0)
               !---------------------------------------------------------------------------!

               !----- Add the boundary condition to the MPI buffer. -----------------------!
               call MPI_Pack(scratch%scr3 ,mtp,MPI_REAL                                    &
                            ,node_buffs_adv(iaflag,nm)%pack_send_buff                      &
                            ,node_buffs_adv(iaflag,nm)%nsend,ipos,MPI_COMM_WORLD,ierr)
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!



         !----- Send out the stuff to the node. -------------------------------------------!
         mpiid = 6000000 + 10*maxgrds*(mchnum-1) + 10*(ngrid-1) + iaflag
         call MPI_Isend(node_buffs_adv(iaflag,nm)%pack_send_buff,ipos,MPI_PACKED           &
                       ,ipaths(5,itype,ngrid,nm),mpiid,MPI_COMM_WORLD,isend_req(nm),ierr)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!
#endif
   return
end subroutine node_sendadv
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine gets the lateral boundary condition information and copy to the     !
! edges of the node sub-domain.                                                            !
!------------------------------------------------------------------------------------------!
subroutine node_getadv(iaflag)
   use mem_grid
   use node_mod

   use var_tables
   use mem_scratch
   use mem_cuparm , only : nclouds  ! ! intent(in)
   use mem_aerad  , only : nwave    ! ! intent(in)

   implicit none
   !----- Module variables. ---------------------------------------------------------------!
#if defined(RAMS_MPI)
   include 'interface.h'
   include 'mpif.h'
#endif
   !----- Arguments. ----------------------------------------------------------------------!
   integer                            , intent(in)   :: iaflag
   !----- Local variables. ----------------------------------------------------------------!
#if defined(RAMS_MPI)
   integer, dimension(MPI_STATUS_SIZE)               :: status1
#endif
   integer                                           :: ierr
   integer                                           :: ipos
   integer                                           :: itype
   integer                                           :: nm
   integer                                           :: ibytes
   integer                                           :: msgid
   integer                                           :: ihostnum
   integer                                           :: i1
   integer                                           :: i2
   integer                                           :: j1
   integer                                           :: j2
   integer                                           :: fdzp
   integer                                           :: fdep
   integer                                           :: nmp
   integer                                           :: nv
   integer                                           :: node_src
   integer                                           :: mtc
   integer                                           :: mtp
   integer                                           :: nptsxy
   logical                                           :: getnow
   !---------------------------------------------------------------------------------------!


#if defined(RAMS_MPI)
   !---------------------------------------------------------------------------------------!
   !     Select case depending on the type of variable to be exchanged.                    !
   !---------------------------------------------------------------------------------------!
   select case (iaflag)
   case (2)
      !----- Zonal wind variables. --------------------------------------------------------!
      itype = 2

   case (3)
      !----- Meridional wind variables. ---------------------------------------------------!
      itype = 3

   case default
      !----- Thermodynamic variables or vertical wind. ------------------------------------!
      itype = 1
   end select
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
         call MPI_Unpack(node_buffs_adv(iaflag,nm)%pack_recv_buff                          &
                        ,node_buffs_adv(iaflag,nm)%nrecv                                   &
                        ,ipos,i1,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_adv(iaflag,nm)%pack_recv_buff                          &
                        ,node_buffs_adv(iaflag,nm)%nrecv                                   &
                        ,ipos,i2,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_adv(iaflag,nm)%pack_recv_buff                          &
                        ,node_buffs_adv(iaflag,nm)%nrecv                                   &
                        ,ipos,j1,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_adv(iaflag,nm)%pack_recv_buff                          &
                        ,node_buffs_adv(iaflag,nm)%nrecv                                   &
                        ,ipos,j2,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         call MPI_Unpack(node_buffs_adv(iaflag,nm)%pack_recv_buff                          &
                        ,node_buffs_adv(iaflag,nm)%nrecv                                   &
                        ,ipos,node_src,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

         !----- Total number of horizontal points. ----------------------------------------!
         nptsxy=(i2-i1+1)*(j2-j1+1)

         varloop: do nv = 1,num_var(ngrid)
            !------------------------------------------------------------------------------!
            !     Decide whether to send the variable or not depending on the value of     !
            ! IAFLAG.                                                                      !
            !------------------------------------------------------------------------------!
            select case (iaflag)
            case (1)
               getnow = vtab_r(nv,ngrid)%iadvt == 1
            case (2)
               getnow = vtab_r(nv,ngrid)%iadvu == 1
            case (3)
               getnow = vtab_r(nv,ngrid)%iadvv == 1
            case (4)
               getnow = vtab_r(nv,ngrid)%iadvw == 1
            case (5)
               getnow = trim(vtab_r(nv,ngrid)%name) == 'SCAL_IN'
            end select
            !------------------------------------------------------------------------------!


            if (getnow) then

               !----- Find the size of the vertical and environmental axes. ---------------!
               call ze_dims(ngrid,vtab_r(nv,ngrid)%idim_type,.false.,fdzp,fdep)
               !---------------------------------------------------------------------------!


               !----- Find the total number of points. ------------------------------------!
               mtp= fdzp * nptsxy * fdep
               !---------------------------------------------------------------------------!


               !----- Unpack the buffer into the scratch array. ---------------------------!
               call MPI_Unpack(node_buffs_adv(iaflag,nm)%pack_recv_buff                    &
                              ,node_buffs_adv(iaflag,nm)%nrecv,ipos,scratch%scr4,mtp       &
                              ,MPI_REAL,MPI_COMM_WORLD,ierr)
               !---------------------------------------------------------------------------!


               !----- Extract the information. --------------------------------------------!
               call ex_adv_buff(fdzp,mxp,myp,fdep,mtp,vtab_r(nv,ngrid)%var_p               &
                               ,scratch%scr4,i1-i0,i2-i0,j1-j0,j2-j0)
               !---------------------------------------------------------------------------!
            end if
         end do varloop
      end if
   end do machloop
#endif
   return
end subroutine node_getadv
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine simply copies a stretch of the sub-domain data into a buffer vector. !
!  X and Y (n2 and n3) are the only dimensions allowed to have sub-domains, and n1 and n4  !
! should always match the original dimension.                                              !
!------------------------------------------------------------------------------------------!
subroutine mk_adv_buff(nz,nx,ny,ne,nadv,arr,buff,iwest,ieast,jsouth,jnorth)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: ne
   integer                     , intent(in)    :: nadv
   integer                     , intent(in)    :: iwest
   integer                     , intent(in)    :: ieast
   integer                     , intent(in)    :: jsouth
   integer                     , intent(in)    :: jnorth
   real, dimension(nz,nx,ny,ne), intent(in)    :: arr
   real, dimension(nadv)       , intent(out)   :: buff
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
   if (ind /= nadv) then
      write (unit=*,fmt='(a)')       '------------------------------------------------'
      write (unit=*,fmt='(a)')       ' Mismatch betweeen buffer and ADV domain sizes.'
      write (unit=*,fmt='(a)')       ' '
      write (unit=*,fmt='(a,1x,i6)') ' - NADV   = ',nadv
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
      call abort_run('Mismatch in the buffer size!','mk_adv_buff','mpass_adv.f90')
   end if

   return
end subroutine mk_adv_buff
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine simply copies the buffer vector back to a 3D subdomain array.  X and !
! Y (n2 and n3) are the only dimensions allowed to have sub-domains, and n1 and n4         !
! should always match the original dimension.                                              !
!------------------------------------------------------------------------------------------!
subroutine ex_adv_buff(nz,nx,ny,ne,nadv,arr,buff,iwest,ieast,jsouth,jnorth)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: ne
   integer                     , intent(in)    :: nadv
   integer                     , intent(in)    :: iwest
   integer                     , intent(in)    :: ieast
   integer                     , intent(in)    :: jsouth
   integer                     , intent(in)    :: jnorth
   real, dimension(nz,nx,ny,ne), intent(inout) :: arr
   real, dimension(nadv)       , intent(in)    :: buff
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
   if (ind /= nadv) then
      write (unit=*,fmt='(a)')       '------------------------------------------------'
      write (unit=*,fmt='(a)')       ' Mismatch betweeen buffer and ADV domain sizes.'
      write (unit=*,fmt='(a)')       ' '
      write (unit=*,fmt='(a,1x,i6)') ' - NADV   = ',nadv
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
      call abort_run('Mismatch in the buffer size!!!','ex_adv_buff','mpass_adv.f90')
   end if

   return
end subroutine ex_adv_buff
!==========================================================================================!
!==========================================================================================!
