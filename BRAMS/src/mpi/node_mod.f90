!==========================================================================================!
!                                                                                          !
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
!     This module contains many variables that are used for the communication among nodes. !
!------------------------------------------------------------------------------------------!
module node_mod

   use grid_dims, only : maxgrds   & ! intent(in)
                       , maxmach   ! ! intent(in)

   !---------------------------------------------------------------------------------------!
   !     The following variables hold the dimensions of the node sub-domain, and their     !
   ! absolute position within the full domain.  These are aliases and will be re-assigned  !
   ! every time a new grid is called.                                                      !
   !---------------------------------------------------------------------------------------!
   integer                                         :: mxp
   integer                                         :: myp
   integer                                         :: mzp
   integer                                         :: ia
   integer                                         :: iz
   integer                                         :: ja
   integer                                         :: jz
   integer                                         :: i0
   integer                                         :: j0
   integer                                         :: ibcon
   integer                                         :: ia_1
   integer                                         :: ia_2
   integer                                         :: ia_3
   integer                                         :: ia1
   integer                                         :: ia2
   integer                                         :: ia3
   integer                                         :: ja_1
   integer                                         :: ja_2
   integer                                         :: ja_3
   integer                                         :: ja1
   integer                                         :: ja2
   integer                                         :: ja3
   integer                                         :: iz_1
   integer                                         :: iz_2
   integer                                         :: iz_3
   integer                                         :: iz1
   integer                                         :: iz2
   integer                                         :: iz3
   integer                                         :: jz_1
   integer                                         :: jz_2
   integer                                         :: jz_3
   integer                                         :: jz1
   integer                                         :: jz2
   integer                                         :: jz3
   integer                                         :: izu
   integer                                         :: jzv
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     The following variables hold the machine number, and the "master" (input/output)  !
   ! node.                                                                                 !
   !---------------------------------------------------------------------------------------!
   integer                                         :: ipara
   integer                                         :: master_num
   integer                                         :: mchnum
   integer                                         :: mynum
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     The following variables contain the full information on all grids, and they will  !
   ! be used to assign the variables at the beginning of this module.                      !
   !---------------------------------------------------------------------------------------!
   integer, dimension(maxgrds)            , target :: mmxp
   integer, dimension(maxgrds)            , target :: mmyp
   integer, dimension(maxgrds)            , target :: mmzp
   integer, dimension(maxgrds)                     :: mia
   integer, dimension(maxgrds)                     :: miz
   integer, dimension(maxgrds)                     :: mja
   integer, dimension(maxgrds)                     :: mjz
   integer, dimension(maxgrds)                     :: mi0
   integer, dimension(maxgrds)                     :: mj0
   integer, dimension(maxgrds)                     :: mibcon
   integer, dimension(maxgrds)                     :: mnestflg
   integer, dimension(maxgrds)                     :: mfeednode
   !---------------------------------------------------------------------------------------!






   !---------------------------------------------------------------------------------------!
   !     The following variables contain the table with the location and size of each of   !
   ! the nodes sub-domain.                                                                 !
   !---------------------------------------------------------------------------------------!
   integer                                         :: nmachs
   integer                                         :: load_bal
   integer, dimension(maxmach)                     :: machs
   integer, dimension(maxmach,maxgrds)             :: nodemxp
   integer, dimension(maxmach,maxgrds)             :: nodemyp
   integer, dimension(maxmach,maxgrds)             :: nodemzp
   integer, dimension(maxmach,maxgrds)             :: nodeia
   integer, dimension(maxmach,maxgrds)             :: nodeiz
   integer, dimension(maxmach,maxgrds)             :: nodeja
   integer, dimension(maxmach,maxgrds)             :: nodejz
   integer, dimension(maxmach,maxgrds)             :: nodei0
   integer, dimension(maxmach,maxgrds)             :: nodej0
   integer, dimension(maxmach,maxgrds)             :: nodeibcon
   integer, dimension(maxmach,maxgrds)             :: nodenestflg
   integer, dimension(maxmach,maxgrds)             :: nodefeednode
   integer, dimension(maxmach,maxgrds,4)           :: nodeconn
   integer, dimension(maxgrds,8)                   :: nodebounds ! Reprod.-Saulo Barros
   integer, dimension(5,7,maxgrds,maxmach)         :: ipaths
   integer, dimension(6,maxgrds,maxmach)           :: iget_paths
   !---------------------------------------------------------------------------------------!





   !---------------------------------------------------------------------------------------!
   !     Other variables, not sure for what they are used.                                 !
   !---------------------------------------------------------------------------------------!
   integer                                         :: newbuff_feed
   integer                                         :: nbuff_feed
   integer                                         :: newbuff_nest
   integer                                         :: nbuff_nest
   integer                                         :: f_ndmd_size
   integer                                         :: nbuff_st
   integer                                         :: nbuff_adv
   !---------------------------------------------------------------------------------------!
   integer, dimension(maxmach)                     :: irecv_req
   integer, dimension(maxmach)                     :: isend_req
   integer, dimension(maxmach)                     :: nsend_buff
   integer, dimension(maxmach)                     :: nrecv_buff
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     This structure contains the buffers that are used for sending/receiving the       !
   ! lateral boundary condition.                                                           !
   !---------------------------------------------------------------------------------------!
   type pack_buffs
      character, dimension(:), pointer :: pack_send_buff
      character, dimension(:), pointer :: pack_recv_buff
      integer                          :: nsend
      integer                          :: nrecv
   end type pack_buffs
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     The buffers.                                                                      !
   !---------------------------------------------------------------------------------------!
   type(pack_buffs), dimension(  maxmach) :: node_buffs_lbc
   type(pack_buffs), dimension(  maxmach) :: node_buffs_feed
   type(pack_buffs), dimension(  maxmach) :: node_buffs_nest
   type(pack_buffs), dimension(6,maxmach) :: node_buffs_st
   type(pack_buffs), dimension(3,maxmach) :: node_buffs_adv
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     The default size of integer and real numbers, used to allocate the buffers.       !
   !---------------------------------------------------------------------------------------!
   integer                                         :: mpi_int_size
   integer                                         :: mpi_real_size
   !---------------------------------------------------------------------------------------!



   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine allocates the buffers.                                            !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_node_buff(this_buff,nbuff,number_size)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(pack_buffs), intent(inout) :: this_buff
      integer         , intent(in)    :: nbuff
      integer         , intent(in)    :: number_size
      !------------------------------------------------------------------------------------!

      call nullify_node_buff(this_buff)

      !------------------------------------------------------------------------------------!
      !      Set the buffer size to the maximum needed size plus 100, which is to account  !
      ! for some integer dimensions that are packed at the beginning of the message.       !
      !------------------------------------------------------------------------------------!
      this_buff%nsend = (nbuff + 100) * number_size
      this_buff%nrecv = (nbuff + 100) * number_size
      
      allocate(this_buff%pack_send_buff(this_buff%nsend))
      allocate(this_buff%pack_recv_buff(this_buff%nrecv))

      return
   end subroutine alloc_node_buff
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine nullifies the buffers, which is important as a first step to make !
   ! sure that the buffers are cleanly allocated.                                          !
   !---------------------------------------------------------------------------------------!
   subroutine nullify_node_buff(this_buff)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(pack_buffs), intent(inout) :: this_buff
      !------------------------------------------------------------------------------------!
      
      nullify(this_buff%pack_send_buff)
      nullify(this_buff%pack_recv_buff)

      return
   end subroutine nullify_node_buff
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine nullifies the buffers, which is important as a first step to make !
   ! sure that the buffers are cleanly allocated.                                          !
   !---------------------------------------------------------------------------------------!
   subroutine dealloc_node_buff(this_buff)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(pack_buffs), intent(inout) :: this_buff
      !------------------------------------------------------------------------------------!
      
      this_buff%nsend = 0
      this_buff%nrecv = 0
      if (associated(this_buff%pack_send_buff)) deallocate(this_buff%pack_send_buff)
      if (associated(this_buff%pack_send_buff)) deallocate(this_buff%pack_recv_buff)

      return
   end subroutine dealloc_node_buff
   !=======================================================================================!
   !=======================================================================================!

end module node_mod
!==========================================================================================!
!==========================================================================================!
