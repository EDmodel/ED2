!
! Copyright (C) 1991-2004  ; All Rights Reserved ; Colorado State University
! Colorado State University Research Foundation ; ATMET, LLC
! 
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this 
! program; if not, write to the Free Software Foundation, Inc., 
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================
!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################


Module ed_node_coms

  use ed_max_dims, only : maxgrds,maxmach

  !---------------------------------------------------------------------------
  integer :: mxp,myp,mzp,ia,iz,ja,jz,i0,j0,master_num,mchnum,ibcon,ipara
  integer :: mynum,sendnum,recvnum
  !---------------------------------------------------------------------------
  integer, target, dimension(maxgrds) :: mmxp,mmyp,mmzp
  integer, dimension(maxgrds) :: mia,miz,mja,mjz  &
                                ,mi0,mj0,mibcon,mnestflg,mfeednode
  integer, dimension(maxgrds) :: iwest,ieast,iskip,jsouth,jnorth,jskip
  !---------------------------------------------------------------------------
  integer                                 :: nmachs,nnodetot
  integer, dimension(maxmach)             :: machs
  integer, dimension(maxmach,maxgrds)     :: nodemxp,nodemyp,nodemzp  &
                                            ,nodeia,nodeiz,nodeja,nodejz  &
                                            ,nodei0,nodej0,nodeibcon
  !---------------------------------------------------------------------------

end module ed_node_coms
