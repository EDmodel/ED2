!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module rpara

  use grid_dims

  !---------------------------------------------------------------------------
  integer                                         :: mainnum,nmachs,iparallel &
       ,load_bal
  integer, dimension(maxmach)                     :: machnum,nbuff_nest1  &
       ,newbuff_nest1
  integer, dimension(maxmach,maxgrds)             :: nxbeg,nxend,nybeg  &
       ,nyend,nxbegc,nxendc  &
       ,nybegc,nyendc,ixoff  &
       ,iyoff,npxy,ibcflg  &
       ,nestflg,ifeednode  &
       ,ixb,ixe,iyb,iye

  integer, dimension(maxgrds,8,maxmach)           :: ibounds ! Saulo Barros

  integer, dimension(maxmach,maxgrds,4)           :: nextnode
  integer, dimension(5,7,maxgrds,maxmach,maxmach) :: inode_paths_master
  integer, dimension(6,maxgrds,maxmach,maxmach)   :: iget_paths_master
  integer, dimension(2,maxmach,maxmach)           :: lbc_buffs
  real, dimension(maxmach)                        :: hperf
  real, dimension(maxmach,maxgrds)                :: perf
  real, dimension(maxmach,2)                      :: ptimes
  !---------------------------------------------------------------------------

end module rpara
