!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module ref_sounding

   use grid_dims

   integer :: iref,jref
   real :: topref
   real, dimension(nzpmax,maxgrds) :: u01dn,v01dn,pi01dn,th01dn,dn01dn,rt01dn,co201dn
   
   !-------------------------------------------------------------------------------
   
   integer, parameter         :: maxsndg=1296
   integer                    :: ipsflg,itsflg,irtsflg,iusflg,nsndg
   real, dimension(maxsndg)   :: us,vs,ts,thds,ps,hs,rts,co2s

end Module ref_sounding
