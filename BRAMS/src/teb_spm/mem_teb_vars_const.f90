!############################# Change Log ##################################
! 5.0.2
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Brazilian Regional Atmospheric Modeling System - BRAMS
!###########################################################################

Module teb_vars_const
use grid_dims
!---------------------------------------------------------------------------
 integer         ::iteb !Flag for urban parameterization using TEB  - EDF
                        !This flag is set up in RAMSIN 1=on - 0=off

 integer         ::nteb !number of roof, road and wall layers used in TEB  - EDF

real, dimension(MAXSTEB) :: D_ROAD,TC_ROAD,HC_ROAD	    &
		     ,D_WALL,TC_WALL,HC_WALL	    &
		     ,D_ROOF,TC_ROOF,HC_ROOF

real             :: TMINBLD   !Minimum internal building temperature
real             :: RUSHH1    !Morning Rush Hour
real             :: RUSHH2    !Afternoon/Evening Rush Hour
real             :: DAYLIGHT  !Daylight saving time (horario de verao)

integer          :: NURBTYPE    !Number of urban types

integer, dimension (maxubtp) :: ileafcod   !Leaf class code to identify 
                                        !each urban type
real, dimension(maxubtp) :: Z0_TOWN,BLD,BLD_HEIGHT         &
			    ,BLD_HL_RATIO,AROOF,EROOF      &
			    ,AROAD,EROAD,AWALL,EWALL,HTRAF &
			    ,HINDU,PLETRAF,PLEINDU     

real, parameter :: XPI= 3.1415                 &
                  ,XKARMAN= 0.4                &
                  ,XSTEFAN = 5.6697E-08        &
                  ,XDAY= 86400.                &
                  ,XG=9.80665                  &
                  ,XRD=287.0597                &
                  ,XRV=461.5250                &
                  ,XCPD=1004.70895             &
                  ,XLVTT=2.5008E+6             &
                  ,XTT=273.16                  !&

end Module
