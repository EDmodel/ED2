!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


module obs_input

! Surface variables

integer, parameter :: max_sfc_vars=12

type ralph_sfc_obs
   integer :: iqflags(max_sfc_vars,3),ihgtflg
   character(len=16) :: id
   real :: lat,lon,elev,hgt,ff,dd,t,td,p
   integer :: jyear,jmonth,jdate,jtime
end type

type(ralph_sfc_obs) :: rsfc_obs

! Upper air variables

integer, parameter :: max_up_levs=6000,max_upa_vars=12

type ralph_upa_obs
   real, dimension(max_up_levs) :: p,t,z,r,zz,dz,fz
   character(len=16) :: id
   real :: lat,lon,elev
   integer :: lp,lz
   integer :: iqflagsp(max_up_levs,4,3),iqflagsz(max_up_levs,3,3)
   integer :: jyear,jmonth,jdate,jtime
end type

type(ralph_upa_obs) :: rupa_obs


! Header info

integer, parameter :: max_head_vars=10

type obs_header
   character(len=80) :: head_string(max_head_vars),sfc_string(max_sfc_vars)  &
                       ,therm_string(max_upa_vars),wind_string(max_upa_vars) &
                       ,sfc_units(max_sfc_vars) 
   integer :: iun,iver,nhead,nvsfc,nvtherm,nvwind
end type

type(obs_header) :: header(1)

end module
