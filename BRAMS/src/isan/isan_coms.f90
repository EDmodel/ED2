!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################


Module isan_coms
  use grid_dims, only : str_len

  !---------------------------------------------------------------------------
  !    Configuration COMMON blocks for RAMS isentropic data analysis package.
  !---------------------------------------------------------------------------
  !MAXPR         Maximum number of vertical levels that can be used in 
  !                 the pressure data.
  !MAXISN        Maximum number of vertical levels that can be used in 
  !                 the isentropic analysis.
  !MAXX          Maximum number of grid points (in the x or east-west direction) 
  !                 in any of the RAMS or pressure grids.
  !MAXY          Maximum number of grid points (in the y or north-south direction) 
  !                 in any of the RAMS or pressure grids.
  !MAXTIMES      Maximum number of data analysis times that can be processed 
  !                 in a single run.
  !MAXAGRDS      Maximum number of RAMS grids that can have varfiles generated.
  !MAXSIGZ       Maximum number of vertical levels that can be used 
  !                 in the  _z analysis.
  !MAXLEV	       Maximum number of levels in an input rawinsonde.
  !MAXSNAME      Maximum number of input observations
  !MAXISFILES    Maximum number of input data times
  !------------------------------------------------------------------------------------

  integer, parameter :: maxpr=100 ,maxisn=100 ,maxx=1000 ,maxy=1000  &
       ,maxtimes=5000 ,maxagrds=10    ,maxsigz=100  &
       ,maxlev=9999   ,maxsname=100000 ,maxisfiles=100000
  !---------------------------------------------------------------------------
  integer :: ioflgisz,ioflgvar,natime,iszstage,ivrstage,iyear,imonth,idate  &
       ,ihour,isan_inc,i1st_flg,iupa_flg,isfc_flg
  !---------------------------------------------------------------------------
  character(len=str_len)  :: innpr,inrawi,insrfce

  character(len=str_len) :: varpfx, iapr, iarawi, iasrfce 
  ! Modif. by ALF

  character(len=8)   :: pdata,guess1st

  !---------------------------------------------------------------------------
  !     Input pressure file header
  !---------------------------------------------------------------------------
  integer :: marker,isversion,iyy,imm,idd,ihh,itinc,inproj,ivertcoord
  real    :: xnelat,xnelon,cntlat,cntlon,secondlat

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Input pressure data memory

  real, allocatable, dimension(:,:,:) :: p_u,p_v,p_t,p_z,p_r,p_ur,p_vr
  real, allocatable, dimension(:,:)   :: p_lat,p_lon
  real, allocatable, dimension(:,:)   :: p_slp,p_sfp,p_sft,p_snow,p_sst
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Polar-stereo/pressure grid memory

  real, allocatable, dimension(:,:,:) :: pp_u,pp_v,pp_t,pp_z,pp_r
  real, allocatable, dimension(:,:)   :: pp_sglob
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Polar-stereo/isentropic grid memory

  real, allocatable, dimension(:,:,:) :: pi_u,pi_v,pi_p,pi_s,pi_r  &
       ,pi_scra,pi_scrb
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Polar-stereo/sigma-z grid memory
  !                         :: 
  real, allocatable, dimension(:,:,:) :: ps_u,ps_v,ps_p,ps_t,ps_r  &
       ,ps_scra,ps_scrb
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Polar-stereo/surface grid memory

  real, allocatable, dimension(:,:) :: rs_u,rs_v,rs_p,rs_t,rs_r,rs_s  &
       ,rs_top,rs_qual  &
       ,rs_slp,rs_sfp,rs_sft,rs_snow,rs_sst
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Data type to replace A array memory use in ISAN. 

  type isan_grids
     real, pointer, dimension(:,:,:) :: rr_u,rr_v,rr_t,rr_p,rr_r  
     real, pointer, dimension(:,:,:) :: rr_ug,rr_vg,rr_tg,rr_pg,rr_rg  
     real, pointer, dimension(:,:,:) :: rr_pi0,rr_th0,rr_dn0,rr_dn0u,rr_dn0v  
     real, pointer, dimension(:,:)   :: rr_slp,rr_sfp,rr_sft,rr_snow,rr_sst  
  end type isan_grids
  real, allocatable, dimension(:)    :: rr_scr1,rr_scr2,rr_vt2da

  type (isan_grids)                  :: is_grids(maxagrds)

  real, dimension(maxsigz,maxagrds)  :: piref,thref,dnref,rtref

  integer                            :: maxix,maxiy,maxiz

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Input observation data memory
  !
  real, allocatable, dimension(:,:)  :: up_uz,up_vz,up_ur,up_vr,up_zz  &
       ,up_p,up_t,up_z,up_r
  real, allocatable, dimension(:)    :: up_lat,up_lon,up_top
  real, allocatable, dimension(:,:)  :: up_topg
  integer, allocatable, dimension(:) :: up_lp, up_lz
  character(len=8), allocatable, dimension(:) :: up_chstid

  real, allocatable, dimension(:)    :: sf_u,sf_v,sf_p,sf_t,sf_s,sf_r
  real, allocatable, dimension(:)    :: sf_ur,sf_vr
  real, allocatable, dimension(:)    :: sf_lat,sf_lon,sf_top,sf_scra
  character(len=8), allocatable, dimension(:) :: sf_chstid
  character(len=14), allocatable, dimension(:) :: sf_date
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !  Upper air-isentropic/sigma-z memory
  !
  real, allocatable, dimension(:,:)  :: upi_u,upi_v,upi_p,upi_s,upi_r
  real, allocatable, dimension(:,:)  :: ups_u,ups_v,ups_p,ups_t,ups_r
  !---------------------------------------------------------------------------

  integer                          :: npdates
  integer, dimension(maxisfiles,4) :: iproc_flag
  !---------------------------------------------------------------------------
  character(len=128), dimension(maxisfiles,4) :: iproc_names
  !---------------------------------------------------------------------------
  character(len=14), dimension(maxisfiles) :: iproc_dates
  !---------------------------------------------------------------------------
  integer                   :: nprx,npry,nprz,idatelin,iglobew,iglobs,iglobn
  integer, dimension(maxpr) :: levpr
  real                      :: xswlon,xswlat,gdatdx,gdatdy
  real, dimension(maxpr)    :: pnpr
  !---------------------------------------------------------------------------
  integer                   :: nsta,nssfc,notsta  &
       ,maxsta,maxsfc,iobswin
  real                      :: stasep
  real, dimension(maxagrds) :: wvlnth,respon,swvlnth
  !---------------------------------------------------------------------------
  character(len=8), dimension(50)       :: notid
  !---------------------------------------------------------------------------
  integer                    :: nisx,nisy,nisn,interp,igridfl,nigrids,nsigz  &
       ,nfeedvar
  integer, dimension(maxisn) :: levth
  real                       :: gobsep,gobrad,topsigz,hybbot,hybtop,sfcinf  &
       ,sigzwt
  real, dimension(maxsigz)   :: sigz
  real, dimension(maxagrds)  :: gridwt
  !---------------------------------------------------------------------------

End module isan_coms
