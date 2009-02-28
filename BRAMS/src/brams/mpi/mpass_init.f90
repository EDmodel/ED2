!===================================== Change Log =========================================!
! 5.0.0                                                                                    !
!                                                                                          !
! MLO - 09/29/08 Including the new Grell related variables.                                !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This routine gives basic processor ID info to the nodes.                              !
!------------------------------------------------------------------------------------------!
subroutine masterput_processid(nproc,taskids,master_num)

   use rpara

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Arguments -----------------------------------------------------------------------!
   integer, dimension(*), intent(in) :: taskids
   integer              , intent(in) :: master_num,nproc
   !----- Local variables -----------------------------------------------------------------!
   integer                           :: nm
   integer                           :: ierr
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Copying the arguments to rpara structure variables, which will be available to the !
   ! head node.                                                                            !
   !---------------------------------------------------------------------------------------!
   mainnum=master_num
   nmachs=nproc
   do nm=1,nmachs
      machnum(nm)=taskids(nm)
   end do

   !----- Sending the ID information to the nodes -----------------------------------------!
   do nm=1,nmachs
     call MPI_Send(mainnum,1,MPI_INTEGER,machnum(nm),11,MPI_COMM_WORLD,ierr)
     call MPI_Send(machnum(nm),1,MPI_INTEGER,machnum(nm),12,MPI_COMM_WORLD,ierr)
     call MPI_Send(nm,1,MPI_INTEGER,machnum(nm),13,MPI_COMM_WORLD,ierr)
     call MPI_Send(nmachs,1,MPI_INTEGER,machnum(nm),14,MPI_COMM_WORLD,ierr)
     call MPI_Send(machnum,nmachs,MPI_INTEGER,machnum(nm),15,MPI_COMM_WORLD,ierr)
   end do
   return
end subroutine masterput_processid
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will transfer the variables read at the namelist (and some close      !
! friends) to the modules at the slave nodes).                                             !
!------------------------------------------------------------------------------------------!
subroutine masterput_nl(master_num)

   use mem_all
   use rpara
   use therm_lib          , only : level                       & ! intent(in)
                                  ,vapour_on                   & ! intent(in)
                                  ,cloud_on                    & ! intent(in)
                                  ,bulk_on                     ! ! intent(in)
   use mem_mass           , only : iexev                       & ! intent(in)
                                  ,imassflx                    ! ! intent(in)
   use grell_coms         , only:  closure_type                & ! intent(in)
                                  ,maxclouds                   & ! intent(in)
                                  ,iupmethod                   & ! intent(in)
                                  ,depth_min                   & ! intent(in)
                                  ,cap_maxs                    & ! intent(in)
                                  ,maxens_lsf                  & ! intent(in)
                                  ,maxens_dyn                  & ! intent(in)
                                  ,maxens_eff                  & ! intent(in)
                                  ,maxens_cap                  & ! intent(in)
                                  ,iupmethod                   & ! intent(in)
                                  ,iupstrm                     & ! intent(in)
                                  ,radius                      & ! intent(in)
                                  ,zkbmax                      & ! intent(in)
                                  ,max_heat                    & ! intent(in)
                                  ,zcutdown                    & ! intent(in)
                                  ,z_detr                      ! ! intent(in)
   use sib_vars           , only : n_co2                       & ! intent(in)
                                  ,co2_init                    ! ! intent(in)
   use catt_start         , only : catt                        ! ! intent(in)
   use mem_globrad        , only : raddatfn                    ! ! intent(in)
   use emission_source_map, only : plumerise                   ! ! intent(in)
   use plume_utils        , only : prfrq                       ! ! intent(in)
   use teb_spm_start      , only : teb_spm                     ! ! intent(in)
   use teb_vars_const     , only : iteb,tminbld,nteb           & ! intent(in)
                                  ,rushh1,rushh2,daylight      & ! intent(in)
                                  ,d_road,tc_road,hc_road      & ! intent(in)
                                  ,d_wall,tc_wall,hc_wall      & ! intent(in)
                                  ,d_roof,tc_roof,hc_roof      & ! intent(in)
                                  ,nurbtype,ileafcod,z0_town   & ! intent(in)
                                  ,bld,bld_height,bld_hl_ratio & ! intent(in)
                                  ,aroof,eroof,aroad,eroad     & ! intent(in)
                                  ,awall,ewall,htraf,hindu     & ! intent(in)
                                  ,pletraf,pleindu             ! ! intent(in)
   use mem_emiss          , only : ichemi                      & ! intent(in)
                                  ,ichemi_in                   & ! intent(in)
                                  ,chemdata_in                 & ! intent(in)
                                  ,isource                     & ! intent(in)
                                  ,weekdayin                   & ! intent(in)
                                  ,efsat                       & ! intent(in)
                                  ,efsun                       & ! intent(in)
                                  ,eindno                      & ! intent(in)
                                  ,eindno2                     & ! intent(in)
                                  ,eindpm                      & ! intent(in)
                                  ,eindco                      & ! intent(in)
                                  ,eindso2                     & ! intent(in)
                                  ,eindvoc                     & ! intent(in)
                                  ,eveino                      & ! intent(in)
                                  ,eveino2                     & ! intent(in)
                                  ,eveipm                      & ! intent(in)
                                  ,eveico                      & ! intent(in)
                                  ,eveiso2                     & ! intent(in)
                                  ,eveivoc                       ! intent(in)
   use ref_sounding, only:         maxsndg                     ! ! intent(in)

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: master_num
   !----- Local variables -----------------------------------------------------------------!
   integer :: nm,ierr
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Now all the intended variables are broadcast, they will be caught by the nodes.    !
   ! Here THE ORDER MATTERS A LOT. If you are including new variables, make sure that the  !
   ! variables listed here are in the very same order as the nodeget_nl subroutine (later  !
   ! in this file), otherwise horrible things will happen to your run.                     !
   !---------------------------------------------------------------------------------------!

   call MPI_Bcast(TIMMAX,1,MPI_DOUBLE_PRECISION,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(if_adap,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(load_bal,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NGRIDS,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNXP,MAXGRDS,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNYP,MAXGRDS,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNZP,MAXGRDS,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NZG,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NZS,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NXTNEST,MAXGRDS,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IHTRAN,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DELTAX,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DELTAY,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DELTAZ,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DZRAT,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DZMAX,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ZZ,NZPMAX,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)


   call MPI_Bcast(IDELTAT,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NACOUST,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IMONTHA,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IDATEA,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IYEARA,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ITIMEA,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IMONTHZ,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IDATEZ,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IYEARZ,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ITIMEZ,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(CENTLAT,MAXGRDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(CENTLON,MAXGRDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(POLELAT,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(POLELON,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(PLATN,MAXGRDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(PLONN,MAXGRDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(NSTRATX,MAXGRDS,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NSTRATY,MAXGRDS,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NSTRATZ1,NZPMAX,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NSTRATZ2,NZPMAX,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NESTZ1,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NESTZ2,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NINEST,MAXGRDS,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NJNEST,MAXGRDS,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NKNEST,MAXGRDS,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(GRIDU,MAXGRDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(GRIDV,MAXGRDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNSTTOP,MAXGRDS,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNSTBOT,MAXGRDS,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)


   call MPI_Bcast(IOUTPUT,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(INITFLD,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(FRQPRT,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(FRQHIS,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(FRQANL,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(FRQLITE,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NLITE_VARS,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(XLITE(1:20),20,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(YLITE(1:20),20,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ZLITE(1:20),20,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   do nm = 1, nlite_vars
      call MPI_Bcast(LITE_VARS(nm),32,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
   enddo


   call MPI_Bcast(FRQMEAN,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(FRQBOTH,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(AVGTIM,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(INITIAL,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NUD_TYPE,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IF_ODA,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NUDLAT,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(TNUDLAT,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(TNUDTOP,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(TNUDCENT,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ZNUDTOP,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(WT_NUDGE_GRID,maxgrds,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(WT_NUDGE_UV,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(WT_NUDGE_TH,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(WT_NUDGE_PI,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(WT_NUDGE_RT,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NUD_COND,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(TCOND_BEG,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(TCOND_END,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(T_NUDGE_RC,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(WT_NUDGEC_GRID,maxgrds,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(IUPDSST,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ISSTFLG,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   
   call MPI_Bcast(IUPDNDVI,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NDVIFLG,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   
   call MPI_Bcast(RUNTYPE(1:16),16,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(DTLONG,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNDTRAT,maxgrds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(NADDSC,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(NNQPARM,MAXGRDS,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NCLOUDS,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NDEEPEST,MAXGRDS,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NSHALLOWEST,MAXGRDS,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(WCLDBS,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(CONFRQ,MAXCLOUDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(CPTIME,MAXCLOUDS,MPI_DOUBLE_PRECISION,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IUPMETHOD,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IUPSTRM,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(RADIUS,MAXCLOUDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DEPTH_MIN,MAXCLOUDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr) 
   call MPI_Bcast(CAP_MAXS,MAXCLOUDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr) 
   call MPI_Bcast(ZKBMAX,MAXCLOUDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)   
   call MPI_Bcast(ZCUTDOWN,MAXCLOUDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr) 
   call MPI_Bcast(Z_DETR,MAXCLOUDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)   
   call MPI_Bcast(MAX_HEAT,MAXCLOUDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr) 
   do nm=1,maxclouds
      call MPI_Bcast(CLOSURE_TYPE(nm),2,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr) 
   end do
   call MPI_Bcast(MAXENS_LSF,maxclouds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(MAXENS_EFF,maxclouds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(MAXENS_CAP,maxclouds,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)


   call MPI_Bcast(SLZ,NZGMAX,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(STGOFF,NZGMAX,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(SLMSTR,NZGMAX,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IDIFFK,MAXGRDS,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IBRUVAIS,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IHORGRAD,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IF_URBAN_CANOPY,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   
   call MPI_Bcast(CATT,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(RADDATFN(1:160),160,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)! 


   if (CATT == 1) then
     call MPI_Bcast(PLUMERISE,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(PRFRQ,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   end if

   call MPI_Bcast(TEB_SPM,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   ! TEB_SPM
   if (TEB_SPM == 1) then
     call MPI_Bcast(RUSHH1,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(RUSHH2,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(DAYLIGHT,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(ITEB,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(TMINBLD,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(NTEB,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(TC_ROOF,MAXSTEB,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(D_ROOF,MAXSTEB,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(HC_ROOF,MAXSTEB,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(HC_ROAD,MAXSTEB,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(TC_ROAD,MAXSTEB,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(D_ROAD,MAXSTEB,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(HC_WALL,MAXSTEB,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(TC_WALL,MAXSTEB,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(D_WALL,MAXSTEB,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(NURBTYPE,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(ILEAFCOD,MAXUBTP,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(Z0_TOWN,MAXUBTP,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(BLD,MAXUBTP,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(BLD_HEIGHT,MAXUBTP,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(BLD_HL_RATIO,MAXUBTP,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(AROOF,MAXUBTP,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EROOF,MAXUBTP,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(AROAD,MAXUBTP,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EROAD,MAXUBTP,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(AWALL,MAXUBTP,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EWALL,MAXUBTP,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(AROAD,MAXUBTP,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(HTRAF,MAXUBTP,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(HINDU,MAXUBTP,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(PLETRAF,MAXUBTP,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(PLEINDU,MAXUBTP,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

     call MPI_Bcast(ICHEMI,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(ICHEMI_IN,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(CHEMDATA_IN(1:80),80,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(ISOURCE,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(WEEKDAYIN(1:3),3,MPI_CHARACTER,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EFSAT,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EFSUN,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EINDNO,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EINDNO2,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EINDPM,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EINDCO,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EINDSO2,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EINDVOC,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EVEINO,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EVEINO2,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EVEIPM,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EVEICO,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EVEISO2,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EVEIVOC,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   end if
   
   call MPI_Bcast(LSFLG,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IBND,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(JBND,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NPATCH,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NVEGPAT,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ISFCL,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(N_CO2,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(CO2_INIT,maxsndg,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   
   call MPI_Bcast(LONRAD,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(CPHAS,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DISTIM,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NFPT,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DTHCON,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DRTCON,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(SEATMP,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(PCTLCON,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NSLCON,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NVGCON,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(RADFRQ,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ZROUGH,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ILWRTYP,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ISWRTYP,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ICUMFDBK,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ICORFLG,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IEXEV,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IMASSFLX,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(AKMIN,MAXGRDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ALBEDO,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(XKHKM,MAXGRDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ZKHKM,MAXGRDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(CSZ,MAXGRDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(CSX,MAXGRDS,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(LEVEL,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ICLOUD,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IRAIN,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IPRIS,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ISNOW,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IAGGR,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IGRAUP,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IHAIL,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(jnmb,7,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(RPARM,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(PPARM,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(SPARM,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(APARM,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(GPARM,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(HPARM,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(CPARM,1,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(GNU,7,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(vapour_on,1,MPI_LOGICAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(cloud_on,1,MPI_LOGICAL,mainnum,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(bulk_on,1,MPI_LOGICAL,mainnum,MPI_COMM_WORLD,ierr)

   return
end subroutine masterput_nl
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will broadcast grid-related initialisation information to the nodes.  !
!------------------------------------------------------------------------------------------!
subroutine masterput_gridinit(master_num)

   use mem_grid
   use rpara

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: master_num
   !----- Local variables -----------------------------------------------------------------!
   integer :: ierr
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Now all the intended variables are broadcast, they will be caught by the nodes.    !
   ! Here THE ORDER MATTERS A LOT. If you are including new variables, make sure that the  !
   ! variables listed here are in the very same order as the nodeget_gridinit subroutine   !
   ! (later in this file), otherwise horrible things will happen to your run.              !
   !---------------------------------------------------------------------------------------!
   call MPI_Bcast(NNX,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNX1,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNX2,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNY,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNY1,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNY2,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNZ,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNXYZP,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNXYSP,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNXYP,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(JDIM,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine masterput_gridinit
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will send most of the grid information to the nodes, including the    !
! sub-domain size, its relative position in the full domain, etc.                          !
!------------------------------------------------------------------------------------------!
subroutine masterput_grid_dimens(master_num)

   use mem_grid
   use cyclic_mod
   use rpara

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'mpif.h'
   include 'interface.h'
   !----- Arguments -----------------------------------------------------------------------!
   integer                      , intent(in) :: master_num
   !----- Local variables -----------------------------------------------------------------!
   integer                                   :: nm,zzz,nxpts,nypts,nzpts,ierr,nmiii
   integer, dimension(2,maxmach)             :: lbc_buffs_tmp
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Now all the intended variables are broadcast, they will be caught by the nodes.    !
   ! Here THE ORDER MATTERS A LOT. If you are including new variables, make sure that the  !
   ! variables listed here are in the very same order as the nodeget_grid_dimens subrou-   !
   ! tine (later in this file), otherwise horrible things will happen to your run.         !
   !---------------------------------------------------------------------------------------!
   do nm=1,nmachs
      do ngrid=1,ngrids
         nxpts=nxend(nm,ngrid)-nxbeg(nm,ngrid)+1
         nypts=nyend(nm,ngrid)-nybeg(nm,ngrid)+1
         nzpts=nnzp(ngrid)
         call MPI_Bcast(nxpts,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(nypts,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(nzpts,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(nxbegc(nm,ngrid),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(nxendc(nm,ngrid),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(nybegc(nm,ngrid),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(nyendc(nm,ngrid),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(ixoff(nm,ngrid),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(iyoff(nm,ngrid),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
         call MPI_Bcast(ibcflg(nm,ngrid),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
      end do
      call MPI_Bcast(machnum(nm),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   end do

   do nm=1,nmachs
      call MPI_Send(ibounds(1,1,nm),8*maxgrds,MPI_INTEGER,machnum(nm),21,MPI_COMM_WORLD    &
                   ,ierr)
      call MPI_Send(inode_paths_master(1,1,1,1,nm),5*7*maxgrds*maxmach,MPI_INTEGER         &
                   ,machnum(nm),22,MPI_COMM_WORLD,ierr)
      call MPI_Send(iget_paths_master(1,1,1,nm),6*maxgrds*maxmach,MPI_INTEGER,machnum(nm)  &
                   ,23,MPI_COMM_WORLD,ierr)
      if (npts_cyc > 0) then
         call MPI_Send(ipathst_cyc,8*npts_cyc,MPI_INTEGER,machnum(nm),24,MPI_COMM_WORLD    &
                      ,ierr)
         call MPI_Send(ipathsu_cyc,8*npts_cyc,MPI_INTEGER,machnum(nm),25,MPI_COMM_WORLD    &
                      ,ierr)
         call MPI_Send(ipathsv_cyc,8*npts_cyc,MPI_INTEGER,machnum(nm),26,MPI_COMM_WORLD    &
                      ,ierr)
      end if
      zzz=1000
      do nmiii=1,nmachs
        zzz=zzz+1
        call MPI_Send(lbc_buffs(1,nmiii,nm),1,MPI_INTEGER,machnum(nm),zzz,MPI_COMM_WORLD   &
                     ,ierr)
        zzz=zzz+1
        call MPI_Send(lbc_buffs(2,nmiii,nm),1,MPI_INTEGER,machnum(nm),zzz,MPI_COMM_WORLD   &
                     ,ierr)
      end do
      zzz=zzz+1
      call MPI_Send(newbuff_nest1(nm),1,MPI_INTEGER,machnum(nm),zzz,MPI_COMM_WORLD,ierr)
      zzz=zzz+1
      call MPI_Send(nbuff_nest1(nm),1,MPI_INTEGER,machnum(nm),zzz,MPI_COMM_WORLD,ierr)
   enddo
   return
end subroutine masterput_grid_dimens
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine sends the grid setting information (e.g. grid resolution, domain     !
! levels) to the nodes.                                                                    !
!------------------------------------------------------------------------------------------!
subroutine masterput_gridset(master_num)

   use mem_grid
   use rpara

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: master_num
   !----- Local variables -----------------------------------------------------------------!
   integer             :: ierr
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Now all the intended variables are broadcast, they will be caught by the nodes.    !
   ! Here THE ORDER MATTERS A LOT. If you are including new variables, make sure that the  !
   ! variables listed here are in the very same order as the nodeget_gridset subroutine    !
   ! (later in this file), otherwise horrible things will happen to your run.              !
   !---------------------------------------------------------------------------------------!
   call MPI_Bcast(nrz,nzpmax*maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ipm,nxpmax*maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(jpm,nypmax*maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(kpm,nzpmax*maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(xmn,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ymn,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(zmn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(xtn,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ytn,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ztn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(dzmn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(dzm2n,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(dztn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(dzt2n,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(deltaxn,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(deltayn,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(deltazn,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(ztop,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine masterput_gridset
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine sends information that did not fit in any of the previous calls, but  !
! it does the same thing: sends information to the nodes.                                  ! 
!------------------------------------------------------------------------------------------!
subroutine masterput_misc(master_num)

   use mem_grid
   use rpara
   use mem_cuparm
   use ref_sounding

   implicit none

   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Arguments -----------------------------------------------------------------------!
   integer                      , intent(in) :: master_num
   !----- Local variables -----------------------------------------------------------------!
   integer             :: ierr
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Now all the intended variables are broadcast, they will be caught by the nodes.    !
   ! Here THE ORDER MATTERS A LOT. If you are including new variables, make sure that the  !
   ! variables listed here are in the very same order as the nodeget_misc subroutine       !
   ! (later in this file), otherwise horrible things will happen to your run.              !
   !---------------------------------------------------------------------------------------!
   call MPI_Bcast(nsubs,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(itopo,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(impl,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iadvl,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iadvf,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(time,1,MPI_DOUBLE_PRECISION,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(u01dn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(v01dn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(pi01dn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(th01dn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(dn01dn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(rt01dn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(htn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(hwn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ht2n,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ht4n,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(hw2n,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(hw4n,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ht,nzpmax,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(hw,nzpmax,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ht2,nzpmax,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ht4,nzpmax,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(hw2,nzpmax,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(hw4,nzpmax,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(if_cuinv,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(tnudcu,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(wt_cu_grid,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(tcu_beg,1,MPI_DOUBLE_PRECISION,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(tcu_end,1,MPI_DOUBLE_PRECISION,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(cu_til,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(cu_tel,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine masterput_misc
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine sends nesting-related information to the nodes.                       ! 
!------------------------------------------------------------------------------------------!
subroutine masterput_cofnest(master_num)

   use mem_grid
   use rpara

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: master_num
   !----- Local variables -----------------------------------------------------------------!
   integer             :: ierr
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Now all the intended variables are broadcast, they will be caught by the nodes.    !
   ! Here THE ORDER MATTERS A LOT. If you are including new variables, make sure that the  !
   ! variables listed here are in the very same order as the nodeget_cofnest subroutine    !
   ! (later in this file), otherwise horrible things will happen to your run.              !
   !---------------------------------------------------------------------------------------!
   call MPI_Bcast(ei1,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ei2,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ei3,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ei4,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ei5,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ei6,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ei7,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(ej1,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ej2,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ej3,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ej4,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ej5,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ej6,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ej7,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(ek1,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ek2,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ek3,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ek4,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ek5,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ek6,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ek7,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(fbcf,nzpmax*maxgrds*4,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine masterput_cofnest
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine sends the microphysics variables that might have been read from the  !
! collection table or derived from this to the nodes.                                      !
!------------------------------------------------------------------------------------------!
subroutine masterput_micphys(master_num)

   use micphys
   use rpara

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: master_num
   !----- Local variables -----------------------------------------------------------------!
   integer             :: ierr
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Now all the intended variables are broadcast, they will be caught by the nodes.    !
   ! Here THE ORDER MATTERS A LOT. If you are including new variables, make sure that the  !
   ! variables listed here are in the very same order as the nodeget_micphys subroutine    !
   ! (later in this file), otherwise horrible things will happen to your run.              !
   !---------------------------------------------------------------------------------------!
   call MPI_Bcast(availcat,ncat,MPI_LOGICAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(progncat,ncat,MPI_LOGICAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(shapefac,nhcat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(cfmas,nhcat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(pwmas,nhcat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(cfvt,nhcat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(pwvt,nhcat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(parm,ncat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(emb0,ncat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(emb1,ncat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(emb2,nhcat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(rxmin,ncat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(cxmin,nhcat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(coltabc,nembc*nembc*npairc,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(coltabr,nembc*nembc*npairr,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(ipairc,nhcat*nhcat,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ipairr,nhcat*nhcat,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine masterput_micphys
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine sends some CARMA radiation variables to the nodes.                   !
!------------------------------------------------------------------------------------------!
subroutine masterput_carma(master_num)

   use mem_globrad, only: &
        nirp, ngauss, ntotal, nwave, nsol, nsolp, &
        o3mixp,pj,ako3,akco2,akh2o,contnm,gangle, &
        gratio,gweight,weight,treal,ttmag,nprob,  &
        pso2,pso3,psh2o,psco2,wave,solfx,xah2o, &
        xaco2,xao2,xao3,ta,tb,wa,wb,ga,gb,tia,tib, &
        wia,wib,gia,gib,alpha,gama,caseE,caseW,caseG, &
        sq3,jdble,jn,tpi,cpcon,fdegday
   use rpara

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: master_num
   !----- Local variables -----------------------------------------------------------------!
   integer             :: ierr
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Now all the intended variables are broadcast, they will be caught by the nodes.    !
   ! Here THE ORDER MATTERS A LOT. If you are including new variables, make sure that the  !
   ! variables listed here are in the very same order as the nodeget_carma subroutine      !
   ! (later in this file), otherwise horrible things will happen to your run.              !
   !---------------------------------------------------------------------------------------!
   call MPI_Bcast(o3mixp,6,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(pj,6,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ako3,4*6,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(akco2,6*6,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(akh2o,54*6,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(contnm,nirp,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gangle,ngauss,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gratio,ngauss,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gweight,ngauss,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(weight,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(treal,2*nwave,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ttmag,2*nwave,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nprob,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(pso2,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(pso3,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(psh2o,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(psco2,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(wave,(nwave+1),MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(solfx,nsol,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(xah2o,nsolp,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(xaco2,nsolp,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(xao2,nsolp,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(xao3,nsolp,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ta,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(tb,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(wa,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(wb,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ga,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gb,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(tia,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(tib,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(wia,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(wib,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gia,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gib,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(alpha,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gama,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(caseE,(9*nwave),MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(caseW,(9*nwave),MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(caseG,(9*nwave),MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(sq3,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(jdble,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(jn,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(cpcon,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(fdegday,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine masterput_carma
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This routine catches the basic processor ID sent from the head node.                  !
!------------------------------------------------------------------------------------------!
subroutine nodeget_processid(init)

   use grid_dims
   use node_mod

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in) :: init
   !----- Local variables -----------------------------------------------------------------!
   integer             :: ierr
   !---------------------------------------------------------------------------------------!

   if(init == 1) then
      
      call MPI_Recv(master_num,1,MPI_INTEGER,0,11,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call MPI_Recv(mchnum,1,MPI_INTEGER,0,12,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call MPI_Recv(mynum,1,MPI_INTEGER,0,13,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call MPI_Recv(nmachs,1,MPI_INTEGER,0,14,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call MPI_Recv(machs,nmachs,MPI_INTEGER,0,15,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

   endif
   write (unit=*,fmt='(a,1x,i5,1x,a)') ' ==== Node',mynum,'got first message!'

   return
end subroutine nodeget_processid
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This routine receives the namelist variables (and their closest friends) from the     !
! head node.                                                                               !
!------------------------------------------------------------------------------------------!
subroutine nodeget_nl


   use mem_all
   use node_mod
   use therm_lib          , only : level                       & ! intent(out)
                                  ,vapour_on                   & ! intent(out)
                                  ,cloud_on                    & ! intent(out)
                                  ,bulk_on                     ! ! intent(out)
   use mem_mass           , only : iexev                       & ! intent(out)
                                  ,imassflx                    ! ! intent(out)
   use grell_coms         , only:  closure_type                & ! intent(out)
                                  ,maxclouds                   & ! intent(out)
                                  ,iupmethod                   & ! intent(out)
                                  ,depth_min                   & ! intent(out)
                                  ,cap_maxs                    & ! intent(out)
                                  ,maxens_lsf                  & ! intent(out)
                                  ,maxens_dyn                  & ! intent(out)
                                  ,maxens_eff                  & ! intent(out)
                                  ,maxens_cap                  & ! intent(out)
                                  ,iupmethod                   & ! intent(out)
                                  ,iupstrm                     & ! intent(out)
                                  ,radius                      & ! intent(out)
                                  ,zkbmax                      & ! intent(out)
                                  ,max_heat                    & ! intent(out)
                                  ,zcutdown                    & ! intent(out)
                                  ,z_detr                      ! ! intent(out)
   use sib_vars           , only : n_co2                       & ! intent(out)
                                  ,co2_init                    ! ! intent(out)
   use catt_start         , only : catt                        ! ! intent(out)
   use mem_globrad        , only : raddatfn                    ! ! intent(out)
   use emission_source_map, only : plumerise                   ! ! intent(out)
   use plume_utils        , only : prfrq                       ! ! intent(out)
   use teb_spm_start      , only : teb_spm                     ! ! intent(out)
   use teb_vars_const     , only : iteb,tminbld,nteb           & ! intent(out)
                                  ,rushh1,rushh2,daylight      & ! intent(out)
                                  ,d_road,tc_road,hc_road      & ! intent(out)
                                  ,d_wall,tc_wall,hc_wall      & ! intent(out)
                                  ,d_roof,tc_roof,hc_roof      & ! intent(out)
                                  ,nurbtype,ileafcod,z0_town   & ! intent(out)
                                  ,bld,bld_height,bld_hl_ratio & ! intent(out)
                                  ,aroof,eroof,aroad,eroad     & ! intent(out)
                                  ,awall,ewall,htraf,hindu     & ! intent(out)
                                  ,pletraf,pleindu             ! ! intent(out)
   use mem_emiss          , only : ichemi                      & ! intent(out)
                                  ,ichemi_in                   & ! intent(out)
                                  ,chemdata_in                 & ! intent(out)
                                  ,isource                     & ! intent(out)
                                  ,weekdayin                   & ! intent(out)
                                  ,efsat                       & ! intent(out)
                                  ,efsun                       & ! intent(out)
                                  ,eindno                      & ! intent(out)
                                  ,eindno2                     & ! intent(out)
                                  ,eindpm                      & ! intent(out)
                                  ,eindco                      & ! intent(out)
                                  ,eindso2                     & ! intent(out)
                                  ,eindvoc                     & ! intent(out)
                                  ,eveino                      & ! intent(out)
                                  ,eveino2                     & ! intent(out)
                                  ,eveipm                      & ! intent(out)
                                  ,eveico                      & ! intent(out)
                                  ,eveiso2                     & ! intent(out)
                                  ,eveivoc                       ! intent(out)
   use ref_sounding, only:         maxsndg                     ! ! intent(out)

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Local variables -----------------------------------------------------------------!
   integer :: nm,ierr
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Now all the broadcasted variables are captured by each slave node.                 !
   ! Here THE ORDER MATTERS A LOT. If you are including new variables, make sure that the  !
   ! variables listed here are in the very same order as the masterput_nl subroutine       !
   ! (earlier in this file), otherwise horrible things will happen to your run.            !
   !---------------------------------------------------------------------------------------!
   call MPI_Bcast(TIMMAX,1,MPI_DOUBLE_PRECISION,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(if_adap,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(load_bal,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NGRIDS,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNXP,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNYP,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNZP,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NZG,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NZS,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NXTNEST,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IHTRAN,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DELTAX,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DELTAY,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DELTAZ,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DZRAT,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DZMAX,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ZZ,NZPMAX,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)


   call MPI_Bcast(IDELTAT,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NACOUST,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IMONTHA,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IDATEA,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IYEARA,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ITIMEA,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IMONTHZ,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IDATEZ,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IYEARZ,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ITIMEZ,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(CENTLAT,MAXGRDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(CENTLON,MAXGRDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(POLELAT,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(POLELON,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(PLATN,MAXGRDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(PLONN,MAXGRDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(NSTRATX,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NSTRATY,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NSTRATZ1,NZPMAX,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NSTRATZ2,NZPMAX,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NESTZ1,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NESTZ2,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NINEST,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NJNEST,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NKNEST,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(GRIDU,MAXGRDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(GRIDV,MAXGRDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNSTTOP,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNSTBOT,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)


   call MPI_Bcast(IOUTPUT,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(INITFLD,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(FRQPRT,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(FRQHIS,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(FRQANL,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(FRQLITE,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NLITE_VARS,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(XLITE(1:20),20,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(YLITE(1:20),20,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ZLITE(1:20),20,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   do nm = 1, nlite_vars
      call MPI_Bcast(LITE_VARS(nm),32,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
   enddo


   call MPI_Bcast(FRQMEAN,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(FRQBOTH,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(AVGTIM,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(INITIAL,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NUD_TYPE,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IF_ODA,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NUDLAT,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(TNUDLAT,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(TNUDTOP,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(TNUDCENT,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ZNUDTOP,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(WT_NUDGE_GRID,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(WT_NUDGE_UV,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(WT_NUDGE_TH,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(WT_NUDGE_PI,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(WT_NUDGE_RT,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NUD_COND,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(TCOND_BEG,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(TCOND_END,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(T_NUDGE_RC,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(WT_NUDGEC_GRID,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(IUPDSST,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ISSTFLG,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  
   call MPI_Bcast(IUPDNDVI,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NDVIFLG,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  
   call MPI_Bcast(RUNTYPE(1:16),16,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(DTLONG,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNDTRAT,maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NADDSC,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(NNQPARM,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NCLOUDS,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NDEEPEST,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NSHALLOWEST,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(WCLDBS,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(CONFRQ,MAXCLOUDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(CPTIME,MAXCLOUDS,MPI_DOUBLE_PRECISION,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IUPMETHOD,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IUPSTRM,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(RADIUS,MAXCLOUDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DEPTH_MIN,MAXCLOUDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr) 
   call MPI_Bcast(CAP_MAXS,MAXCLOUDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr) 
   call MPI_Bcast(ZKBMAX,MAXCLOUDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)   
   call MPI_Bcast(ZCUTDOWN,MAXCLOUDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr) 
   call MPI_Bcast(Z_DETR,MAXCLOUDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)   
   call MPI_Bcast(MAX_HEAT,MAXCLOUDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr) 
   do nm=1,maxclouds
      call MPI_Bcast(CLOSURE_TYPE(nm),2,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr) 
   end do
   call MPI_Bcast(MAXENS_LSF,maxclouds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(MAXENS_EFF,maxclouds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(MAXENS_CAP,maxclouds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(SLZ,NZGMAX,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(STGOFF,NZGMAX,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(SLMSTR,NZGMAX,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IDIFFK,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IBRUVAIS,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IHORGRAD,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IF_URBAN_CANOPY,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
  
   call MPI_Bcast(CATT,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(RADDATFN(1:160),160,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)


   if (CATT == 1) then
     call MPI_Bcast(PLUMERISE,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(PRFRQ,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   endif

   call MPI_Bcast(TEB_SPM,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   ! TEB_SPM
   if (TEB_SPM == 1) then
     call MPI_Bcast(RUSHH1,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(RUSHH2,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(DAYLIGHT,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(ITEB,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(TMINBLD,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(NTEB,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(TC_ROOF,MAXSTEB,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(D_ROOF,MAXSTEB,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(HC_ROOF,MAXSTEB,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(HC_ROAD,MAXSTEB,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(TC_ROAD,MAXSTEB,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(D_ROAD,MAXSTEB,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(HC_WALL,MAXSTEB,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(TC_WALL,MAXSTEB,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(D_WALL,MAXSTEB,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(NURBTYPE,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(ILEAFCOD,MAXUBTP,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(Z0_TOWN,MAXUBTP,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(BLD,MAXUBTP,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(BLD_HEIGHT,MAXUBTP,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(BLD_HL_RATIO,MAXUBTP,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(AROOF,MAXUBTP,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EROOF,MAXUBTP,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(AROAD,MAXUBTP,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EROAD,MAXUBTP,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(AWALL,MAXUBTP,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EWALL,MAXUBTP,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(AROAD,MAXUBTP,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(HTRAF,MAXUBTP,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(HINDU,MAXUBTP,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(PLETRAF,MAXUBTP,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(PLEINDU,MAXUBTP,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

     call MPI_Bcast(ICHEMI,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(ICHEMI_IN,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(CHEMDATA_IN(1:80),80,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(ISOURCE,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(WEEKDAYIN(1:3),3,MPI_CHARACTER,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EFSAT,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EFSUN,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EINDNO,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EINDNO2,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EINDPM,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EINDCO,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EINDSO2,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EINDVOC,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EVEINO,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EVEINO2,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EVEIPM,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EVEICO,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EVEISO2,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
     call MPI_Bcast(EVEIVOC,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   end if
  
   call MPI_Bcast(LSFLG,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IBND,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(JBND,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NPATCH,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NVEGPAT,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ISFCL,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(N_CO2,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(CO2_INIT,maxsndg,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
  
   call MPI_Bcast(LONRAD,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(CPHAS,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DISTIM,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NFPT,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DTHCON,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(DRTCON,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(SEATMP,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(PCTLCON,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NSLCON,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NVGCON,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(RADFRQ,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ZROUGH,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ILWRTYP,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ISWRTYP,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ICUMFDBK,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ICORFLG,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IEXEV,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IMASSFLX,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(AKMIN,MAXGRDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ALBEDO,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(XKHKM,MAXGRDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ZKHKM,MAXGRDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(CSZ,MAXGRDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(CSX,MAXGRDS,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(LEVEL,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ICLOUD,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IRAIN,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IPRIS,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ISNOW,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IAGGR,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IGRAUP,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(IHAIL,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(jnmb,7,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(RPARM,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(PPARM,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(SPARM,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(APARM,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(GPARM,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(HPARM,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(CPARM,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(GNU,7,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(vapour_on,1,MPI_LOGICAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(cloud_on,1,MPI_LOGICAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(bulk_on,1,MPI_LOGICAL,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine nodeget_nl
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine receives grid-related initialisation information from the head node.  !
!------------------------------------------------------------------------------------------!
subroutine nodeget_gridinit

   use mem_grid
   use node_mod

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Local variables -----------------------------------------------------------------!
   integer :: nm,ierr
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Now all the broadcasted variables are captured by each slave node.                 !
   ! Here THE ORDER MATTERS A LOT. If you are including new variables, make sure that the  !
   ! variables listed here are in the very same order as the masterput_gridinit subroutine !
   ! (earlier in this file), otherwise horrible things will happen to your run.            !
   !---------------------------------------------------------------------------------------!
   call MPI_Bcast(NNX,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNX1,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNX2,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNY,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNY1,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNY2,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNZ,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNXYZP,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNXYSP,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(NNXYP,MAXGRDS,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(JDIM,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine nodeget_gridinit
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine receives most of the grid information from the head node, including   !
! the sub-domain size, its relative position in the full domain, etc.                      !
!------------------------------------------------------------------------------------------!
subroutine nodeget_grid_dimens()

   use mem_grid
   use node_mod
   use cyclic_mod 

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Local variables -----------------------------------------------------------------!
   integer                              :: nm,ierr,ng,zzz
   integer, allocatable, dimension(:,:) :: node_buffs_tmp
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Now all the broadcasted variables are captured by each slave node.                 !
   ! Here THE ORDER MATTERS A LOT. If you are including new variables, make sure that the  !
   ! variables listed here are in the very same order as the masterput_grid_dimens subrou- !
   ! tine (earlier in this file), otherwise horrible things will happen to your run.       !
   !---------------------------------------------------------------------------------------!
   do nm=1,nmachs
     do ng=1,ngrids
       call MPI_Bcast(nodemxp(nm,ng),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(nodemyp(nm,ng),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(nodemzp(nm,ng),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)


       call MPI_Bcast(nodeia(nm,ng),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(nodeiz(nm,ng),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(nodeja(nm,ng),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(nodejz(nm,ng),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(nodei0(nm,ng),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(nodej0(nm,ng),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(nodeibcon(nm,ng),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
     end do
     call MPI_Bcast(machs(nm),1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   end do


   call MPI_Recv(nodebounds,8*maxgrds,MPI_INTEGER,master_num,21,MPI_COMM_WORLD             &
                ,MPI_STATUS_IGNORE,ierr)
   call MPI_Recv(ipaths,5*7*maxgrds*maxmach,MPI_INTEGER,master_num,22,MPI_COMM_WORLD       &
                ,MPI_STATUS_IGNORE,ierr)
   call MPI_Recv(iget_paths,6*maxgrds*maxmach,MPI_INTEGER,master_num,23,MPI_COMM_WORLD     &
                ,MPI_STATUS_IGNORE,ierr)
   if(npts_cyc > 0) then
      call MPI_Recv(ipathst_cyc,8*npts_cyc,MPI_INTEGER,master_num,24,MPI_COMM_WORLD        &
                   ,MPI_STATUS_IGNORE,ierr)
      call MPI_Recv(ipathsu_cyc,8*npts_cyc,MPI_INTEGER,master_num,25,MPI_COMM_WORLD        &
                   ,MPI_STATUS_IGNORE,ierr)
      call MPI_Recv(ipathsv_cyc,8*npts_cyc,MPI_INTEGER,master_num,26,MPI_COMM_WORLD        &
                   ,MPI_STATUS_IGNORE,ierr)
   end if

   zzz=1000
   do nm=1,nmachs
     zzz=zzz+1
     call MPI_Recv(node_buffs(nm)%nsend,1,MPI_INTEGER,master_num,zzz,MPI_COMM_WORLD        &
                  ,MPI_STATUS_IGNORE,ierr)
     zzz=zzz+1
     call MPI_Recv(node_buffs(nm)%nrecv,1,MPI_INTEGER,master_num,zzz,MPI_COMM_WORLD        &
                  ,MPI_STATUS_IGNORE,ierr)
   end do
   zzz=zzz+1
   call MPI_Recv(newbuff_nest,1,MPI_INTEGER,master_num,zzz,MPI_COMM_WORLD                  &
                ,MPI_STATUS_IGNORE,ierr)
   zzz=zzz+1
   call MPI_Recv(nbuff_nest,1,MPI_INTEGER,master_num,zzz,MPI_COMM_WORLD,MPI_STATUS_IGNORE  &
                ,ierr)

   !----- Assigning the node_mod shortcuts with the appropiate values. --------------------!
   do ng=1,ngrids
      mmxp(ng)=nodemxp(mynum,ng)
      mmyp(ng)=nodemyp(mynum,ng)
      mmzp(ng)=nodemzp(mynum,ng)
      mia(ng)=nodeia(mynum,ng)
      miz(ng)=nodeiz(mynum,ng)
      mja(ng)=nodeja(mynum,ng)
      mjz(ng)=nodejz(mynum,ng)
      mi0(ng)=nodei0(mynum,ng)
      mj0(ng)=nodej0(mynum,ng)
      mibcon(ng)=nodeibcon(mynum,ng)
   end do
   write (unit=*,fmt='(a)') '---------------------------------------------------------'
   write (unit=*,fmt='(a,1x,i5)') 'In nodeget_grid_dimens, mynum=',mynum
   write (unit=*,fmt='(a,1x,8(i5,1x))') 'mmzp=',mmzp(1:ngrids) 
   write (unit=*,fmt='(a,1x,8(i5,1x))') 'mmxp=',mmxp(1:ngrids) 
   write (unit=*,fmt='(a,1x,8(i5,1x))') 'mmyp=',mmyp(1:ngrids) 
   write (unit=*,fmt='(a)') '---------------------------------------------------------'


   return
end subroutine nodeget_grid_dimens
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine receives the grid setting information (e.g. grid resolution, domain  !
! levels) from the head node.                                                              !
!------------------------------------------------------------------------------------------!
subroutine nodeget_gridset

   use mem_grid
   use node_mod
   implicit none

   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Local variables -----------------------------------------------------------------!
   integer             :: ierr
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Now all the broadcasted variables are captured by each slave node.                 !
   ! Here THE ORDER MATTERS A LOT. If you are including new variables, make sure that the  !
   ! variables listed here are in the very same order as the masterput_gridset subroutine  !
   ! (earlier in this file), otherwise horrible things will happen to your run.            !
   !---------------------------------------------------------------------------------------!
   call MPI_Bcast(nrz,nzpmax*maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ipm,nxpmax*maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(jpm,nypmax*maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(kpm,nzpmax*maxgrds,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(xmn,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ymn,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(zmn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(xtn,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ytn,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ztn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(dzmn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(dzm2n,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(dztn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(dzt2n,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(deltaxn,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(deltayn,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(deltazn,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(ztop,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine nodeget_gridset
!==========================================================================================!
!==========================================================================================!
!    This subroutine receives information that did not fit in any of the previous calls,   !
! but it does the same thing: receives data from the head node.                            ! 
!------------------------------------------------------------------------------------------!
subroutine nodeget_misc

   use mem_grid
   use mem_cuparm
   use ref_sounding
   use node_mod, only : master_num

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Local variables -----------------------------------------------------------------!
   integer             :: ierr
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Now all the broadcasted variables are captured by each slave node.                 !
   ! Here THE ORDER MATTERS A LOT. If you are including new variables, make sure that the  !
   ! variables listed here are in the very same order as the masterput_misc subroutine     !
   ! (earlier in this file), otherwise horrible things will happen to your run.            !
   !---------------------------------------------------------------------------------------!
   call MPI_Bcast(nsubs,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(itopo,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(impl,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iadvl,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(iadvf,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(time,1,MPI_DOUBLE_PRECISION,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(u01dn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(v01dn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(pi01dn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(th01dn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(dn01dn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(rt01dn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(htn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(hwn,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ht2n,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ht4n,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(hw2n,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(hw4n,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ht,nzpmax,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(hw,nzpmax,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ht2,nzpmax,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ht4,nzpmax,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(hw2,nzpmax,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(hw4,nzpmax,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(if_cuinv,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(tnudcu,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(wt_cu_grid,maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(tcu_beg,1,MPI_DOUBLE_PRECISION,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(tcu_end,1,MPI_DOUBLE_PRECISION,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(cu_til,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(cu_tel,1,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine nodeget_misc
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine receives nesting-related information from the head node.              ! 
!------------------------------------------------------------------------------------------!
subroutine nodeget_cofnest

   use mem_grid
   use node_mod, only : master_num

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Local variables -----------------------------------------------------------------!
   integer             :: ierr
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Now all the broadcasted variables are captured by each slave node.                 !
   ! Here THE ORDER MATTERS A LOT. If you are including new variables, make sure that the  !
   ! variables listed here are in the very same order as the masterput_cofnest subroutine  !
   ! (earlier in this file), otherwise horrible things will happen to your run.            !
   !---------------------------------------------------------------------------------------!
   call MPI_Bcast(ei1,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ei2,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ei3,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ei4,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ei5,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ei6,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ei7,nxpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(ej1,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ej2,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ej3,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ej4,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ej5,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ej6,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ej7,nypmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(ek1,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ek2,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ek3,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ek4,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ek5,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ek6,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ek7,nzpmax*maxgrds,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(fbcf,nzpmax*maxgrds*4,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine nodeget_cofnest
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine receives the microphysics variables that might have been read from   !
! the collection table or derived from this from the head node.                            !
!------------------------------------------------------------------------------------------!
subroutine nodeget_micphys

   use micphys
   use node_mod, only : master_num

   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Local variables -----------------------------------------------------------------!
   integer             :: ierr
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Now all the broadcasted variables are captured by each slave node.                 !
   ! Here THE ORDER MATTERS A LOT. If you are including new variables, make sure that the  !
   ! variables listed here are in the very same order as the masterput_micphys subroutine  !
   ! (earlier in this file), otherwise horrible things will happen to your run.            !
   !---------------------------------------------------------------------------------------!
   call MPI_Bcast(availcat,ncat,MPI_LOGICAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(progncat,ncat,MPI_LOGICAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(shapefac,nhcat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(cfmas,nhcat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(pwmas,nhcat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(cfvt,nhcat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(pwvt,nhcat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   call MPI_Bcast(parm,ncat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(emb0,ncat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(emb1,ncat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(emb2,nhcat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(rxmin,ncat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(cxmin,nhcat,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(coltabc,nembc*nembc*npairc,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(coltabr,nembc*nembc*npairr,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ipairc,nhcat*nhcat,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ipairr,nhcat*nhcat,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine nodeget_micphys
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine receives some CARMA radiation variables from the nodes.              !
!------------------------------------------------------------------------------------------!
subroutine nodeget_carma()


   use mem_globrad, only: &
        nirp, ngauss, ntotal, nwave, nsol, nsolp, &
        o3mixp,pj,ako3,akco2,akh2o,contnm,gangle, &
        gratio,gweight,weight,treal,ttmag,nprob,  &
        pso2,pso3,psh2o,psco2,wave,solfx,xah2o, &
        xaco2,xao2,xao3,ta,tb,wa,wb,ga,gb,tia,tib, &
        wia,wib,gia,gib,alpha,gama,caseE,caseW,caseG, &
        sq3,jdble,jn,tpi,cpcon,fdegday
   use rpara
   use node_mod, only: master_num
   implicit none
   !----- External variable declaration ---------------------------------------------------!
   include 'interface.h'
   include 'mpif.h'
   !----- Local variables -----------------------------------------------------------------!
   integer             :: ierr
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Now all the broadcasted variables are captured by each slave node.                 !
   ! Here THE ORDER MATTERS A LOT. If you are including new variables, make sure that the  !
   ! variables listed here are in the very same order as the masterput_carma subroutine    !
   ! (earlier in this file), otherwise horrible things will happen to your run.            !
   !---------------------------------------------------------------------------------------!
   call MPI_Bcast(o3mixp,6,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(pj,6,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ako3,4*6,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(akco2,6*6,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(akh2o,54*6,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(contnm,nirp,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gangle,ngauss,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gratio,ngauss,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gweight,ngauss,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(weight,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(treal,2*nwave,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ttmag,2*nwave,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(nprob,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(pso2,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(pso3,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(psh2o,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(psco2,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(wave,(nwave+1),MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(solfx,nsol,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(xah2o,nsolp,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(xaco2,nsolp,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(xao2,nsolp,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(xao3,nsolp,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ta,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(tb,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(wa,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(wb,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(ga,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gb,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(tia,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(tib,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(wia,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(wib,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gia,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gib,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(alpha,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(gama,ntotal,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(caseE,(9*nwave),MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(caseW,(9*nwave),MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(caseG,(9*nwave),MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(sq3,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(jdble,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(jn,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(cpcon,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
   call MPI_Bcast(fdegday,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)

   return
end subroutine nodeget_carma
!==========================================================================================!
!==========================================================================================!
