!============================= Change Log =================================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This sub-routine is the main time-step driver in BRAMS.                             !
!------------------------------------------------------------------------------------------!
subroutine timestep()

   use mem_basic          , only : basic_g         ! ! intent(inout)
   use node_mod           , only : mzp             & ! intent(in)
                                 , mxp             & ! intent(in)
                                 , myp             & ! intent(in)
                                 , ia              & ! intent(in)
                                 , iz              & ! intent(in)
                                 , ja              & ! intent(in)
                                 , jz              & ! intent(in)
                                 , i0              & ! intent(in)
                                 , j0              & ! intent(in)
                                 , izu             & ! intent(in)
                                 , jzv             & ! intent(in)
                                 , mynum           & ! intent(in)
                                 , ibcon           & ! intent(in)
                                 , ipara           ! ! intent(in)
   use mem_cuparm         , only : nnqparm         & ! intent(in)
                                 , if_cuinv        ! ! intent(in)
   use mem_varinit        , only : nud_type        ! ! intent(in)
   use mem_turb           , only : if_urban_canopy & ! intent(in)
                                 , ihorgrad        ! ! intent(in)
   use mem_oda            , only : if_oda          ! ! intent(in)
   use mem_grid           , only : ngrid           & ! intent(in)
                                 , time            & ! intent(in)
                                 , dtlong          & ! intent(in)
                                 , dtlongn         & ! intent(in)
                                 , iyeara          & ! intent(in)
                                 , imontha         & ! intent(in)
                                 , idatea          & ! intent(in)
                                 , grid_g          & ! intent(inout)
                                 , nxtnest         & ! intent(in)
                                 , if_adap         & ! intent(in)
                                 , dtlt            & ! intent(in)
                                 , istp            & ! intent(in)
                                 , jdim            & ! intent(in)
                                 , nzp             & ! intent(in)
                                 , f_thermo_e      & ! intent(in)
                                 , f_thermo_w      & ! intent(in)
                                 , f_thermo_s      & ! intent(in)
                                 , f_thermo_n      ! ! intent(in)
   use mem_scalar         , only : scalar_g        ! ! intent(in)
   use mem_leaf           , only : isfcl           ! ! intent(in)
   use catt_start         , only : catt            ! ! intent(in)
   use emission_source_map, only : burns           ! ! sub-routine
   use teb_spm_start      , only : TEB_SPM         ! ! intent(in)
   use mem_emiss          , only : ichemi          & ! intent(in)
                                 , isource         ! ! intent(in)
   use therm_lib          , only : bulk_on         ! ! intent(in)
   use advect_kit         , only : calc_advec      ! ! sub-routine
   use mem_mass           , only : iexev           & ! intent(in)
                                 , imassflx        ! ! intent(in)
   use mem_mnt_advec      , only : iadvec          ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   real              :: t1
   real              :: w1
   !----- External functions. -------------------------------------------------------------!
   real   , external :: cputime
   !----- Locally saved variables. --------------------------------------------------------!
   integer, save :: ncall=0
   !----- Local constants. ----------------------------------------------------------------!
   logical, parameter :: acct     = .false. ! To Activate acctimes
   !---------------------------------------------------------------------------------------!



   if (acct) call acctimes('init',0,' ',t1,w1)


   !------ Zero out all tendency arrays. --------------------------------------------------!
   t1 = cputime(w1)
   call tend0()          
   if (acct) call acctimes('accu',1,'TEND0',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Update the biomass burning fields if this is a CATT run.                          !
   !---------------------------------------------------------------------------------------!
   if (CATT == 1) then
      t1 = cputime(w1)
      call burns(ngrid,mzp,mxp,myp,ia,iz,ja,jz,scalar_g,time,iyeara,imontha,idatea)
      if (acct) call acctimes('accu',2,'BURNS',t1,w1)
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !  Thermodynamic diagnosis                                                              !
   !---------------------------------------------------------------------------------------!
   if (.not. bulk_on) then
      t1 = cputime(w1)
      call thermo(mzp,mxp,myp,1,mxp,1,myp)
      if (acct) call acctimes('accu',3,'THERMO',t1,w1)
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Medvigy's mass conservation fix.                                                  !
   !---------------------------------------------------------------------------------------!
   if (iexev == 2) then
      t1 = cputime(w1)
      call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'ADV')
      if (acct) call acctimes('accu',4,'EXEVOLVE_ADV',t1,w1)
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Radiation.                                                                           !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   call radiate(mzp,mxp,myp,ia,iz,ja,jz,mynum) 
   if (acct) call acctimes('accu',5,'RADIATE',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Surface layer, soil and vegetation or ecosystem demography model.                 !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   select case (isfcl)
   case (0:2,4)
      call leaf3_timestep()
   case (5)
      call ed_timestep()
   end select
   if (acct) call acctimes('accu',6,'SFCLYR',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     The dry deposition driver in case this is a CATT run.                             !
   !---------------------------------------------------------------------------------------!
   if (CATT==1) then
      t1 = cputime(w1)
      call drydep_driver(mzp,mxp,myp,ia,iz,ja,jz)
      if (acct) call acctimes('accu',7,'DryDep',t1,w1)
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Send boundaries to adjoining nodes.                                                  !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   call mpilbc_driver('sendlbc',1)
   if (acct) call acctimes('accu',8,'LBC_1st',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Coriolis terms.                                                                      !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   call corlos(mzp,mxp,myp,i0,j0,ia,iz,ja,jz,izu,jzv) 
   if (acct) call acctimes('accu',9,'CORLOS',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Velocity advection.                                                                  !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   call advectc('V',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
   if (acct) call acctimes('accu',10,'ADVECTv',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Urban canopy parameterization.                                                       !
   !---------------------------------------------------------------------------------------!
   if (if_urban_canopy == 1) then
      t1 = cputime(w1)
      call urban_canopy()      
      if (acct) call acctimes('accu',11,'URBAN_CANOPY',t1,w1)
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Analysis nudging and boundary condition.                                             !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   if (nud_type > 0) call datassim()  
   if (acct) call acctimes('accu',12,'DATASSIM',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Observation data assimilation.                                                       !
   !---------------------------------------------------------------------------------------!
   if (if_oda == 1) then
      t1 = cputime(w1)
      call oda_nudge()  
      if (acct) call acctimes('accu',13,'ODA_NUDGE',t1,w1)
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Nested grid boundaries.                                                              !
   !---------------------------------------------------------------------------------------!
   if (nxtnest(ngrid) >= 1) then
      t1 = cputime(w1)
      call nstbdriv()  
      if (acct) call acctimes('accu',14,'NSTBDRIV',t1,w1)
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Rayleigh friction for theta.                                                         !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   call rayft()           
   if (acct) call acctimes('accu',15,'RAYFT',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Get the overlap region between parallel nodes.                                    !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   call mpilbc_driver('getlbc',1)
   if (acct) call acctimes('accu',16,'LBC_2nd',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! Exner function correction.                                                            !
   !---------------------------------------------------------------------------------------!
   if (iexev == 2) then
      t1 = cputime(w1)
      call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'THA')
      if (acct) call acctimes('accu',17,'EXEVOLVE_THA',t1,w1)
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Sub-grid diffusion terms.                                                            !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   if (if_adap == 0 .and. ihorgrad == 2) then
      call diffuse_brams31()
   else
      call diffuse()
   end if
   if (acct) call acctimes('accu',18,'DIFFUSE',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Thermodynamic advection.                                                             !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   select case (iadvec)
   case (1)
      call advectc('T',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
   case (2)
      call radvc_mnt_driver(mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
   end select
   if (acct) call acctimes('accu',19,'ADVECTs',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Save mass flux into permanent arrays for use by Lagrangian models.                !
   !---------------------------------------------------------------------------------------!
   if (imassflx == 1) then
     t1 = cputime(w1)
     call prep_advflx_to_mass(mzp,mxp,myp,ia,iz,ja,jz,ngrid)
     if (acct) call acctimes('accu',20,'ADVEC_TO_MASS',t1,w1)
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Main driver for cumulus parametrisation.  The sub-routine will decide which       !
   ! closure shall be used.  This must be called after the tendencies of CO2, ice-liquid   !
   ! potential temperature, mixing ratio, and TKE have been updated because the convection !
   ! scheme uses the tendency of these variables.                                          !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   call rconv_driver()
   if (acct) call acctimes('accu',21,'CUPARM',t1,w1)
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Update urban emissions.                                                           !
   !---------------------------------------------------------------------------------------!
   if (teb_spm == 1) then
      if (isource == 1) then
         t1 = cputime(w1)
         call sources_teb(mzp,mxp,myp,ia,iz,ja,jz,ngrid,dtlt)
         if (acct) call acctimes('accu',22,'EMISS',t1,w1)
      end if

      !---- Update chemistry. -------------------------------------------------------------!
      if (ichemi==1) then
         t1 = cputime(w1)
         call ozone(mzp,mxp,myp,ia,iz,ja,jz,ngrid,dtlt)
         if (acct) call acctimes('accu',23,'OZONE',t1,w1)
      end if
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Update scalars.                                                                      !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   call predtr()
   call mpilbc_driver('fulllbc',1)
   if (acct) call acctimes('accu',24,'PREDTR',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Make moisture variables positive definite.                                           !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   call negadj1(mzp,mxp,myp,1,mxp,1,myp) 
   if (acct) call acctimes('accu',25,'NEGADJ1',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Microphysics.                                                                        !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   if (bulk_on) then
      call micro_driver()
   end if
   if (acct) call acctimes('accu',26,'MICRO',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !  Thermodynamic diagnosis.                                                             !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   if (.not. bulk_on) then
      call thermo(mzp,mxp,myp,1,mxp,1,myp) 
   end if
   if (acct) call acctimes('accu',27,'THERMO',t1,w1)
   !---------------------------------------------------------------------------------------!


   !----- Medvigy's mass conservation fix -------------------------------------------------!
   if (iexev == 2) then
       t1 = cputime(w1)
       call exevolve(mzp,mxp,myp,ngrid,ia,iz,ja,jz,izu,jzv,jdim,mynum,dtlt,'THS')
       if (acct) call acctimes('accu',28,'EXEVOLVE_THS',t1,w1)
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Apply scalar boundary conditions.                                                !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   call trsets()
   if (acct) call acctimes('accu',29,'TRSETS',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Lateral velocity boundaries - radiative.                                         !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   call latbnd()
   if (acct) call acctimes('accu',30,'LATBND',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      First stage of the Asselin filter.                                               !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   call hadvance(1)
   if (acct) call acctimes('accu',31,'HADVANCE1',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Buoyancy term for w equation.                                                    !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   call buoyancy()
   if (acct) call acctimes('accu',32,'BUOYANCY',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Acoustic small timesteps.                                                        !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   call acoustic_new()
   if (acct) call acctimes('accu',33,'ACOUSTIC',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Last stage of Asselin filter.                                                    !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   call hadvance(2)
   if (acct) call acctimes('accu',34,'HADVANCE2',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Velocity/pressure boundary conditions.                                           !
   !---------------------------------------------------------------------------------------!
   t1 = cputime(w1)
   call vpsets()          
   if (acct) call acctimes('accu',35,'VPSETS',t1,w1)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Print the table with the time spent by each step.                                 !
   !---------------------------------------------------------------------------------------!
   if (acct) then
      if (mod(istp,4) == 0) then
         call acctimes('prin',-1,' ',t1,w1)
      end if
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Call THERMO at the true boundaries.                                               !
   !---------------------------------------------------------------------------------------!
   call thermo_boundary_driver((time+dtlongn(ngrid)),dtlongn(ngrid),f_thermo_e(ngrid)      &
                              ,f_thermo_w(ngrid),f_thermo_s(ngrid),f_thermo_n(ngrid)       &
                              ,nzp,mxp,myp,jdim)
   !---------------------------------------------------------------------------------------!

   return
end subroutine timestep
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine is used to depict the time spent by each of the stages of one full  !
! time step.  This should be used for profiling only, and should be normally turned off    !
! when running "normal" simulations.                                                       !
!------------------------------------------------------------------------------------------!
subroutine acctimes(action,num,string,t1,w1)

   use mem_all
   use node_mod

   implicit none
   !----- Included commons. ---------------------------------------------------------------!
   include 'interface.h'
   !----- Local constants. ----------------------------------------------------------------!
   integer                               , parameter   :: num_times=100
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=*)                      , intent(in)  :: action
   character(len=*)                      , intent(in)  :: string
   integer                               , intent(in)  :: num
   real                                  , intent(in)  :: t1
   real                                  , intent(in)  :: w1
   !----- Local variables. ----------------------------------------------------------------!
   integer                                             :: npts
   integer                                             :: i
   integer                                             :: j
   real                                                :: sumtime
   real                                                :: pcpu
   real                                                :: cpuinc
   real                                                :: ww
   !----- External functions. -------------------------------------------------------------!
   real                                  , external    :: cputime
   real                                  , external    :: walltime
   real                                  , external    :: valugp
   !----- Locally saved variables. --------------------------------------------------------!
   character(len=8), dimension(num_times), save        :: crtimes
   real            , dimension(num_times), save        :: rtimes
   real            , dimension(num_times), save        :: wtimes
   !---------------------------------------------------------------------------------------!


   select case (trim(action))
   case ('init')
      crtimes(1:num_times)                = ' '
      rtimes (1:num_times)                = 0.
      wtimes (1:num_times)                = 0.
      basic_g(ngrid)%cputime(1:mxp,1:myp) = 0.

   case ('prin')
      if (mynum /= 2 .and. mynum /= 0) return
      sumtime = sum(rtimes)

      write(unit=*,fmt='(a,1x,f10.3,1x,a)') '======= total CPU =====',sumtime,'=========='
      do i=1,35
         write (unit=*,fmt='(a,1x,i4,2(1x,a),4(1x,f10.3))')                                &
              ' Timings:',mynum,'-',crtimes(i),rtimes(i),rtimes(i)/sumtime*100.,wtimes(i)  &
             ,wtimes(i)-rtimes(i)
      end do
      sumtime=sum(basic_g(ngrid)%cputime)

      write(unit=*,fmt='(a,2(a,1x,i5),a,1x,f10.5)') 'Total CPU secs','. Ngrid:',ngrid      &
                                                   ,'. Mynum:',mynum,'. Secs:',sumtime

   case ('null')
      !----- Accumulate full times into tables. -------------------------------------------!
      pcpu         = (cputime(ww)-t1)
      rtimes (num) = rtimes(num) + pcpu
      crtimes(num) = string
      wtimes (num) = wtimes(num)+(ww-w1)

   case ('accu')
      !----- Accumulate full times into tables. -------------------------------------------!
      pcpu         = (cputime(ww)-t1)
      rtimes (num) = rtimes(num) + pcpu
      crtimes(num) = string
      wtimes (num) = wtimes(num)+(ww-w1)

      !------------------------------------------------------------------------------------!
      !     Divide cpu time equally amoung columns and accumulate in                       !
      ! basic_g(ngrid)%cputime.                                                            !
      !------------------------------------------------------------------------------------!
      npts   = (iz-ia+1)*(jz-ja+1)
      cpuinc = pcpu / real(npts)
      do j=ja,jz
         do i=ia,iz
            basic_g(ngrid)%cputime(i,j) = basic_g(ngrid)%cputime(i,j) + cpuinc
         end do
      end do
   end select

   return
end subroutine acctimes
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will update the lateral boundary conditions.                         !
!------------------------------------------------------------------------------------------!
subroutine mpilbc_driver(action,izzflag)
   use node_mod, only : ipara ! ! intent(in)
   use mem_grid, only : ngrid ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=*), intent(in) :: action
   integer         , intent(in) :: izzflag
   !----- Local variables. ----------------------------------------------------------------!
   integer                      :: ist
   integer                      :: iadv
   !---------------------------------------------------------------------------------------!



   !----- Nothing to do if this is a serial run. ------------------------------------------!
   if (ipara /= 1) return
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Decide whether to send or receive 
   !---------------------------------------------------------------------------------------!
   select case(trim(action))
   case ('sendlbc')
      !------------------------------------------------------------------------------------!
      !     Send thermodynamic lateral conditions.                                         !
      !------------------------------------------------------------------------------------!
      call node_sendlbc()

      !------------------------------------------------------------------------------------!
      !     If it is the first grid, check whether we need to exchange cyclic conditions.  !
      !------------------------------------------------------------------------------------!
      if (ngrid  ==  1) call node_sendcyclic(1)
      !------------------------------------------------------------------------------------!

   case ('getlbc')
      !------------------------------------------------------------------------------------!
      !     Get thermodynamic lateral conditions.                                          !
      !------------------------------------------------------------------------------------!
      call node_getlbc()

      !------------------------------------------------------------------------------------!
      !     If it is the first grid, check whether we need to exchange cyclic conditions.  !
      !------------------------------------------------------------------------------------!
      if (ngrid  ==  1) call node_getcyclic(1)
      !------------------------------------------------------------------------------------!

   case ('fulllbc')
      !------------------------------------------------------------------------------------!
      !     Make a full exchange of the lateral conditions.                                !
      !------------------------------------------------------------------------------------!
      call node_sendlbc()
      call node_getlbc()



      !------------------------------------------------------------------------------------!
      !     If it is the first grid, check whether we need to exchange cyclic conditions.  !
      !------------------------------------------------------------------------------------!
      if (ngrid  ==  1) then
         call node_sendcyclic(1)
         call node_getcyclic(1)
      end if
      !------------------------------------------------------------------------------------!

   case ('sendst')
      !------------------------------------------------------------------------------------!
      !     Send staggered lateral conditions.  Here we also use the flag telling which    !
      ! field should be exchanged.                                                         !
      !------------------------------------------------------------------------------------!
      call node_sendst(izzflag)

      !------------------------------------------------------------------------------------!
      !     If it is the first grid, check whether we need to exchange cyclic conditions.  !
      !------------------------------------------------------------------------------------!
      if (ngrid  ==  1) call node_sendcyclic(izzflag)
      !------------------------------------------------------------------------------------!

   case ('getst')
      !------------------------------------------------------------------------------------!
      !     Get staggered lateral conditions.  Here we also use the flag telling which     !
      ! field should be exchanged.                                                         !
      !------------------------------------------------------------------------------------!
      call node_getst(izzflag)

      !------------------------------------------------------------------------------------!
      !     If it is the first grid, check whether we need to exchange cyclic conditions.  !
      !------------------------------------------------------------------------------------!
      if (ngrid  ==  1) call node_getcyclic(izzflag)
      !------------------------------------------------------------------------------------!

   case ('fullst')
      !------------------------------------------------------------------------------------!
      !     Full exchange of staggered lateral conditions.  Here we also use the flag      !
      ! telling which field should be exchanged.                                           !
      !------------------------------------------------------------------------------------!
      call node_sendst(izzflag)
      call node_getst (izzflag)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     If it is the first grid, check whether we need to exchange cyclic conditions.  !
      !------------------------------------------------------------------------------------!
      if (ngrid  ==  1) then
         call node_sendcyclic(izzflag)
         call node_getcyclic(izzflag)
      end if
      !------------------------------------------------------------------------------------!

   case ('sendadv')

      !------------------------------------------------------------------------------------!
      !      Send the advection boundary conditions variables.                             !
      !------------------------------------------------------------------------------------!
      select case (izzflag)
      case (0)
         !----- Send most variables (iaflag from 1 to 4). ---------------------------------!
         do iadv=1,4
            call node_sendadv(iadv)
         end do
      case default
         !----- Send only variabes corresponding to izzflag. ------------------------------!
         call node_sendadv(izzflag)
      end select
      !------------------------------------------------------------------------------------!

   case ('getadv')

      !------------------------------------------------------------------------------------!
      !      Get the advection boundary conditions variables.                              !
      !------------------------------------------------------------------------------------!
      select case (izzflag)
      case (0)
         !----- Get most variables (iaflag from 1 to 4). ----------------------------------!
         do iadv=1,4
            call node_getadv(iadv)
         end do
      case default
         !----- Get only variabes corresponding to izzflag. -------------------------------!
         call node_getadv(izzflag)
      end select
      !------------------------------------------------------------------------------------!

   case ('fulladv')

      !------------------------------------------------------------------------------------!
      !     Full exchange of advection lateral conditions.  Here we also use the flag      !
      ! telling which field should be exchanged.                                           !
      !------------------------------------------------------------------------------------!
      select case (izzflag)
      case (0)
         !----- Get most variables (iaflag from 1 to 4). ----------------------------------!
         do iadv=1,4
            call node_sendadv(iadv)
            call node_getadv (iadv)
         end do
      case default
         !----- Get only variabes corresponding to izzflag. -------------------------------!
         call node_sendadv(izzflag)
         call node_getadv (izzflag)
      end select
      !------------------------------------------------------------------------------------!
   end select

   return
end subroutine mpilbc_driver
!==========================================================================================!
!==========================================================================================!
