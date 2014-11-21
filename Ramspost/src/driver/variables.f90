!==========================================================================================!
!==========================================================================================!
subroutine RAMS_varlib(cvar,nx,ny,nz,nsl,npat,ncld,ngrd,flnm,cdname,cdunits,ivar_type &
                      ,a,b)
   use rout_coms   , only : undefflg        ! ! intent(in)
   use rpost_dims  , only : nwave           ! ! intent(in)
   use leaf_coms   , only : ustmin          & ! intent(in)
                          , ubmin           ! ! intent(in)
   use misc_coms   , only : memsize4        & ! intent(inout)
                          , ierr_getvar     & ! intent(inout)
                          , ifound          ! ! intent(inout)
   use scratch_coms, only : scr             & ! intent(inout)
                          , alloc_scratch   & ! subroutine
                          , zero_scratch    ! ! subroutine
   use micro_coms  , only : cfmas           & ! intent(in)
                          , pwmas           ! ! intent(in)
   use rconstants  , only : mmdryi          ! ! intent(in)
   use somevars    , only : co2_on          & ! intent(in)
                          , myco2con        & ! intent(in)
                          , mynzs           ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                       , intent(in)    :: nx
   integer                       , intent(in)    :: ny
   integer                       , intent(in)    :: nz
   integer                       , intent(in)    :: ngrd
   integer                       , intent(in)    :: nsl
   integer                       , intent(in)    :: npat
   integer                       , intent(in)    :: ncld
   integer                       , intent(out)   :: ivar_type
   character(len=*)              , intent(in)    :: cvar
   character(len=*)              , intent(in)    :: flnm
   character(len=*)              , intent(out)   :: cdname
   character(len=*)              , intent(out)   :: cdunits
   real            , dimension(*), intent(inout) :: a
   real            , dimension(*), intent(inout) :: b
   !----- Local variables. ----------------------------------------------------------------!
   integer                                       :: idim_type
   integer                                       :: irecind
   integer                                       :: irecsize
   integer                                       :: ind
   integer                                       :: ispec
   integer                                       :: ierr
   integer                                       :: ihyd
   real                                          :: fmult
   integer                       , save          :: memsave4 = 0
   !----- External functions. -------------------------------------------------------------!
   integer                       , external      :: RAMS_getvar
   !---------------------------------------------------------------------------------------!

   if (memsave4 == 0) then
      write (unit=*,fmt='(a,1x,i12)') '   - Allocating scratch variables; Size=',memsize4
      call alloc_scratch(scr,memsize4)
      memsave4 = memsize4
   elseif (memsize4 > memsave4) then
      write (unit=*,fmt='(a,1x,i12)') '   - Re-allocating scratch variables; New size='    &
                                     ,memsize4
      call alloc_scratch(scr,memsize4)
      memsave4 = memsize4
   end if

   !----- Initialise some variables. ------------------------------------------------------!
   ivar_type=0
   ierr_getvar=0
   ierr=0
   ifound=0


   !---------------------------------------------------------------------------------------!
   !     Reset all scratch variables.                                                      !
   !---------------------------------------------------------------------------------------!
   call zero_scratch(scr)
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !     The huge select case starts here...                                               !
   !---------------------------------------------------------------------------------------!
   select case (trim(cvar))
   case ('u','ue_avg','v','ve_avg')
      ivar_type   = 3
      ierr        = RAMS_getvar('UP',idim_type,ngrd,scr%u,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr        = RAMS_getvar('VP',idim_type,ngrd,scr%v,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_rotate(nx,ny,nz,scr%u,scr%v,ngrd)

      select case (trim(cvar))
      case ('u','ue_avg')
         call RAMS_comp_clone(nx,ny,nz,a,scr%u)
         call RAMS_comp_avgu(nx,ny,nz,a)
         cdname      = 'True zonal wind'
         cdunits     = 'm/s'

      case ('v','ve_avg')
         call RAMS_comp_clone(nx,ny,nz,a,scr%v)
         call RAMS_comp_avgv(nx,ny,nz,a)
         cdname      = 'True meridional wind'
         cdunits     = 'm/s'
      end select

   case ('w','w_avg')
      ivar_type   = 3
      ierr        = RAMS_getvar('WP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_avgw(nx,ny,nz,a)
      cdname      = 'True vertical velocity'
      cdunits     = 'm/s'

   case ('speed','direction')
      ivar_type   = 3
      ierr        = RAMS_getvar('UP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr        = RAMS_getvar('VP',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      select case (trim(cvar))
      case ('speed')
         call RAMS_comp_speed(nx,ny,nz,a,scr%c)
         cdname      = 'total wind speed'
         cdunits     = 'm/s'
      case ('direction')
         call RAMS_comp_dir(nx,ny,nz,a,scr%c,ngrd)
         cdname='direction'
         cdunits='deg'
      end select

   case ('pvap')
      ivar_type=3
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_press(nx,ny,nz,scr%c)               ! c is pressure [hPa]
      ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)    ! a is mixing ratio
      call RAMS_comp_noneg(nx,ny,nz,a)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_pvap(nx,ny,nz,scr%c,a)              ! a is vapour pressure [hPa]
      cdname='water vapour pressure'
      cdunits='hPa'

   case ('alpha')
      ivar_type=3
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm) ! c is Exner function
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,a,b,flnm)  ! a is potential temperature [K]
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,a,scr%c)              ! a is temperature [K]
      ierr= RAMS_getvar('RV',idim_type,ngrd,scr%d,b,flnm) ! d is mixing ratio [kg/kg]
      call RAMS_comp_noneg(nx,ny,nz,scr%d)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_press(nx,ny,nz,scr%c)                ! c is pressure [hPa]
      call RAMS_comp_mults(nx,ny,nz,scr%c,100.)           ! c is pressure [Pa]
      call RAMS_comp_pvap(nx,ny,nz,scr%c,scr%d)           ! d is vapour pressure [hPa]
      call RAMS_comp_spvol(nx,ny,nz,a,scr%d,scr%c)        ! a is specific volume [m3/kg]
      cdname='specific volume'
      cdunits='m3/kg'

   case ('solenoidx','solenoidy')
      ivar_type=3
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)    ! c is Exner function
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%e,b,flnm) ! e is pot. temperature [K]
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,scr%e,scr%c)             ! e is temperature [K]
      ierr= RAMS_getvar('RV',idim_type,ngrd,scr%d,b,flnm)    ! d is mixing ratio [kg/kg]
      call RAMS_comp_noneg(nx,ny,nz,scr%d)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_press(nx,ny,nz,scr%c)                   ! c is pressure [hPa]
      call RAMS_comp_mults(nx,ny,nz,scr%c,100.)              ! c is pressure [Pa]
      call RAMS_comp_pvap(nx,ny,nz,scr%c,scr%d)              ! d is vapour pressure [hPa]
      call RAMS_comp_spvol(nx,ny,nz,scr%e,scr%d,scr%c)       ! e is specific volume [m3/kg]
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%h,b,flnm)  ! h is topography [m]
      ierr_getvar = ierr_getvar + ierr
      select case (trim(cvar))
      case ('solenoidx')
         call RAMS_comp_solenoidx(nx,ny,nz,scr%e,scr%c,a,scr%h,ngrd)
         cdname='x-solenoid term'
         cdunits='rad/s2'
      case ('solenoidy')
         call RAMS_comp_solenoidy(nx,ny,nz,scr%e,scr%c,a,scr%h,ngrd)
         cdname='y-solenoid term'
         cdunits='rad/s2'
      end select

   case ('relvortx','relvorty')
      ivar_type=3

      ierr= RAMS_getvar('UP',idim_type,ngrd,scr%u,b,flnm)      ! u is U
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%u,1.e-6)
      ierr_getvar = ierr_getvar + ierr

      ierr= RAMS_getvar('VP',idim_type,ngrd,scr%v,b,flnm)      ! v is V
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%v,1.e-6)
      ierr_getvar = ierr_getvar + ierr

      ierr= RAMS_getvar('WP',idim_type,ngrd,scr%w,b,flnm)      ! w is W
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%w,1.e-6)
      ierr_getvar = ierr_getvar + ierr

      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%h,b,flnm)    ! h is height
      ierr_getvar = ierr_getvar + ierr

      ! x is the zonal vorticity and y is the meridional vorticity
      call RAMS_comp_relvortx(nx,ny,nz,scr%x,scr%v,scr%w,scr%m,scr%n,scr%o,scr%h,ngrd)
      call RAMS_comp_relvorty(nx,ny,nz,scr%y,scr%u,scr%w,scr%m,scr%n,scr%o,scr%h,ngrd)


      select case (trim(cvar))
      case ('relvortx')
         call RAMS_comp_clone (nx,ny,nz,a,scr%x)
         call RAMS_comp_rotate(nx,ny,nz,a,scr%y,ngrd)
         call RAMS_comp_avgu  (nx,ny,nz,a)
         cdname  = 'x-vorticity'
         cdunits = '1/s'
      case ('relvorty')
         call RAMS_comp_clone (nx,ny,nz,a,scr%y)
         call RAMS_comp_rotate(nx,ny,nz,scr%x,a,ngrd)
         call RAMS_comp_avgv  (nx,ny,nz,a)
         cdname  = 'y-vorticity'
         cdunits = '1/s'
      end select

   case ('relvortz','absvortz')
      ivar_type=3 

      ierr= RAMS_getvar('UP',idim_type,ngrd,scr%u,b,flnm) ! u is u
      call RAMS_flush_to_zero(nx,ny,nz,1,a,1.e-6)
      ierr_getvar = ierr_getvar + ierr

      ierr= RAMS_getvar('VP',idim_type,ngrd,scr%v,b,flnm) ! v is v
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%v,1.e-6)
      ierr_getvar = ierr_getvar + ierr

      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%h,b,flnm) ! h is topo
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_relvortz(nx,ny,nz,a,scr%u,scr%v,scr%m,scr%n,scr%o,scr%h,ngrd)

      select case (trim(cvar))
      case ('relvortz')
         cdname='relative z-vorticity'
         cdunits='1/s'
      case ('absvortz')

         ierr= RAMS_getvar('GLAT',idim_type,ngrd,scr%l,b,flnm) ! l is latitude
         ierr_getvar = ierr_getvar + ierr
         ! a is absolute vorticity
         call RAMS_comp_totvortz(nx,ny,nz,a,scr%l)

         cdname='absolute z-vorticity'
         cdunits='1/s'
      end select

   case ('potvort')
      ivar_type=3

      ierr= RAMS_getvar('UP',idim_type,ngrd,scr%u,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%u,1.e-6)      ! u is U


      ierr= RAMS_getvar('VP',idim_type,ngrd,scr%v,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%v,1.e-6)      ! v is V

      ierr= RAMS_getvar('WP',idim_type,ngrd,scr%w,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%w,1.e-6)      ! w is W

      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%h,b,flnm) ! h is topo
      ierr_getvar = ierr_getvar + ierr

      ierr= RAMS_getvar('GLAT',idim_type,ngrd,scr%l,b,flnm)  ! l is latitude
      ierr_getvar = ierr_getvar + ierr

      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%t,b,flnm) ! t is potential temperature
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_dn0(nx,ny,nz,scr%m,scr%n,scr%d,scr%o,ngrd)  ! d is density

      ! x, y, and z are the vorticity components. 
      call RAMS_comp_relvortx(nx,ny,nz,scr%x,scr%u,scr%w,scr%m,scr%n,scr%o,scr%h,ngrd)
      call RAMS_comp_relvortx(nx,ny,nz,scr%y,scr%v,scr%w,scr%m,scr%n,scr%o,scr%h,ngrd)
      call RAMS_comp_relvortx(nx,ny,nz,scr%z,scr%u,scr%v,scr%m,scr%n,scr%o,scr%h,ngrd)

      ! z is the absolute vertical vorticity
      call RAMS_comp_totvortz(nx,ny,nz,scr%z,scr%l)

      ! a is the potential vorticity
      call RAMS_comp_potvort(nx,ny,nz,a,scr%x,scr%y,scr%z,scr%t,scr%d,scr%m,scr%n,scr%o    &
                            ,scr%h,ngrd)

      cdname='potential vorticity'
      cdunits='1/s'

   case ('xydiv','moistdiv')
      ivar_type=3


      ierr= RAMS_getvar('UP',idim_type,ngrd,scr%u,b,flnm) ! u is U
      ierr_getvar = ierr_getvar + ierr
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%u,1.e-6)

      ierr= RAMS_getvar('VP',idim_type,ngrd,scr%v,b,flnm) ! v is V
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%v,1.e-6)
      ierr_getvar = ierr_getvar + ierr


      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%h,b,flnm) ! d is topo
      ierr_getvar = ierr_getvar + ierr


      select case(trim(cvar))
      case ('moistdiv')
         ierr= RAMS_getvar('RV',idim_type,ngrd,scr%r,b,flnm) ! r is rv
         ierr_getvar = ierr_getvar + ierr
         call RAMS_flush_to_zero(nx,ny,nz,1,scr%r,1.e-6)
         call RAMS_comp_mult(nx,ny,nz,scr%u,scr%r) ! x is rv*u
         call RAMS_comp_mult(nx,ny,nz,scr%v,scr%r) ! y is rv*v
      end select

      call RAMS_comp_xydiv(nx,ny,nz,a,scr%u,scr%v,scr%m,scr%n,scr%o,scr%h,ngrd)


      select case(trim(cvar))
      case ('moistdiv')
         cdname='Moisture divergence XY plane'
         cdunits='kg/kg/s'

      case ('xydiv')
         cdname='Divergence XY plane'
         cdunits='1/s'
      end select


   case ('mass_moistdiv')
      ivar_type=3


      ierr= RAMS_getvar('UP',idim_type,ngrd,scr%u,b,flnm) ! u is U
      ierr_getvar = ierr_getvar + ierr
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%u,1.e-6)

      ierr= RAMS_getvar('VP',idim_type,ngrd,scr%v,b,flnm) ! v is V
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%v,1.e-6)
      ierr_getvar = ierr_getvar + ierr


      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%h,b,flnm) ! h is topo
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_dn0(nx,ny,nz,scr%e,scr%t,scr%d,scr%h,ngrd) ! d is density

      ierr= RAMS_getvar('RV',idim_type,ngrd,scr%r,b,flnm) ! r is rv
      ierr_getvar = ierr_getvar + ierr
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%r,1.e-6)
      call RAMS_comp_mult(nz,ny,nz,scr%r,scr%d) ! r is rho*rv
      call RAMS_comp_mult(nx,ny,nz,scr%u,scr%r) ! x is u*rho*rv
      call RAMS_comp_mult(nx,ny,nz,scr%v,scr%r) ! y is v*rho*rv

      call RAMS_comp_xydiv(nx,ny,nz,a,scr%u,scr%v,scr%m,scr%n,scr%o,scr%h,ngrd)
      cdname='Moisture divergence XY plane'
      cdunits='kg/m3/s'

   case ('prwtr')
      ivar_type=2

      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,scr%c,scr%z,scr%d,scr%e,ngrd)

      ierr= RAMS_getvar('RTP',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_flush_to_zero(nx,ny,nz,1,scr%c,1.e-6)

      call RAMS_comp_prwtr(nx,ny,nz,a,scr%c,scr%d,scr%e,ngrd)

      cdname='Precipitable water'
      cdunits='kg/m2'

   case ('pi')
      ivar_type=3
      ierr= RAMS_getvar('PI',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Exner function'
      cdunits='J/(kg K)'

   case ('press')
      ivar_type=3
      ierr= RAMS_getvar('PI',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_press(nx,ny,nz,a)
      cdname='pressure'
      cdunits='mb'

   case ('theta')
      ivar_type=3
      ierr= RAMS_getvar('THETA',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='potential temp'
      cdunits='K'

   case ('thil')
      ivar_type=3
      ierr= RAMS_getvar('THP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='ice-liquid potential temp'
      cdunits='K'

   case ('co2')
      ivar_type=3
      if (co2_on /= 0) then
         ierr= RAMS_getvar('CO2P',idim_type,ngrd,a,b,flnm)
         ierr_getvar = ierr_getvar + ierr
      else
         write (unit=*,fmt='(a,1x,es12.5)') '       # Assigning constant CO2P =',myco2con(1)
         call ae0(nx*ny*nz,a,myco2con(1))
      end if
      cdname='carbon dioxide mixing ratio'
      cdunits='umol/mol'

   case ('dn0','pi0','th0')
      ivar_type=3
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%h,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,scr%e,scr%t,scr%d,scr%h,ngrd)
      select case (trim(cvar))
      case ('dn0')
         call RAMS_comp_clone(nx,ny,nz,a,scr%d)
         cdname='ref density'
         cdunits='kg/m3'
      case ('pi0')
         call RAMS_comp_clone(nx,ny,nz,a,scr%e)
         cdname='ref Exner func'
         cdunits='J/(kg K)'
      case ('th0')
         call RAMS_comp_clone(nx,ny,nz,a,scr%t)
         cdname='reference virtual potential temp'
         cdunits='K'
      end select

   case ('tempk','tempc')
      ivar_type=3
      ierr= RAMS_getvar('THETA',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,a,scr%c)
      select case (trim(cvar))
      case ('tempk')
         cdname='temperature'
         cdunits='K'
      case ('tempc')
         call RAMS_comp_tempC(nx,ny,nz,1,a)
         cdname='temperature'
         cdunits='C'
      end select


   case ('thetae_iv','theiv')
      ivar_type=3
      ierr= RAMS_getvar('THP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%f,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,scr%f,scr%c)
      call RAMS_comp_press(nx,ny,nz,scr%c)
      call RAMS_comp_mults(nx,ny,nz,scr%c,100.)
      ierr= RAMS_getvar('RV',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('RTP',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_thetaeiv(nx,ny,nz,a,scr%f,scr%c,scr%d,scr%e)
      cdname='Ice-vapour equivt pot temp'
      cdunits='K'

   case ('thetae_ivs','theivs')
      ivar_type=3
      ierr= RAMS_getvar('THP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%f,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,scr%f,scr%c)
      call RAMS_comp_press(nx,ny,nz,scr%c)
      call RAMS_comp_mults(nx,ny,nz,scr%c,100.) ! Pressure is in Pa.


      !----- Get ice and liquid mixing ratio. ---------------------------------------------!
      call RAMS_comp_zero(nx,ny,nz,scr%l)
      call RAMS_comp_zero(nx,ny,nz,scr%i)
      ierr= RAMS_getvar('RCP',idim_type,ngrd,scr%x,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%l,scr%x) ! Cloud (liquid)
      ierr= RAMS_getvar('RRP',idim_type,ngrd,scr%x,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%l,scr%x) ! Rain (liquid)
      ierr= RAMS_getvar('RPP',idim_type,ngrd,scr%x,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%i,scr%x) ! Pristine ice (ice)
      ierr= RAMS_getvar('RSP',idim_type,ngrd,scr%x,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%i,scr%x) ! Snow (ice)
      ierr= RAMS_getvar('RAP',idim_type,ngrd,scr%x,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%i,scr%x) ! Aggregates (ice)
      ierr= RAMS_getvar('RGP',idim_type,ngrd,scr%x,b,flnm)
      if(ierr == 0) then
         ierr= RAMS_getvar('Q6',idim_type,ngrd,scr%q,b,flnm)
         if(ierr == 0) then
            call RAMS_comp_iceliq(nx,ny,nz,scr%q,scr%x,scr%y,scr%z)
            call RAMS_comp_accum (nx,ny,nz,scr%i,scr%y)
            call RAMS_comp_accum (nx,ny,nz,scr%l,scr%z)
         end if
      end if
      ierr= RAMS_getvar('RHP',idim_type,ngrd,scr%x,b,flnm)
      if(ierr == 0) then
         ierr= RAMS_getvar('Q7',idim_type,ngrd,scr%q,b,flnm)
         if(ierr == 0) then
            call RAMS_comp_iceliq(nx,ny,nz,scr%q,scr%x,scr%y,scr%z)
            call RAMS_comp_accum (nx,ny,nz,scr%i,scr%y)
            call RAMS_comp_accum (nx,ny,nz,scr%l,scr%z)
         end if
      end if


      call RAMS_comp_thetaeivs(nx,ny,nz,a,scr%f,scr%c,scr%l,scr%i)
      cdname='Sat Ice-vapour equivt pot temp'
      cdunits='K'

   case ('theta_v')
      ivar_type=3
      ierr= RAMS_getvar('THETA',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('RV',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('RTP',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_thetv(nx,ny,nz,a,scr%c,scr%d)
      cdname='virtual pot temp'
      cdunits='K'

   case ('rv')
      ivar_type=3
      ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      if(ierr == 0) then
         call RAMS_comp_mults(nx,ny,nz,a,1.e3)
         call RAMS_comp_noneg(nx,ny,nz,a)
      endif
      cdname='vapour mixing ratio'
      cdunits='g/kg'

   case ('cloud','rain','pristine','snow','aggregates','aggr','graupel','hail')
      ivar_type=3
      select case (trim(cvar))
      case ('cloud')
         ierr= RAMS_getvar('RCP',idim_type,ngrd,a,b,flnm)
         ierr_getvar = ierr_getvar + ierr
         cdname='cloud mixing ratio'
         cdunits='g/kg'
      case ('rain')
         ierr= RAMS_getvar('RRP',idim_type,ngrd,a,b,flnm)
         cdname='rain mixing ratio'
         cdunits='g/kg'
      case ('pristine')
         ierr= RAMS_getvar('RPP',idim_type,ngrd,a,b,flnm)
         cdname='pristine ice mixing ratio'
         cdunits='g/kg'

      case ('snow')
         ierr= RAMS_getvar('RSP',idim_type,ngrd,a,b,flnm)
         cdname='snow mixing ratio'
         cdunits='g/kg'

      case ('aggregates','agg')
         ierr= RAMS_getvar('RAP',idim_type,ngrd,a,b,flnm)
         cdname='aggregate mixing ratio'
         cdunits='g/kg'

      case ('graupel')
         ierr= RAMS_getvar('RGP',idim_type,ngrd,a,b,flnm)
         cdname='graupel mixing ratio'
         cdunits='g/kg'

      case ('hail')
         ierr= RAMS_getvar('RHP',idim_type,ngrd,a,b,flnm)
         cdname='hail mixing ratio'
         cdunits='g/kg'
      end select


      if(ierr == 0) then
         call RAMS_comp_mults(nx,ny,nz,a,1.e3)
         call RAMS_comp_noneg(nx,ny,nz,a)
      endif

   case ('liquid')
      ivar_type=3
      call RAMS_comp_zero(nx,ny,nz,a)
      ierr= RAMS_getvar('RCP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RRP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)

      ierr= RAMS_getvar('RGP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) then
         ierr= RAMS_getvar('Q6',idim_type,ngrd,scr%d,b,flnm)
         if(ierr == 0) then
            call RAMS_comp_fracliq(nx,ny,nz,scr%d)
            call RAMS_comp_mult(nx,ny,nz,scr%c,scr%d)
         endif
         call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      endif

      ierr= RAMS_getvar('RHP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) then
         ierr= RAMS_getvar('Q7',idim_type,ngrd,scr%d,b,flnm)
         if(ierr == 0) then
            call RAMS_comp_fracliq(nx,ny,nz,scr%d)
            call RAMS_comp_mult(nx,ny,nz,scr%c,scr%d)
         endif
         call RAMS_comp_accum(nx,ny,nz,a,scr%c)
       endif

      call RAMS_comp_mults(nx,ny,nz,a,1.e3)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='liquid mixing ratio'
      cdunits='g/kg'

   case ('ice')
      ivar_type=3
      call RAMS_comp_zero(nx,ny,nz,a)
      ierr= RAMS_getvar('RPP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RSP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      ierr= RAMS_getvar('RAP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,a,scr%c)

      ierr= RAMS_getvar('RGP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) then
         ierr= RAMS_getvar('Q6',idim_type,ngrd,scr%d,b,flnm)
         if(ierr == 0) then
            call RAMS_comp_fracice(nx,ny,nz,scr%d)
            call RAMS_comp_mult(nx,ny,nz,scr%c,scr%d)
         endif
         call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      endif

      ierr= RAMS_getvar('RHP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) then
         ierr= RAMS_getvar('Q7',idim_type,ngrd,scr%d,b,flnm)
         if(ierr == 0) then
            call RAMS_comp_fracice(nx,ny,nz,scr%d)
            call RAMS_comp_mult(nx,ny,nz,scr%c,scr%d)
         endif
         call RAMS_comp_accum(nx,ny,nz,a,scr%c)
      endif

      call RAMS_comp_mults(nx,ny,nz,a,1.e3)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='ice mixing ratio'
      cdunits='g/kg'

   case ('rtp')
      ivar_type=3
      ierr= RAMS_getvar('RTP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,1.e3)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='total water mixing ratio'
      cdunits='g/kg'


   case ('fog','low_cloud','mid_cloud','high_cloud')
      ivar_type=2

      ierr= RAMS_getvar('PBLHGT',idim_type,ngrd,scr%z,b,flnm) ! z is PBL height
      ierr_getvar = ierr_getvar + ierr

      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%h,b,flnm)   ! h is topography
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_zero(nx,ny,nz,scr%e)                      ! e is total condensed
      ierr= RAMS_getvar('RCP',idim_type,ngrd,scr%x,b,flnm)     ! x is cloud
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%e,scr%x) !
      ierr= RAMS_getvar('RRP',idim_type,ngrd,scr%x,b,flnm)     ! x is rain
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%e,scr%x) !
      ierr= RAMS_getvar('RPP',idim_type,ngrd,scr%x,b,flnm)     ! x is pristine ice
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%e,scr%x) !
      ierr= RAMS_getvar('RSP',idim_type,ngrd,scr%x,b,flnm)     ! x is snow
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%e,scr%x) !
      ierr= RAMS_getvar('RAP',idim_type,ngrd,scr%x,b,flnm)     ! x is aggregates
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%e,scr%x) !
      ierr= RAMS_getvar('RGP',idim_type,ngrd,scr%x,b,flnm)     ! x is graupel
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%e,scr%x) !
      ierr= RAMS_getvar('RHP',idim_type,ngrd,scr%x,b,flnm)     ! x is hail
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz,scr%e,scr%x) !

      call RAMS_comp_dn0(nx,ny,nz,scr%x,scr%y,scr%d,scr%h,ngrd)    ! d is density


      call RAMS_comp_zero(nx,ny,1,scr%u) ! u is bottom
      call RAMS_comp_zero(nx,ny,1,scr%t) ! t is Top

      select case(trim(cvar))
      case ('fog')
         call RAMS_comp_clone(nx,ny,1,scr%t,scr%z)
         cdname = 'Fog'

      case ('low_cloud')
         call RAMS_comp_clone(nx,ny,1,scr%u,scr%z)
         call RAMS_comp_adds(nx,ny,1,scr%t,3500.)
         cdname = 'Low clouds'
      case ('mid_cloud')
         call RAMS_comp_adds(nx,ny,1,scr%u,3500.)
         call RAMS_comp_adds(nx,ny,1,scr%t,7000.)
         cdname = 'Middle clouds'
      case ('high_cloud')
         call RAMS_comp_adds(nx,ny,1,scr%u,7000.)
         call RAMS_comp_adds(nx,ny,1,scr%t,100000.)
         cdname = 'High clouds'
      end select


      call RAMS_comp_massint(nx,ny,nz,1,a,scr%d,scr%e,scr%u,scr%t,scr%h,ngrd)
      cdunits = 'kg/m2'


   case ('dewptk','dewptc')
      ivar_type=3
      ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dewK(nx,ny,nz,a,scr%c,scr%d)
      select case (trim(cvar))
      case ('dewptk')
         cdname='dewpoint temp'
         cdunits='K'
      case ('dewptc')
         call RAMS_comp_tempC(nx,ny,nz,1,a)
         cdname='dewpoint temp'
         cdunits='C'
      end select

   case ('rh')
      ivar_type=3
      ierr= RAMS_getvar('RV',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_rh(nx,ny,nz,a,scr%c,scr%d)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='relative humidity'
      cdunits='pct'



   case ('cloud_concen_mg')
      ivar_type=3
   ! variable 18 is iccp
      ierr= RAMS_getvar('CCP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='cloud concen'
      cdunits='#/mg'

   case ('rain_concen_kg','pris_concen_kg','snow_concen_kg','agg_concen_kg'                &
        ,'graup_concen_kg','hail_concen_kg')
      ivar_type=3
      select case (trim(cvar))
      case ('rain_concen_kg')
         ierr= RAMS_getvar('CRP',idim_type,ngrd,a,b,flnm)
         cdname='rain concen'
         cdunits='#/kg'
      case ('pris_concen_kg')
         ierr= RAMS_getvar('CPP',idim_type,ngrd,a,b,flnm)
         cdname='pristine concen'
         cdunits='#/kg'
      case ('snow_concen_kg')
         ierr= RAMS_getvar('CSP',idim_type,ngrd,a,b,flnm)
         cdname='snow concen'
         cdunits='#/kg'
      case ('agg_concen_kg')
         ierr= RAMS_getvar('CAP',idim_type,ngrd,a,b,flnm)
         cdname='aggregate concen'
         cdunits='#/kg'
      case ('graup_concen_kg')
         ierr= RAMS_getvar('CGP',idim_type,ngrd,a,b,flnm)
         cdname='graupel concen'
         cdunits='#/kg'
      case ('hail_concen_kg')
         ierr= RAMS_getvar('CHP',idim_type,ngrd,a,b,flnm)
         cdname='hail concen'
         cdunits='#/kg'
      end select
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)

   case ('cloud_concen_cm3')
      ivar_type=3
      ierr= RAMS_getvar('CCP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,scr%z,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,nz,a,scr%d)
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='cloud concen'
      cdunits='#/cm3'

   case ('rain_concen_m3','pris_concen_m3','snow_concen_m3','agg_concen_m3'                &
        ,'graup_concen_m3','hail_concen_m3')
      ivar_type=3
      select case (trim(cvar))
      case ('rain_concen_m3')
         ierr= RAMS_getvar('CRP',idim_type,ngrd,a,b,flnm)
         cdname='rain concen'
         cdunits='#/m3'

      case ('pris_concen_m3')
         ierr= RAMS_getvar('CPP',idim_type,ngrd,a,b,flnm)
         cdname='pristine concen'
         cdunits='#/m3'

      case ('snow_concen_m3')
         ierr= RAMS_getvar('CSP',idim_type,ngrd,a,b,flnm)
         cdname='snow concen'
         cdunits='#/m3'

      case ('agg_concen_m3')
         ierr= RAMS_getvar('CAP',idim_type,ngrd,a,b,flnm)
         cdname='aggregates concen'
         cdunits='#/m3'

      case ('graup_concen_m3')
         ierr= RAMS_getvar('CGP',idim_type,ngrd,a,b,flnm)
         cdname='graupel concen'
         cdunits='#/m3'

      case ('hail_concen_m3')
         ierr= RAMS_getvar('CHP',idim_type,ngrd,a,b,flnm)
         cdname='hail concen'
         cdunits='#/m3'
      end select

      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,scr%z,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,nz,a,scr%d)
      call RAMS_comp_noneg(nx,ny,nz,a)

   case ('ccn_concen')
      ivar_type=3
      ierr= RAMS_getvar('CCCNP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      call RAMS_comp_mults(nx,ny,nz,a,1.e-6)
      cdname='ccnx concen'
      cdunits='#/mg'

   case ('ifn_conc')
      ivar_type=3
      ierr= RAMS_getvar('CIFNP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='CN mixing ratio'
      cdunits='#/kg'

   case ('cloud_diam','rain_diam','pristine_diam','snow_diam','agg_diam'                   &
        ,'graup_diam','hail_diam')
      ivar_type=3

      select case (trim(cvar))
      case ('cloud_diam')
         ihyd    = 1
         fmult   = 1.e6
         cdname  = 'cloud diam'
         cdunits = 'microns'
         ierr= RAMS_getvar('RCP',idim_type,ngrd,a,b,flnm)
         ierr_getvar = ierr_getvar + ierr
         ierr= RAMS_getvar('CCP',idim_type,ngrd,scr%c,b,flnm)
         ierr_getvar = ierr_getvar + ierr
      case ('rain_diam')
         ihyd    = 2
         fmult   = 1.e3
         cdname  = 'rain diam'
         cdunits = 'mm'
         ierr= RAMS_getvar('RRP',idim_type,ngrd,a,b,flnm)
         ierr_getvar = ierr_getvar + ierr
         ierr= RAMS_getvar('CRP',idim_type,ngrd,scr%c,b,flnm)
         ierr_getvar = ierr_getvar + ierr

      case ('pris_diam')
         ihyd    = 3
         fmult   = 1.e6
         cdname  = 'pristine diam'
         cdunits = 'microns'
         ierr= RAMS_getvar('RPP',idim_type,ngrd,a,b,flnm)
         ierr_getvar = ierr_getvar + ierr
         ierr= RAMS_getvar('CPP',idim_type,ngrd,scr%c,b,flnm)
         ierr_getvar = ierr_getvar + ierr

      case ('snow_diam')
         ihyd    = 4
         fmult   = 1.e3
         cdname  = 'snow diam'
         cdunits = 'mm'
         ierr= RAMS_getvar('RSP',idim_type,ngrd,a,b,flnm)
         ierr_getvar = ierr_getvar + ierr
         ierr= RAMS_getvar('CSP',idim_type,ngrd,scr%c,b,flnm)
         ierr_getvar = ierr_getvar + ierr

      case ('agg_diam')
         ihyd    = 5
         fmult   = 1.e3
         cdname  = 'aggregates diam'
         cdunits = 'mm'
         ierr= RAMS_getvar('RAP',idim_type,ngrd,a,b,flnm)
         ierr_getvar = ierr_getvar + ierr
         ierr= RAMS_getvar('CAP',idim_type,ngrd,scr%c,b,flnm)
         ierr_getvar = ierr_getvar + ierr

      case ('graup_diam')
         ihyd    = 6
         fmult   = 1.e3
         cdname  = 'graupel diam'
         cdunits = 'mm'
         ierr= RAMS_getvar('RGP',idim_type,ngrd,a,b,flnm)
         ierr_getvar = ierr_getvar + ierr
         ierr= RAMS_getvar('CGP',idim_type,ngrd,scr%c,b,flnm)
         ierr_getvar = ierr_getvar + ierr

      case ('hail_diam')
         ihyd    = 7
         fmult   = 1.e3
         cdname  = 'hail diam'
         cdunits = 'mm'
         ierr= RAMS_getvar('RHP',idim_type,ngrd,a,b,flnm)
         ierr_getvar = ierr_getvar + ierr
         ierr= RAMS_getvar('CHP',idim_type,ngrd,scr%c,b,flnm)
         ierr_getvar = ierr_getvar + ierr
      end select

      call RAMS_comp_hydrodiam(nx,ny,nz,a,scr%c,cfmas(ihyd),pwmas(ihyd))
      call RAMS_comp_mults(nx,ny,nz,a,fmult)
      call RAMS_comp_noneg(nx,ny,nz,a)

   case ('rain_temp')
      ivar_type=3
      ierr= RAMS_getvar('Q2',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_raintemp(nx,ny,nz,a)
      cdname='rain temperature'
      cdunits='K'

   case ('graup_temp')
      ivar_type=3
      ierr= RAMS_getvar('Q6',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_qtcpcp(nx,ny,nz,a)
      cdname='graupel temperature'
      cdunits='C'

   case ('hail_temp')
      ivar_type=3
      ierr= RAMS_getvar('Q7',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_qtcpcp(nx,ny,nz,a)
      cdname='hail temperature'
      cdunits='C'

   case ('graup_fracliq')
      ivar_type=3
      ierr= RAMS_getvar('Q6',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_fracliq(nx,ny,nz,a)
      cdname='graupel liq frac'
      cdunits=' '

   case ('hail_fracliq')
      ivar_type=3
      ierr= RAMS_getvar('Q7',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_fracliq(nx,ny,nz,a)
      cdname='hail liq frac'
      cdunits=' '

   case ('geo')
      ivar_type=3
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_z(nx,ny,nz,a,scr%c,ngrd)
      cdname='geopotential height'
      cdunits='m'

   case ('tke')
      ivar_type=3
      ierr= RAMS_getvar('TKEP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='turb kinetic energy'
      cdunits='m2/s2'

   case ('thsrc') 
      ivar_type=6
      ierr= RAMS_getvar('THSRC',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz*ncld,a,86400.)
      cdname='deep conv heat rate'
      cdunits='K/day'

   case ('rtsrc') 
      ivar_type=6
      ierr= RAMS_getvar('RTSRC',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz*ncld,a,86400.)
      call RAMS_comp_mults(nx,ny,nz*ncld,a,1000.)
      cdname='deep conv moist rate'
      cdunits='g/kg/day'

   case ('co2src') 
      ivar_type=6
      if (co2_on /= 0) then
         ierr= RAMS_getvar('CO2SRC',idim_type,ngrd,a,b,flnm)
         ierr_getvar = ierr_getvar + ierr
         call RAMS_comp_mults(nx,ny,nz*ncld,a,86400.)
      else
         write (unit=*,fmt='(a,1x,es12.5)') '       # Assigning zero CO2SRC = ',0.
         call ae0(nx*ny*nz*ncld,a,0.)
      end if
      cdname='deep conv co2 rate'
      cdunits='umol/mol/day'

   case ('fthrd') 
      ivar_type=3
      ierr= RAMS_getvar('FTHRD',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,86400.)
      cdname='rad heat rate'
      cdunits='K/day'

   case ('fthrd_lw') 
      ivar_type=3
      ierr= RAMS_getvar('FTHRD_LW',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,nz,a,86400.)
      cdname='rad heat rate (longwave)'
      cdunits='K/day'

   case ('khh') 
      ivar_type=3
      ierr= RAMS_getvar('HKH',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='horiz diffusion coeff'
      cdunits='m2/s'

   case ('khv') 
      ivar_type=3
      ierr= RAMS_getvar('VKH',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='vert diffusion coeff'
      cdunits='m2/s'

   case ('totpcp','totpcp_in','precip','precip_in')
      ivar_type=2
      call RAMS_comp_zero(nx,ny,1,a)
      ierr= RAMS_getvar('ACCPR',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('ACCPP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('ACCPS',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('ACCPA',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('ACCPG',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('ACCPH',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)

      select case (trim(cvar))
      case ('precip','precip_in')
         ierr= RAMS_getvar('ACONPR',idim_type,ngrd,scr%c,b,flnm)
         if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
         cdname='total accum precip'
      case ('totpcp','totpcp_in')
         cdname='total resolved precip'
      end select

      select case (trim(cvar))
      case ('totpcp','precip')
         cdunits='kg/m2'
      case ('precip_in','totpcp_in')
         call RAMS_comp_mults(nx,ny,nz,a,.03937)
         cdunits='in liq'
      end select
      call RAMS_comp_noneg(nx,ny,1,a)

   case ('icepcp')
      ivar_type=2
      call RAMS_comp_zero(nx,ny,1,a)
      ierr= RAMS_getvar('ACCPP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('ACCPS',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('ACCPA',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)

      cdname='purely ice precip'
      cdunits='kg/m2'
      call RAMS_comp_noneg(nx,ny,1,a)

   case ('mixpcp')
      ivar_type=2
      call RAMS_comp_zero(nx,ny,1,a)
      ierr= RAMS_getvar('ACCPG',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('ACCPH',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)

      cdname='mixed (ice/liq) precip'
      cdunits='kg/m2'
      call RAMS_comp_noneg(nx,ny,1,a)

   case ('pcprate','pcprate_in')
      ivar_type=2
      call RAMS_comp_zero(nx,ny,1,a)
      ierr= RAMS_getvar('PCPRR',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('PCPRP',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('PCPRS',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('PCPRA',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('PCPRG',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      ierr= RAMS_getvar('PCPRH',idim_type,ngrd,scr%c,b,flnm)
      if(ierr == 0) call RAMS_comp_accum(nx,ny,1,a,scr%c)
      call RAMS_comp_noneg(nx,ny,1,a)
      cdname='resolved precip rate'

      select case (trim(cvar))
      case ('pcprate','precipr')
         call RAMS_comp_mults(nx,ny,1,a,3600.)
         cdunits='mm/hr'
      case ('pcprate_in','precipr_in')
         call RAMS_comp_mults(nx,ny,1,a,141.732)
         cdunits='in/hr'
      end select

   case ('conpcp','conprr')
      ivar_type=9
      ierr= RAMS_getvar('CONPRR',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,ncld,a,3600.)
      call RAMS_comp_noneg(nx,ny,ncld,a)
      cdname='convective pcp rate'
      cdunits='mm/hr'

   case ('acccon')
      ivar_type=2
      ierr= RAMS_getvar('ACONPR',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,1,a)
      cdname='accum convective pcp'
      cdunits='mm'

   case ('cape','cine','showalter')
      ivar_type=2

      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%i,b,flnm)  ! i is topt
      ierr_getvar = ierr_getvar + ierr

      ierr= RAMS_getvar('THP',idim_type,ngrd,scr%h,b,flnm)   ! h is thil
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%c,b,flnm)    ! c is Exner
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%f,b,flnm) ! f is theta
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_tempK(nx,ny,nz,scr%f,scr%c)             ! f is temperature
      call RAMS_comp_clone(nx,ny,nz,scr%g,scr%c)             ! g is Exner
      call RAMS_comp_press(nx,ny,nz,scr%c)                   ! c is pressure
      call RAMS_comp_mults(nx,ny,nz,scr%c,100.)              ! c is pressure
      ierr= RAMS_getvar('RV',idim_type,ngrd,scr%d,b,flnm)    ! d is rv
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('RTP',idim_type,ngrd,scr%e,b,flnm)   ! e is rtp
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_thetaeiv(nx,ny,nz,scr%h,scr%f,scr%c,scr%d,scr%e) ! h is thetae_iv

      select case(trim(cvar))
      case ('cine')
         call RAMS_comp_cine_cape(nx,ny,nz,2,trim(cvar),a,scr%c,scr%g,scr%h,scr%e,scr%d    &
                                 ,scr%f,scr%i,ngrd)
         cdname  = 'Convective inhibition'
         cdunits = 'J/kg'
      case ('cape')
         call RAMS_comp_cine_cape(nx,ny,nz,2,trim(cvar),a,scr%c,scr%g,scr%h,scr%e,scr%d    &
                                 ,scr%f,scr%i,ngrd)
         cdname  = 'Convective Available Potential Energy'
         cdunits = 'J/kg'
      case ('showalter')
         call RAMS_comp_showalter(nx,ny,nz,a,scr%c,scr%g,scr%h,scr%e,scr%d,scr%f,scr%i,ngrd)
      end select
      cdname='Showalter index'
      cdunits='K'

   case ('sweat')
      ivar_type=2

      ierr= RAMS_getvar('RV',idim_type,ngrd,scr%r,b,flnm)     ! r is rv
      ierr_getvar = ierr_getvar + ierr

      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%t,b,flnm)  ! t is theta
      ierr_getvar = ierr_getvar + ierr

      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%e,b,flnm)     ! e is Exner
      ierr_getvar = ierr_getvar + ierr

      ierr        = RAMS_getvar('UP',idim_type,ngrd,scr%u,b,flnm) ! u is U
      ierr_getvar = ierr_getvar + ierr
      ierr        = RAMS_getvar('VP',idim_type,ngrd,scr%v,b,flnm) ! v is V
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_rotate(nx,ny,nz,scr%u,scr%v,ngrd)
      call RAMS_comp_avgu(nx,ny,nz,scr%u)
      call RAMS_comp_avgv(nx,ny,nz,scr%v)


      call RAMS_comp_clone(nx,ny,nz,scr%d,scr%r)              ! d is rv
      call RAMS_comp_dewK(nx,ny,nz,scr%d,scr%e,scr%t)         ! d is dewpoint theta
      call RAMS_comp_tempK(nx,ny,nz,scr%t,scr%e)              ! t is temperature

      call RAMS_comp_clone(nx,ny,nz,scr%p,scr%e)              ! p is pressure
      call RAMS_comp_press(nx,ny,nz,scr%p)                    ! p is pressure
      call RAMS_comp_mults(nx,ny,nz,scr%p,100.)               ! p is pressure

      call RAMS_comp_sweat(nx,ny,nz,a,scr%p,scr%t,scr%d,scr%u,scr%v,scr%o)
      cdname='SWEAT index'
      cdunits='idx'

   case ('kindex')
      ivar_type=2

      ierr= RAMS_getvar('RV',idim_type,ngrd,scr%r,b,flnm)     ! r is rv
      ierr_getvar = ierr_getvar + ierr

      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%t,b,flnm)  ! t is theta
      ierr_getvar = ierr_getvar + ierr

      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%e,b,flnm)     ! e is Exner
      ierr_getvar = ierr_getvar + ierr


      call RAMS_comp_clone(nx,ny,nz,scr%d,scr%r)              ! d is rv
      call RAMS_comp_dewK(nx,ny,nz,scr%d,scr%e,scr%t)         ! d is dewpoint theta
      call RAMS_comp_tempK(nx,ny,nz,scr%t,scr%e)              ! t is temperature

      call RAMS_comp_clone(nx,ny,nz,scr%p,scr%e)              ! p is pressure
      call RAMS_comp_press(nx,ny,nz,scr%p)                    ! p is pressure
      call RAMS_comp_mults(nx,ny,nz,scr%p,100.)               ! p is pressure

      call RAMS_comp_kindex(nx,ny,nz,a,scr%p,scr%t,scr%d)
      cdname='K index'
      cdunits='idx'

   case ('cfxup')
      ivar_type=6
      ierr= RAMS_getvar('CFXUP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz*ncld,a)
      cdname='Convective upward flux'
      cdunits='kg/m2/s'

   case ('cfxdn')
      ivar_type=6
      ierr= RAMS_getvar('CFXDN',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_nopos(nx,ny,nz*ncld,a)
      cdname='Convective downward flux'
      cdunits='kg/m2/s'

   case ('dfxup')
      ivar_type=6
      ierr= RAMS_getvar('DFXUP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz*ncld,a)
      cdname='Detrainment upward flux'
      cdunits='kg/m2/s'

   case ('dfxdn')
      ivar_type=6
      ierr= RAMS_getvar('DFXDN',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz*ncld,a)
      cdname='Detrainment upward flux'
      cdunits='kg/m2/s'

   case ('efxup')
      ivar_type=6
      ierr= RAMS_getvar('EFXUP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz*ncld,a)
      cdname='Entrainment upward flux'
      cdunits='kg/m2/s'

   case ('efxdn')
      ivar_type=6
      ierr= RAMS_getvar('EFXDN',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz*ncld,a)
      cdname='Entrainment upward flux'
      cdunits='kg/m2/s'

   ! Extra turbulence parameters

   case ('tkem')
      ivar_type=3
      ierr= RAMS_getvar('TKEPB',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   !   call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='Mean turbulent kinetic energy'
      cdunits='m2/s2'

   case ('ltscale')
      ivar_type=3
      ierr= RAMS_getvar('TL',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   !   call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='Mean Lagrangean time scale'
      cdunits='s'

   case ('sigw')
      ivar_type=3
      ierr= RAMS_getvar('SIGW',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
   !   call RAMS_comp_noneg(nx,ny,nz,a)
      cdname='Mean vertical velocity standard deviation'
      cdunits='m/s'

   case ('snowdepth')
       ivar_type=2
       ierr= RAMS_getvar('SNOW_DEPTH',idim_type,ngrd,a,b,flnm)
       ierr_getvar = ierr_getvar + ierr
       cdname='Depth of the snow layer'
       cdunits='m'

   case ('pblhgt')
      ivar_type=2
      ierr= RAMS_getvar('PBLHGT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='PBL height'
      cdunits='m'

   case ('lmo')
       ivar_type=2
       ierr= RAMS_getvar('LMO',idim_type,ngrd,a,b,flnm)
       ierr_getvar = ierr_getvar + ierr
       cdname='Obukhov lenght scale'
       cdunits='m'

   case ('umom_flx')
      ivar_type=2
      ierr= RAMS_getvar('UW',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,scr%z,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,1,a,scr%d)
      cdname='sfc u-momentum flx'
      cdunits='Pa'

   case ('vmom_flx')
      ivar_type=2
      ierr= RAMS_getvar('VW',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,scr%z,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,1,a,scr%d)
      cdname='sfc v-momentum flx'
      cdunits='Pa'

   case ('wmom_flx')
      ivar_type=2
      ierr= RAMS_getvar('SFLUX_W',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_dn0(nx,ny,nz,scr%z,scr%c,scr%d,scr%e,ngrd)
      call RAMS_comp_mult(nx,ny,1,a,scr%d)
      cdname='sfc w-momentum flx'
      cdunits='Pa'

   case ('cosz')
      ivar_type=2
      ierr= RAMS_getvar('COSZ',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='co-sine of zenith angle'
      cdunits=' '

   case ('zen')
      ivar_type=2
      ierr= RAMS_getvar('COSZ',idim_type,ngrd,scr%c,b,flnm)
      call RAMS_comp_zenith(nx,ny,scr%c,a)
      ierr_getvar = ierr_getvar + ierr
      cdname='zenith angle'
      cdunits='deg'

   case ('par')
      ivar_type=2
      ierr= RAMS_getvar('PAR_BEAM',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PAR_DIFFUSE',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_accum(nx,ny,1,a,scr%e)
      call RAMS_comp_noneg(nx,ny,1,a)

      cdname='sfc PAR rad.'
      cdunits='W/m2'

   case ('nir')
      ivar_type=2
      ierr= RAMS_getvar('NIR_BEAM',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('NIR_DIFFUSE',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_accum(nx,ny,1,a,scr%e)
      call RAMS_comp_noneg(nx,ny,1,a)

      cdname='sfc NIR rad.'
      cdunits='W/m2'

   case ('rshort')
      ivar_type=2
      ierr= RAMS_getvar('RSHORT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='sfc SW. rad.'
      cdunits='W/m2'

   case ('par_beam')
      ivar_type=2
      ierr= RAMS_getvar('PAR_BEAM',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,1,a)

      cdname='direct sfc PAR rad.'
      cdunits='W/m2'

   case ('nir_beam')
      ivar_type=2
      ierr= RAMS_getvar('NIR_BEAM',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,1,a)

      cdname='direct sfc NIR rad.'
      cdunits='W/m2'

   case ('rshort_beam')
      ivar_type=2
      ierr= RAMS_getvar('RSHORT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('RSHORT_DIFFUSE',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_subt(nx,ny,1,a,scr%e)
      call RAMS_comp_noneg(nx,ny,1,a)
      cdname='direct sfc. SW rad.'
      cdunits='W/m2'

   case ('par_diff')
      ivar_type=2
      ierr= RAMS_getvar('PAR_DIFFUSE',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,1,a)

      cdname='diffuse sfc PAR rad.'
      cdunits='W/m2'

   case ('nir_diff')
      ivar_type=2
      ierr= RAMS_getvar('NIR_DIFFUSE',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,1,a)

      cdname='diffuse sfc NIR rad.'
      cdunits='W/m2'

   case ('rshort_diff')
      ivar_type=2
      ierr= RAMS_getvar('RSHORT_DIFFUSE',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,1,a)
      cdname='diffuse sfc. SW rad.'
      cdunits='W/m2'

   case ('rshorttoa')
      ivar_type=2
      ierr= RAMS_getvar('RSHORT_TOP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='TOA SW rad.'
      cdunits='W/m2'

   case ('rlong')
      ivar_type=2
      ierr= RAMS_getvar('RLONG',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='sfc LW. rad.'
      cdunits='W/m2'

   case ('rlongup')
      ivar_type=2
      ierr= RAMS_getvar('RLONGUP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='rlongup'
      cdunits='W/m2'

   case ('rshortuptoa')
      ivar_type=2
      ierr= RAMS_getvar('RSHORTUP_TOP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='rlongup at TOA'
      cdunits='W/m2'

   case ('rlonguptoa')
      ivar_type=2
      ierr= RAMS_getvar('RLONGUP_TOP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='rlongup at TOA'
      cdunits='W/m2'

   case ('albedt')
      ivar_type=2
      ierr= RAMS_getvar('ALBEDT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='albedt'
      cdunits=' '

   ! 2D TOPOGRAPHY AND GEOGRAPHIC VALUES

   case ('topoa')
      ivar_type=2
      ierr= RAMS_getvar('TOPTA',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='topo'
      cdunits='m'

   case ('topo')
      ivar_type=2
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='topo'
      cdunits='m'

   case ('lon','longitude')
      ivar_type=2
      ierr= RAMS_getvar('GLON',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='longitude'
      cdunits='deg'

   case ('lat','latitude')
      ivar_type=2
      ierr= RAMS_getvar('GLAT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='latitude'
      cdunits='deg'

   case ('mynum')
      ivar_type=2
      ierr= RAMS_getvar('MYNUM',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='node ID'
      cdunits=' '

   ! 2D MISCELLANEOUS FIELDS

   case ('sea_press')
      ivar_type=2
      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('PI',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('THETA',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_slpmm5(nx,ny,nz,scr%e,scr%d,scr%c,a)
      cdname='sea level pressure;'
      cdunits='mb;'

   ! Special use of sst: acquired for patch #1 even where no water exists

   case ('sst')
      ivar_type=2
      ierr = RAMS_getvar('SOIL_ENERGY',idim_type,ngrd   &
           ,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call rams_fill_sst(nx,ny,nsl*npat,nsl,a,scr%c)
      cdname='water temperature'
      cdunits='C'


   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! LEAF2 variables section

   ! If want a horiz plot, specify a string like 'tgpatch'; it will
   !   return i,j,ip array.
   ! Specify a new ivar_type, not corresponding to anal file var type.  With
   !   horiz plot, get back into iplt.  If have this var type, don't slice.
   ! Need replacement for rams3to2d because windowing is done in there.
   ! Replacement would window but not slice.
   ! Then, if want xz (vert cross section) have name like tgpatch_vert.
   ! This would return entire 4d array from hvlib.f.
   ! Then we have to slice and window with yet another replacement to rams3to2d.

   ! nkk is the record number, where n is the LEAF field number (1, 2, 3, or 4)
   ! and kk is the k level.
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


   case ('pfarea')

      ivar_type = 7
      irecind = 1
      irecsize = nx * ny * npat
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      irecind = irecind + irecsize
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='patch fractional area'
      cdunits=''

   case ('soil_z0_p','soil_z0_ps')

      irecind = 1
      irecsize = nx * ny * npat
      select case(trim(cvar))
      case ('soil_z0_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('SOIL_Z0',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('soil_z0_p')
         ivar_type = 7
      case ('soil_z0_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nx,ny,1,npat,a)
      end select

      cdname='soil roughness'
      cdunits='m'

   case ('vtype','vtype_bp')

      irecind = 1
      irecsize = nx * ny * npat

      select case (trim(cvar))
      case ('vtype_bp')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
             ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('LEAF_CLASS',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_vegclass(irecsize,a(irecind))
      select case (trim(cvar))
      case ('vtype')
         ivar_type = 7
      case ('vtype_bp')
         ivar_type = 2
         call RAMS_comp_bigpatch(nx,ny,1,npat  &
            ,a(irecind),a(1),b)
      end select
      cdname='vegetation class'
      cdunits='#'

   case ('scolour','scolour_bp')

      irecind = 1
      irecsize = nx * ny * npat

      select case (trim(cvar))
      case ('scolour_bp')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
             ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('SOIL_COLOR',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_vegclass(irecsize,a(irecind))
      select case (trim(cvar))
      case ('scolour')
         ivar_type = 7
      case ('scolour_bp')
         ivar_type = 2
         call RAMS_comp_bigpatch(nx,ny,1,npat,a(irecind),a(1),b)
      end select
      cdname='soil colour class'
      cdunits='#'

   case ('ndvi','ndvi_bp')

      irecind = 1
      irecsize = nx * ny * npat

      select case (trim(cvar))
      case ('ndvi_bp')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
             ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('VEG_NDVIC',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('ndvi')
         ivar_type = 7
      case ('ndvi_bp')
         ivar_type = 2
         call RAMS_comp_bigpatch(nx,ny,1,npat,a(irecind),a(1),b)
      end select

      cdname='ndvi'
      cdunits='#'

   case ('vegfrac','veg_fracarea_ps')

      irecind = 1
      irecsize = nx * ny * npat

      select case (trim(cvar))
      case ('veg_fracarea_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('VEG_FRACAREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('vegfrac')
         ivar_type = 7
      case ('veg_fracarea_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nx,ny,1,npat,a)
      end select

      cdname='vegetation frac area'
      cdunits=''

   case ('land')

      ivar_type = 2
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_1minus(nx,ny,1,a)
      cdname='land frac area'
      cdunits=''

   case ('agb','vegagb','agb_ps')

      irecind = 1
      irecsize = nx * ny * npat
      select case (trim(cvar))
      case ('agb_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('VEG_AGB',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('agb','vegagb')
         ivar_type = 7
      case ('agb_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='above ground biomass'
      cdunits='kgC/m2'

   case ('lai','veglai','lai_ps')

      irecind = 1
      irecsize = nx * ny * npat
      select case (trim(cvar))
      case ('lai_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('VEG_LAI',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('lai','veglai')
         ivar_type = 7
      case ('lai_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='green leaf area index'
      cdunits=''


   case ('tai','tai_ps')

      irecind = 1
      irecsize = nx * ny * npat

      select case (trim(cvar))
      case ('tai_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('VEG_TAI',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('tai')
         ivar_type = 7
      case ('tai_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname=' total leaf area index'
      cdunits=''


   case ('vegalb','vegalb_ps')

      irecind = 1
      irecsize = nx * ny * npat

      select case (trim(cvar))
      case ('tai_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('VEG_ALBEDO',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('vegalb')
         ivar_type = 7
      case ('vegalb_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='vegetation albedo'
      cdunits=''

   case ('ustar')
      ivar_type = 7
      ierr = RAMS_getvar('USTAR',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='friction velocity'
      cdunits='m/s'

   case ('tstar')
      ivar_type = 7
      ierr = RAMS_getvar('TSTAR',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='potential temperature scale'
      cdunits='K'

   case ('rstar')
      ivar_type = 7
      ierr = RAMS_getvar('RSTAR',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a,1000.)
      cdname='water vapour mixing ratio scale'
      cdunits='g/kg'

   case ('cstar')
      ivar_type = 7
      ierr = RAMS_getvar('TSTAR',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='CO2 mixing ratio scale'
      cdunits='umol/mol'

   case ('z0','z0_ps')

      irecind = 1
      irecsize = nx * ny * npat

      select case (trim(cvar))
      case ('z0_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
      end select

      irecind = irecind + irecsize
      ierr = RAMS_getvar('PATCH_ROUGH',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('z0_p')
         ivar_type = 7
      case ('z0_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='roughness'
      cdunits='m'


   case ('net_z0_p','net_z0_ps')

      irecind = 1
      irecsize = nx * ny * npat
      select case (trim(cvar))
      case ('net_z0_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
      end select

      irecind = irecind + irecsize
      ierr = RAMS_getvar('NET_Z0',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('net_z0_p')
         ivar_type = 7
      case ('net_z0_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='net roughness'
      cdunits='m'

   case ('vegz0','vegz0_ps')

      irecind = 1
      irecsize = nx * ny * npat

      select case (trim(cvar))
      case ('vegz0_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select


      ierr = RAMS_getvar('VEG_ROUGH',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('vegz0')
         ivar_type = 7
      case ('vegz0_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nx,ny,1,npat,a)
      end select

      cdname='vegetation roughness'
      cdunits='m'

   case ('vegdisp','vegdisp_ps')

      irecind = 1
      irecsize = nx * ny * npat

      select case (trim(cvar))
      case ('veg_disp_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('VEG_DISPLACE',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('vegdisp')
         ivar_type = 7
      case ('vegdisp_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nx,ny,1,npat,a)
      end select

      cdname='vegetation displacement height'
      cdunits='m'

   case ('veghgt','veghgt_ps')

      irecind = 1
      irecsize = nx * ny * npat

      select case (trim(cvar))
      case ('veg_disp_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('VEG_HEIGHT',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('veghgt')
         ivar_type = 7
      case ('veghgt_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nx,ny,1,npat,a)
      end select

      cdname='vegetation height'
      cdunits='m'

   case ('patch_wetind')

      ivar_type = 7
      irecind = 1
      irecsize = nx * ny * npat
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      irecind = irecind + irecsize
      ierr = RAMS_getvar('WET_INDEX',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='patch wetness index'
      cdunits=''

   case ('snowlevels')

      ivar_type = 7
      irecind = 1
      irecsize = nx * ny * npat
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      irecind = irecind + irecsize
      ierr = RAMS_getvar('KSNOW',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='number of snow levels'
      cdunits='#'

   case ('lwater_p','lwater_ps')

      irecind = 1
      irecsize = nx * ny * npat
      select case (trim(cvar))
      case ('lwater_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
      irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('VEG_WATER',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('lwater_p')
         ivar_type = 7
      case ('lwater_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nx,ny,1,npat,a)
      end select

       cdname='leaf water'
      cdunits='kg/m2'



   case ('rvcan','rvcan_ps')

      irecind = 1
      irecsize = nx * ny * npat
      select case (trim(cvar))
      case ('rvcan_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('CAN_RVAP',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),1.e3)

      select case (trim(cvar))
      case ('rvcan')
         ivar_type = 7
      case ('rvcan_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='canopy mixing ratio'
      cdunits='g/kg'

   case ('co2can','co2can_ps')

      irecind = 1
      irecsize = nx * ny * npat
      select case (trim(cvar))
      case ('co2can_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('CAN_CO2',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('co2can')
         ivar_type = 7
      case ('co2can_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='CO2 mixing ratio'
      cdunits='umol/mol'

   case ('gpp_p','gpp')

      irecind = 1
      irecsize = nx * ny * npat
      select case (trim(cvar))
      case ('gpp')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('GPP',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('gpp_p')
         ivar_type = 7
      case ('gpp')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Gross Primary Production'
      cdunits='umol/m2/s'

   case ('plresp_p','plresp')

      irecind = 1
      irecsize = nx * ny * npat
      select case (trim(cvar))
      case ('plresp')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('PLRESP',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('plresp_p')
         ivar_type = 7
      case ('plresp')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Plant respiration'
      cdunits='umol/m2/s'

   case ('resphet_p','resphet')

      irecind = 1
      irecsize = nx * ny * npat
      select case (trim(cvar))
      case ('resphet')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('RESPHET',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('resphet_p')
         ivar_type = 7
      case ('resphet')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Heterotrophic respiration'
      cdunits='umol/m2/s'


   case ('nee_p','nee','cflxca_p','cflxca')

      irecind = 1
      select case (trim(cvar))
      case ('nee','cflxca')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('CFLXAC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),-1.0)

      select case (trim(cvar))
      case ('nee_p','cflxca_p')
         ivar_type = 7
      case ('nee','cflxca')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Net ecosystem exchange'
      cdunits='umol/m2/s'

   case ('hflxca_p','hflxca')

      irecind = 1
      select case (trim(cvar))
      case ('hflxca')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('HFLXAC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),-1.0)

      select case (trim(cvar))
      case ('hflxca_p')
         ivar_type = 7
      case ('hflxca')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Sensible heat flux'
      cdunits='W/m2'

   case ('hflxvc_p','hflxvc')

      irecind = 1
      select case (trim(cvar))
      case ('hflxvc')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('HFLXVC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('hflxvc_p')
         ivar_type = 7
      case ('hflxvc')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Sensible heat (Veg->CAS)'
      cdunits='W/m2'

   case ('hflxgc_p','hflxgc')

      irecind = 1
      select case (trim(cvar))
      case ('hflxgc')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('HFLXGC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('hflxgc_p')
         ivar_type = 7
      case ('hflxgc')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Sensible heat (Gnd->CAS)'
      cdunits='W/m2'

   case ('qwflxca_p','qwflxca')

      irecind = 1
      select case (trim(cvar))
      case ('qwflxca')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('QWFLXAC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),-1.0)

      select case (trim(cvar))
      case ('qwflxca_p')
         ivar_type = 7
      case ('qwflxca')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Latent heat flux'
      cdunits='W/m2'

   case ('qwflxvc_p','qwflxvc')

      irecind = 1
      select case (trim(cvar))
      case ('qwflxvc')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('QWFLXVC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('qwflxvc_p')
         ivar_type = 7
      case ('qwflxvc')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Latent heat (Veg->CAS)'
      cdunits='W/m2'

   case ('qwflxgc_p','qwflxgc')

      irecind = 1
      select case (trim(cvar))
      case ('qwflxgc')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('QWFLXGC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('qwflxgc_p')
         ivar_type = 7
      case ('qwflxgc')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Latent heat (Gnd->CAS)'
      cdunits='W/m2'

   case ('qtransp_p','qtransp')

      irecind = 1
      select case (trim(cvar))
      case ('qtransp')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('QTRANSP',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('qtransp_p')
         ivar_type = 7
      case ('qtransp')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Latent heat (Transp.)'
      cdunits='W/m2'

   case ('qrunoff_p','qrunoff')

      irecind = 1
      select case (trim(cvar))
      case ('qrunoff')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('QRUNOFF',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('qrunoff_p')
         ivar_type = 7
      case ('qrunoff')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Energy loss: sfc runoff'
      cdunits='W/m2'

   case ('qdrainage_p','qdrainage')

      irecind = 1
      select case (trim(cvar))
      case ('qdrainage')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('QDRAINAGE',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('qdrainage_p')
         ivar_type = 7
      case ('qdrainage')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Energy loss: Sub-sfc runoff'
      cdunits='W/m2'

   case ('qthroughfall_p','qthroughfall')

      irecind = 1
      select case (trim(cvar))
      case ('qthroughfall')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('QTHROUGHFALL',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('qthroughfall_p')
         ivar_type = 7
      case ('qthroughfall')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Energy input: throughfall'
      cdunits='W/m2'

   case ('qintercepted_p','qintercepted')

      irecind = 1
      select case (trim(cvar))
      case ('qintercepted')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('QINTERCEPTED',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('qintercepted_p')
         ivar_type = 7
      case ('qintercepted')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Energy input: Intercepted water'
      cdunits='W/m2'

   case ('qwshed_p','qwshed')

      irecind = 1
      select case (trim(cvar))
      case ('qwshed')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('QWSHED',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('qwshed_p')
         ivar_type = 7
      case ('qwshed')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Energy input: Water shedding'
      cdunits='W/m2'



   case ('eflxca_p','eflxca')

      irecind = 1
      select case (trim(cvar))
      case ('eflxca')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('EFLXAC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),-1.0)

      select case (trim(cvar))
      case ('eflxca_p')
         ivar_type = 7
      case ('eflxca')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Eddy flux for enthalpy'
      cdunits='W/m2'

   case ('wflxca_p','wflxca')

      irecind = 1
      select case (trim(cvar))
      case ('wflxca')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('WFLXAC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),-86400.)

      select case (trim(cvar))
      case ('wflxca_p')
         ivar_type = 7
      case ('wflxca')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Eddy flux for water'
      cdunits='kg/m2/day'

   case ('wflxvc_p','wflxvc')

      irecind = 1
      select case (trim(cvar))
      case ('wflxvc')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('WFLXVC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),86400.)

      select case (trim(cvar))
      case ('wflxvc_p')
         ivar_type = 7
      case ('wflxvc')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Veg. sfc. evaporation'
      cdunits='kg/m2/day'

   case ('wflxgc_p','wflxgc')

      irecind = 1
      select case (trim(cvar))
      case ('wflxgc')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('WFLXGC',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),86400.)

      select case (trim(cvar))
      case ('wflxgc_p')
         ivar_type = 7
      case ('wflxgc')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Ground evaporation'
      cdunits='kg/m2/day'

   case ('transp_p','transp')

      irecind = 1
      select case (trim(cvar))
      case ('transp')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('TRANSP',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),86400.)

      select case (trim(cvar))
      case ('transp_p')
         ivar_type = 7
      case ('transp')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Transpiration'
      cdunits='kg/m2/day'

   case ('runoff_p','runoff')

      irecind = 1
      select case (trim(cvar))
      case ('runoff')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('RUNOFF',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),86400.)

      select case (trim(cvar))
      case ('runoff_p')
         ivar_type = 7
      case ('runoff')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Surface runoff'
      cdunits='kg/m2/day'

   case ('drainage_p','drainage')

      irecind = 1
      select case (trim(cvar))
      case ('drainage')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('DRAINAGE',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),86400.)

      select case (trim(cvar))
      case ('drainage_p')
         ivar_type = 7
      case ('drainage')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Sub-surface runoff'
      cdunits='kg/m2/day'

   case ('throughfall_p','throughfall')

      irecind = 1
      select case (trim(cvar))
      case ('throughfall')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('THROUGHFALL',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),86400.)

      select case (trim(cvar))
      case ('throughfall_p')
         ivar_type = 7
      case ('throughfall')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Throughfall'
      cdunits='kg/m2/day'

   case ('intercepted_p','intercepted')

      irecind = 1
      select case (trim(cvar))
      case ('intercepted')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('INTERCEPTED',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),86400.)

      select case (trim(cvar))
      case ('intercepted_p')
         ivar_type = 7
      case ('intercepted')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Intercepted water'
      cdunits='kg/m2/day'

   case ('wshed_p','wshed')

      irecind = 1
      select case (trim(cvar))
      case ('wshed')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('WSHED',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),86400.)

      select case (trim(cvar))
      case ('wshed_p')
         ivar_type = 7
      case ('wshed')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Water shedding'
      cdunits='kg/m2/day'


   case ('tveg','tveg_ps')

      irecind = 1
      irecsize = nx * ny * npat
      select case (trim(cvar))
      case ('tveg_ps')
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)

         irecind = irecind + irecsize
      end select
      
      ierr = RAMS_getvar('CAN_THETA',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('CAN_PRSS',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_theta2temp(nx,ny,npat,scr%e,scr%d)
      call RAMS_comp_tempC(nx,ny,1,npat,scr%e)
      
      ierr = RAMS_getvar('VEG_ENERGY',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('VEG_WATER',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('VEG_HCAP',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_tvegc(nx,ny,npat,a(irecind),scr%c,scr%d,scr%e)

      !----- Filling first patch with SST. -------------------------------!
      ierr = RAMS_getvar('SOIL_ENERGY',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call rams_fill_sst(nx,ny,nsl*npat,nsl,a(irecind),scr%e)

      select case (trim(cvar))
      case ('tveg')
         ivar_type = 7
      case ('tveg_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='vegetation temperature'
      cdunits='C'

   case ('tcan','tcan_ps')

      irecind = 1
      select case (trim(cvar))
      case ('tcan_ps')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('CAN_THETA',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('CAN_PRSS',idim_type,ngrd,scr%c,b,flnm)
      call RAMS_comp_theta2temp(nx,ny,npat,a(irecind),scr%c)
      call RAMS_comp_tempC(nx,ny,1,npat,a(irecind))

      select case (trim(cvar))
      case ('tcan')
         ivar_type = 7
      case ('tcan_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='canopy temperature'
      cdunits='C'


   case ('pcan','pcan_ps')

      irecind = 1
      select case (trim(cvar))
      case ('pcan_ps')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('CAN_PRSS',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_mults(nx,ny,npat,a(irecind),1.e-2)

      select case (trim(cvar))
      case ('pcan')
         ivar_type = 7
      case ('pcan_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='canopy press'
      cdunits='hPa'

   case ('thcan','thcan_ps')

      irecind = 1
      select case (trim(cvar))
      case ('thcan_ps')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('CAN_THETA',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('thcan')
         ivar_type = 7
      case ('thcan_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='canopy potential temperature'
      cdunits='K'

   case ('thecan','thecan_ps')

      irecind = 1
      select case (trim(cvar))
      case ('thecan_ps')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('CAN_THEIV',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('thecan')
         ivar_type = 7
      case ('thecan_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='canopy equiv. pot. temperature'
      cdunits='K'

   case ('rshort_gnd','rshort_gnd_ps')

      irecind = 1
      select case (trim(cvar))
      case ('rshort_gnd_ps')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('RSHORT_GND',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('rshort_gnd')
         ivar_type = 7
      case ('rshort_gnd_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Ground shortwave radiation'
      cdunits='W/m2'

   case ('rlong_gnd','rlong_gnd_ps')

      irecind = 1
      select case (trim(cvar))
      case ('rlong_gnd_ps')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('RLONG_GND',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('rlong_gnd')
         ivar_type = 7
      case ('rlong_gnd_ps')
         ivar_type = 2
         call RAMS_comp_patchsum(nx,ny,1,npat,a)
      end select

      cdname='Ground longwave radiation'
      cdunits='W/m2'

   case ('snow_depth_p','snow_depth_ps')

      irecind = 1

      select case (trim(cvar))
      case ('snow_depth_ps')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('SNOW_DEPTH',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_sum_snowlayers(nx,ny,mynzs,npat,a(irecind))

      select case (trim(cvar))
      case ('snow_depth_p')
         ivar_type = 7
      case ('snow_depth_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nx,ny,1,npat,a)
      end select

      cdname='snow depth'
      cdunits='m'

   case ('snowcover_p','snowcover_ps')

      irecind = 1
      select case (trim(cvar))
      case ('snowcover_ps')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('SNOW_MOIST',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_sum_snowlayers(nx,ny,mynzs,npat,a(irecind))

      select case (trim(cvar))
      case ('snow_depth_p')
         ivar_type = 7
      case ('snowcover_ps')
         ivar_type = 2
         call RAMS_comp_patchsum_l(nx,ny,1,npat,a)
      end select

      cdname='snowcover'
      cdunits='kg/m2'

   case ('sltex','sltex_bp')

      irecind = 1
      select case (trim(cvar))
      case ('sltex_bp')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('SOIL_TEXT',idim_type,ngrd   &
           ,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      select case (trim(cvar))
      case ('sltex')
         ivar_type = 8
      case ('sltex_bp')
         ivar_type = 10
         call RAMS_comp_bigpatch(nx,ny,nsl,npat,a(irecind),a(1),b)
      end select

      cdname='soil textural class'
      cdunits='#'

   case ('smoist','smoist_ps')

      irecind = 1
      select case (trim(cvar))
      case ('smoist_ps')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select
      ierr = RAMS_getvar('SOIL_WATER',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr


      select case (trim(cvar))
      case ('smoist')
         ivar_type = 8
      case ('smoist_ps')
         ivar_type = 10
         call RAMS_comp_patchsum_l(nx,ny,nsl,npat,a)
      end select

      cdname='soil moisture'
      cdunits='m3/m3'


   case ('tsoil','tsoil_ps')

      irecind   = 1

      select case (trim(cvar))
      case ('tsoil_ps')
         irecsize  = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,a(irecind),b,flnm)
         irecind = irecind + irecsize
      end select

      ierr = RAMS_getvar('SOIL_ENERGY',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr

      ierr = RAMS_getvar('SOIL_WATER',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      ierr = RAMS_getvar('SOIL_TEXT',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr

      call RAMS_comp_copysst(nx,ny,nsl,a(irecind))

      call RAMS_comp_uextcm2tl(nx,ny,nsl,npat,a(irecind),scr%c,scr%d)


      select case (trim(cvar))
      case ('tsoil')
         call RAMS_comp_tempC(nx,ny,nsl,npat,a)
         ivar_type = 8
      case ('tsoil_ps')
         ivar_type = 10
         call RAMS_comp_patchsum(nx,ny,nsl,npat,a)
         call RAMS_comp_tempC(nx,ny,nsl,1,a)
      end select

      cdname='soil/sea temp'
      cdunits='C'

   case ('smfrac','smfrac_ps')

      irecind = 1
      select case (trim(cvar))
      case ('smfrac_ps')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
      end select
      irecind = irecind + irecsize
      ierr = RAMS_getvar('SOIL_WATER',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SOIL_TEXT',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      
      
      call rams_comp_slmstf(nx,ny,nsl,npat,a(irecind),scr%c)


      select case (trim(cvar))
      case ('smfrac')
         ivar_type = 8
      case ('smfrac_ps')
         ivar_type = 10
         call RAMS_comp_patchsum_l(nx,ny,nsl,npat,a)
      end select

      cdname='soil moisture frac'
      cdunits='m3/m3'

   case ('smpot','smpot_ps')

      irecind = 1
      select case (trim(cvar))
      case ('smpot_ps')
         irecsize = nx * ny * npat
         ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd   &
              ,a(irecind),b,flnm)
         ierr_getvar = ierr_getvar + ierr
      end select
      irecind = irecind + irecsize
      ierr = RAMS_getvar('SOIL_WATER',idim_type,ngrd,a(irecind),b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SOIL_TEXT',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      
      
      call rams_comp_smpot(nx,ny,nsl,npat,a(irecind),scr%c)


      select case (trim(cvar))
      case ('smpot')
         ivar_type = 8
      case ('smpot_ps')
         ivar_type = 10
         call RAMS_comp_patchsum_l(nx,ny,nsl,npat,a)
      end select

      cdname='soil matric potential'
      cdunits='MPa'

   case ('sfcw_temp')
      ivar_type=2
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SFCWATER_ENERGY',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SFCWATER_MASS',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SFCWATER_NLEV',idim_type,ngrd,scr%f,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_sfcwmeantemp(nx,ny,mynzs,npat,scr%c,scr%d,scr%e,scr%f,a)
      ierr_getvar = ierr_getvar + ierr
      cdname='Pond/snow mean temperature'
      cdunits='C'

   case ('sfcw_mass')
      ivar_type=2
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SFCWATER_MASS',idim_type,ngrd,scr%e,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SFCWATER_NLEV',idim_type,ngrd,scr%f,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_sfcwinteg(nx,ny,mynzs,npat,scr%c,scr%e,scr%f,a)
      call RAMS_comp_noneg(nx,ny,1,a)
      ierr_getvar = ierr_getvar + ierr
      cdname='Pond/snow mass'
      cdunits='kg/m2'

   case ('sfcw_depth')
      ivar_type=2
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,scr%c,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SFCWATER_DEPTH',idim_type,ngrd,scr%d,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('SFCWATER_NLEV',idim_type,ngrd,scr%f,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_sfcwinteg(nx,ny,mynzs,npat,scr%c,scr%d,scr%f,a)
      call RAMS_comp_noneg(nx,ny,1,a)
      ierr_getvar = ierr_getvar + ierr
      cdname='Pond/snow depth'
      cdunits='m'

   ! ------------------------ Stilt-RAMS coupling------------
   case ('afxu')
      ivar_type=3
      ierr= RAMS_getvar('AFXU',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='adv u flux'
      cdunits='kg/m^2s'

   case ('afxub')
      ivar_type=3
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('AFXUB',idim_type,ngrd,a,b,flnm)
      cdname='averaged adv u flux'
      cdunits='kg/m^2s'

   case ('afxv')
      ivar_type=3
      ierr= RAMS_getvar('AFXV',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='adv v flux'
      cdunits='kg/m^2s'

   case ('afxvb')
      ivar_type=3
      ierr= RAMS_getvar('AFXVB',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='averaged adv v flux'
      cdunits='kg/m^2s'

   case ('afxw')
      ivar_type=3
      ierr= RAMS_getvar('AFXW',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='adv w flux'
      cdunits='kg/m^2s'

   case ('afxwb')
      ivar_type=3
      ierr= RAMS_getvar('AFXWB',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='averaged adv W flux'
      cdunits='kg/m^2s'

   case ('sigwb')
      ivar_type=3
      ierr= RAMS_getvar('SIGWB',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='averaged sigma W'
      cdunits='m/s'



   !------------Grell cumulus scheme --------------------------


   case ('cuprliq')
      ivar_type=6
      ierr= RAMS_getvar('CUPRLIQ',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz*ncld,a)
      call RAMS_comp_mults(nx,ny,nz*ncld,a,1000.)
      cdname='Conv. water mixing ratio'
      cdunits='g/kg'

   case ('cuprice')
      ivar_type=6
      ierr= RAMS_getvar('CUPRICE',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz*ncld,a)
      call RAMS_comp_mults(nx,ny,nz*ncld,a,1000.)
      cdname='Conv. water mixing ratio'
      cdunits='g/kg'


   !----- Column integrated convective cloud. ---------------------------------------------!
   case ('conv_cloud')
      ivar_type=9


      ierr= RAMS_getvar('TOPT',idim_type,ngrd,scr%h,b,flnm)
      ierr_getvar = ierr_getvar + ierr


      call RAMS_comp_zero(nx,ny,nz*ncld,scr%e)                      ! e is total condensed
      ierr= RAMS_getvar('CUPRLIQ',idim_type,ngrd,scr%x,b,flnm)      ! x is liquid
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz*ncld,scr%e,scr%x) !
      ierr_getvar = ierr_getvar + ierr
      ierr= RAMS_getvar('CUPRICE',idim_type,ngrd,scr%x,b,flnm)      ! x is ice
      if(ierr == 0) call RAMS_comp_accum(nx,ny,nz*ncld,scr%e,scr%x) !
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_noneg(nx,ny,nz*ncld,scr%e)                     ! Ensure positive

      call RAMS_comp_dn0(nx,ny,nz,scr%x,scr%z,scr%d,scr%h,ngrd)     ! d is density

      !----- Find the bottom and top (entire column). -------------------------------------!
      call RAMS_comp_zero(nx,ny,1,scr%u) ! u is bottom
      call RAMS_comp_zero(nx,ny,1,scr%t) ! t is Top
      call RAMS_comp_adds(nx,ny,1,scr%t,100000.)
      !------------------------------------------------------------------------------------!


      call RAMS_comp_massint(nx,ny,nz,ncld,a,scr%d,scr%e,scr%u,scr%t,scr%h,ngrd)
      cdname  = 'Convective clouds'
      cdunits = 'kg/m2'
      !------------------------------------------------------------------------------------!

   case ('areadn')
      ivar_type=9
      ierr= RAMS_getvar('AREADN',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Downdraft relative area'
      cdunits=''

   case ('areaup')
      ivar_type=9
      ierr= RAMS_getvar('AREAUP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Updraft relative area'
      cdunits=''


   case ('ierr')
      ivar_type=9
      ierr= RAMS_getvar('XIERR',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Grell''s error flag:'
      cdunits=' '

   case ('upmf')
      ivar_type=9
      ierr= RAMS_getvar('UPMF',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='updraft mass flux'
      cdunits='kg/m2/s'

   case ('dnmf')
      ivar_type=9
      ierr= RAMS_getvar('DNMF',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='downdraft mass flux'
      cdunits='kg/m2/s'

   case ('upmx')
      ivar_type=9
      ierr= RAMS_getvar('UPMX',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='potential updraft mass flux'
      cdunits='kg/m2/s'

   case ('dnmx')
      ivar_type=9
      ierr= RAMS_getvar('DNMX',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='potential downdraft mass flux'
      cdunits='kg/m2/s'

   case ('wdndraft')
      ivar_type=9
      ierr= RAMS_getvar('WDNDRAFT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='downdraft velocity at origin of downdrafts'
      cdunits='m/s'

   case ('wupdraft')
      ivar_type=9
      ierr= RAMS_getvar('WUPDRAFT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='updraft velocity at origin of updrafts'
      cdunits='m/s'

   case ('wbuoymin')
      ivar_type=9
      ierr= RAMS_getvar('WBUOYMIN',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='minimum velocity to reach LCL'
      cdunits='m/s'

   case ('zklnb')
      ivar_type=9
      ierr= RAMS_getvar('ZKLNB',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Level of neutral buoyancy'
      cdunits='m'

   case ('zklfc')
      ivar_type=9
      ierr= RAMS_getvar('ZKLFC',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Level of free convection'
      cdunits='m'

   case ('zklcl')
      ivar_type=9
      ierr= RAMS_getvar('ZKLCL',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Lifting condensation level'
      cdunits='m'

   case ('zklod')
      ivar_type=9
      ierr= RAMS_getvar('ZKLOD',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Height of origin of downdraft'
      cdunits='m'

   case ('zklou')
      ivar_type=9
      ierr= RAMS_getvar('ZKLOU',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Height of origin of updraft'
      cdunits='m'

   case ('zkdet')
      ivar_type=9
      ierr= RAMS_getvar('ZKDET',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Top of downdraft detrainment'
      cdunits='m'

   case ('edt')
      ivar_type=9
      ierr= RAMS_getvar('EDT',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Downdraft/updraft ratio'
      cdunits=' '

   case ('aadn')
      ivar_type=9
      ierr= RAMS_getvar('AADN',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Dndraft cloud work function'
      cdunits='J/kg'

   case ('aaup')
      ivar_type=9
      ierr= RAMS_getvar('AAUP',idim_type,ngrd,a,b,flnm)
      ierr_getvar = ierr_getvar + ierr
      cdname='Updraft cloud work function'
      cdunits='J/kg'

   !---------------------------------------------------------------------------------------!
   !     Variables that are reduced to a given height, using the Monin-Obukhov similarity  !
   ! theory.  We also re-compute the "perceived" stars here in case we are running with    !
   ! ED-2.1.                                                                               !
   !---------------------------------------------------------------------------------------!
   case ('tempc2m','theta2m','rv2m','tdewc2m','rhum2m','co22m','u10m','rib_ps','zeta_ps'    &
        ,'ustar_ps','tstar_ps','rstar_ps','cstar_ps','usfc','vsfc','theiv2m')
      ivar_type = 2

      !----- Topography. ------------------------------------------------------------------!
      ierr = RAMS_getvar('TOPT',idim_type,ngrd,scr%c,b,flnm)  ! c = topography

      !----- Winds. -----------------------------------------------------------------------!
      ierr = RAMS_getvar('UP',idim_type,ngrd,scr%d,b,flnm)    ! d = zonal wind
      ierr_getvar = ierr_getvar + ierr
      ierr=RAMS_getvar('VP',idim_type,ngrd,scr%e,b,flnm)      ! e = meridional wind
      ierr_getvar = ierr_getvar + ierr
      call RAMS_comp_clone(nx,ny,nz,scr%y,scr%d)              ! y = zonal wind
      call RAMS_comp_clone(nx,ny,nz,scr%z,scr%e)              ! z = meridional wind
      call RAMS_comp_rotate(nx,ny,nz,scr%y,scr%z,ngrd)        ! rotate
      call RAMS_comp_avgu(nx,ny,nz,scr%y)                     ! y = true zonal wind
      call RAMS_comp_avgv(nx,ny,nz,scr%z)                     ! z = true meridional wind
      call RAMS_comp_speed(nx,ny,nz,scr%d,scr%e)              ! d = wind magnitude
      call RAMS_flush_to_zero(nx,ny,nz,npat,scr%d,ubmin)      ! d = wind magnitude

      !----- Atmospheric properties. ------------------------------------------------------!
      ierr = RAMS_getvar('THETA',idim_type,ngrd,scr%e,b,flnm) ! e = potential temperature.
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('RV',idim_type,ngrd,scr%f,b,flnm)    ! f = H2O mixing ratio
      ierr_getvar = ierr_getvar + ierr
      if (co2_on /= 0) then
         ierr= RAMS_getvar('CO2P',idim_type,ngrd,scr%g,b,flnm)
         ierr_getvar = ierr_getvar + ierr
      else
         write (unit=*,fmt='(a,1x,es12.5)') '       # Assigning constant CO2P =',myco2con(1)
         call ae0(nx*ny*nz,scr%g,myco2con(1))
      end if

      !----- Roughness. -------------------------------------------------------------------!
      ierr = RAMS_getvar('PATCH_ROUGH', idim_type,ngrd,scr%h,b,flnm) ! h = patch roughness.
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('RIBULK',idim_type,ngrd,scr%t,b,flnm)       ! t = bulk Ri
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('ZETA',idim_type,ngrd,scr%u,b,flnm)         ! u = z/L
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('PATCH_AREA',idim_type,ngrd,scr%m,b,flnm)   ! m = patch area.
      ierr_getvar = ierr_getvar + ierr

      !----- Canopy air space (CAS) properties. ----------------------------------------------!
      ierr = RAMS_getvar('CAN_THETA',idim_type,ngrd,scr%n,b,flnm) ! n = CAS potential temp.
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('CAN_RVAP',idim_type,ngrd,scr%o,b,flnm)  ! o = CAS mixing ratio.
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('CAN_PRSS',idim_type,ngrd,scr%s,b,flnm)  ! s = canopy pressure
      ierr_getvar = ierr_getvar + ierr
      ierr = RAMS_getvar('CAN_CO2',idim_type,ngrd,scr%w,b,flnm)   ! w = CAS CO2 mixing r.
      ierr_getvar = ierr_getvar + ierr

      !----- Characteristic scales (aka stars). ----------------------------------------------!
      ierr = RAMS_getvar('USTAR',idim_type,ngrd,scr%p,b,flnm)     ! p = ustar
      ierr_getvar = ierr_getvar + ierr
      call RAMS_flush_to_zero(nx,ny,1,npat,scr%p,ustmin)          ! p = ustar
      ierr = RAMS_getvar('TSTAR',idim_type,ngrd,scr%q,b,flnm)     ! q = tstar
      ierr_getvar = ierr_getvar + ierr
      call RAMS_flush_to_zero(nx,ny,1,npat,scr%q,1.e-6)           ! q = tstar
      ierr = RAMS_getvar('RSTAR',idim_type,ngrd,scr%r,b,flnm)     ! r = rstar
      ierr_getvar = ierr_getvar + ierr
      call RAMS_flush_to_zero(nx,ny,1,npat,scr%r,1.e-6)           ! r = rstar
      ierr = RAMS_getvar('CSTAR',idim_type,ngrd,scr%x,b,flnm)     ! x = cstar
      ierr_getvar = ierr_getvar + ierr
      call RAMS_flush_to_zero(nx,ny,1,npat,scr%x,1.e-6)           ! x = cstar
      select case (trim(cvar))
      case ('theta2m')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'THET',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'Potential temperature at 2m AGL'
         cdunits = 'K'

      case ('tempc2m')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'TEMP',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)
         call RAMS_comp_tempC(nx,ny,1,1,a)

         cdname  = 'Temperature at 2m AGL'
         cdunits = 'C'

      case ('theiv2m')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'THEE',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'Theta_Eiv at 2m AGL'
         cdunits = 'K'

      case ('tdewc2m')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'TDEW',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)
         call RAMS_comp_tempC(nx,ny,1,1,a)

         cdname  = 'Dew/frost point at 2m AGL'
         cdunits = 'C'

      case ('rv2m')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'RVAP',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)
         call RAMS_comp_mults(nx,ny,1,a,1000.)

         cdname  = 'Vapour mixing ratio at 2m AGL'
         cdunits = 'g/kg'

      case ('co22m')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'CO_2',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'CO2 mixing ratio at 2m AGL'
         cdunits = 'umol/mol'

      case ('rhum2m')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'TDEW',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'TEMP',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,scr%v)
         call RAMS_comp_relhum(nx,ny,1,a,scr%v)
         call RAMS_comp_mults(nx,ny,1,a,100.)

         cdname  = 'Relative humidity at 2m AGL'
         cdunits = '%'

      case ('zeta_ps')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'ZETA',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'Dimensionless height'
         cdunits = '---'

      case ('rib_ps')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'RICH',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'Bulk Richardson number'
         cdunits = '---'

      case ('u10m')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'WIND',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,10.,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'Wind speed at 10m AGL'
         cdunits = 'm/s'

      case ('usfc')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'WIND',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,10.,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)
         call RAMS_comp_uvsfc(nx,ny,'U',a,scr%y,scr%z)
         cdname  = 'Zonal wind at 10m AGL'
         cdunits = 'm/s'

      case ('vsfc')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'WIND',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,10.,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)
         call RAMS_comp_uvsfc(nx,ny,'V',a,scr%y,scr%z)
         cdname  = 'Zonal wind at 10m AGL'
         cdunits = 'm/s'

      case ('ustar_ps')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'USTR',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'Friction velocity'
         cdunits = 'm/s'

      case ('tstar_ps')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'TSTR',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'Potential temperature characteristic scale'
         cdunits = 'K'

      case ('rstar_ps')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'RSTR',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'Vapour mixing ratio characteristic scale'
         cdunits = 'g/kg'
         call RAMS_comp_mults(nx,ny,1,a,1.e3)

      case ('cstar_ps')
         call RAMS_reduced_prop(nx,ny,nz,npat,ngrd,'CSTR',scr%c,scr%e,scr%f,scr%g,scr%d    &
                               ,scr%n,scr%o,scr%w,scr%s,2.0,scr%h,scr%t,scr%u,scr%m,scr%p  &
                               ,scr%q,scr%r,scr%x,a)

         cdname  = 'CO2 mixing ratio characteristic scale'
         cdunits = 'umol/mol'

      end select


   case default

      write (unit=*,fmt='(2(a,1x))') '       # Variable name not found in RAMS_varlib -'   &
                                    ,trim(cvar)
      ivar_type = 0

   end select

   if (ierr_getvar > 0) then
     write (unit=*,fmt='(3(a,1x))') '       # WARNING! Not all the variables needed for'   &
                                   ,trim(cvar),' are available...' 
     ivar_type=0
   end if


   select case (ivar_type)
   case (2)
      call RAMS_show_range(nx,ny,1,1,a,cvar,'RAMS_varlib')
   case (3)
      call RAMS_show_range(nx,ny,nz,1,a,cvar,'RAMS_varlib')
   case (6)
      call RAMS_show_range(nx,ny,nz,ncld,a,cvar,'RAMS_varlib')
   case (7)
      call RAMS_show_range(nx,ny,npat,1,a,cvar,'RAMS_varlib')
   case (8)
      call RAMS_show_range(nx,ny,nsl,npat,a,cvar,'RAMS_varlib')
   case (9)
      call RAMS_show_range(nx,ny,ncld,1,a,cvar,'RAMS_varlib')
   case (10)
      call RAMS_show_range(nx,ny,nsl,1,a,cvar,'RAMS_varlib')
   end select

   return
end subroutine RAMS_varlib
!==========================================================================================!
!==========================================================================================!
