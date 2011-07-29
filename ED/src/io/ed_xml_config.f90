!
!! Code for reading XML format config files for ED2
!!
!! Original code by Mike Dietze
!!
!! Main Functions:
!! 
!! read_ed_xml_config(filename,data)
!!    filename = character string
!!    data = ed_data object where parameters are set
!!   Intended to be called by ED_init.f90::EDinit(data) after hardwired defaults are read
!!   in order to set different/site-specific values and/or 
!!   when running batch/optimization/ensemble runs to reset parameters between runs
!!
!! TODO:
!! * write_ed_xml_config(filename,data) -> produce a record of parameters used in a run

subroutine count_pft_xml_config(filename,maxpft)
  integer(4) :: npft,i,myPFT
  integer :: maxpft
  logical(4) :: texist = .false.
  character*(*) :: filename
  character(1),allocatable    ::  readstring(:)
  allocate(readstring(10))
  deallocate(readstring)
print*,'count xml pft : ',trim(filename)
  !! Open File and Init
  call libxml2f90__setformat(1) !set to pure XML but with blank removal
  call libxml2f90__readin_file(trim(filename),'CONFIG')

  call libxml2f90__ll_selectlist('CONFIG')       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','pft',npft)    !get number of pft tags
  print*,"NPFT = ",npft
  maxpft = 0
  if(npft .ge. 1) then
     do i=1,npft
        call libxml2f90__ll_selecttag('DOWN','pft',i) !select pft
        call getConfigINT   ('num','pft',i,myPFT,texist)
        if(texist .and. myPFT > maxpft) maxpft = myPFT
        call libxml2f90__ll_selecttag('UP','config',1) !move back up to top level
     enddo
  endif

  !! next, check for a MAXPFT flag so don't have to read recursive files
  call libxml2f90__ll_exist('DOWN','maxpft',npft)    !get number of pft tags
  if(npft .ge. 1) then
     do i= 1,npft
        call getConfigINT   ('maxpft','config',i,myPFT,texist)  
        if(texist .and. myPFT > maxpft) maxpft = myPFT
     enddo
  endif

  print*,"MAXPFT = ",maxpft

  return

end subroutine count_pft_xml_config

recursive subroutine read_ed_xml_config(filename)

  use pft_coms
  use hydrology_coms
  use met_driver_coms
  use canopy_radiation_coms
  use disturb_coms
  use phenology_coms
  use physiology_coms
  use decomp_coms
  use fusion_fission_coms
  use grid_coms, only : ngrids
  use ed_max_dims, only : str_len

  use soil_coms  !, only: infiltration_method, dewmax, water_stab_thresh
!  use ed_data
  use ed_misc_coms!, only: ied_init_mode,ffilout,integration_scheme,sfilin,sfilout,thsums_database
  use rk4_coms, only : rk4min_veg_temp
  implicit none
  integer(4) :: i,npft,ntag,myPFT,len,ival = 0
  logical(4) :: texist = .false.
  real(8) :: rval
  character*(*) :: filename
  character(len=str_len)  :: cval
  integer             :: ng
!  type(eddata) :: data

  print*,"Opening ED2 XML Config file",trim(filename)

  !! Open File and Init
  call libxml2f90__setformat(1) !set to pure XML but with blank removal
  call libxml2f90__readin_file(trim(filename),trim(filename))!'CONFIG')

  !!********* EXTERN
  !! FIRST, check if any subfiles have to be loaded

  call libxml2f90__ll_selectlist(TRIM(FILENAME))       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','extern',ntag)    !get number of pft tags
  print*,"EXTERN READ FROM FILE ::",ntag

  if(ntag .ge. 1) then
     do i=1,ntag
!        call getConfigSTRING  ('extern','config',i,cval,texist)
!print*,cval,texist
        call libxml2f90__ll_selecttag('DOWN','extern',i)
        call libxml2f90__existid('extern',texist)
        if(texist) then 
           call libxml2f90__ll_getsize('extern',len)
           !----- MLO. Changed this to scalar so the interface check will work. -----------!
           call libxml2f90__ll_getch_scal('extern',len,cval)
           cval = cval(1:len)
           print*,"XML recursively loading ",trim(cval)
           call read_ed_xml_config(trim(cval))
        endif
        !! reset to original list
        call libxml2f90__ll_selectlist(TRIM(FILENAME))       
        call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
        call libxml2f90__ll_exist('DOWN','extern',ntag)    !get number of pft tags
     enddo
!!     stop
  endif
  


  !*******  MISC
  call libxml2f90__ll_selectlist(TRIM(FILENAME))       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','misc',ntag)    !get number of pft tags
  print*,"MISC READ FROM FILE ::",ntag
  if(ntag .ge. 1) then
     do i=1,ntag
        call libxml2f90__ll_selecttag('DOWN','misc',i)

         call getConfigINT  ('restart_mode','misc',i,ival,texist)
         if(texist) ied_init_mode = ival
         call getConfigSTRING  ('output_filepath','misc',i,cval,texist)
         if(texist) ffilout = trim(cval)
         call getConfigSTRING  ('input_filepath','misc',i,cval,texist)
         if (texist) thsums_database = trim(cval)

         call getConfigSTRING  ('history_out_filepath','misc',i,cval,texist)
         if(texist) sfilout = trim(cval)
         call getConfigINT  ('integration_scheme','misc',i,ival,texist)
         if(texist) integration_scheme = ival

         call libxml2f90__ll_selecttag('UP','config',1) !move back up to top level

      enddo
  endif


  !*******  ED_MISC
  call libxml2f90__ll_selectlist(TRIM(FILENAME))       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','ed_misc',ntag)    !get number of pft tags
  print*,"ED_MISC READ FROM FILE ::",ntag
  if(ntag .ge. 1) then
     do i=1,ntag
        call libxml2f90__ll_selecttag('DOWN','ed_misc',i)
        
        call getConfigINT  ('restart_target_year','ed_misc',i,ival,texist)
        if(texist) then
           restart_target_year = ival
           use_target_year = 1
        endif
        call getConfigINT  ('outputMonth','ed_misc',i,ival,texist)
        if(texist) outputMonth = ival
        call getConfigINT  ('burnin','ed_misc',i,ival,texist)
        if(texist) burnin = ival

        !! simple flags to turn on/off edaphic factors
        call getConfigINT  ('vary_elev','ed_misc',i,ival,texist)
        if(texist) vary_elev = ival
        call getConfigINT  ('vary_rad','ed_misc',i,ival,texist)
        if(texist) vary_rad = ival
        call getConfigINT  ('vary_hyd','ed_misc',i,ival,texist)
        if(texist) vary_hyd = ival
        
        
        
        call libxml2f90__ll_selecttag('UP','config',1) !move back up to top level
        
      enddo
  endif




  
  !! read PFT data
  call libxml2f90__ll_selectlist(TRIM(FILENAME))       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','pft',npft)    !get number of pft tags
  print*,"NPFT = ",npft
  print*,"PFT's READ FROM FILE ::",npft
  if(npft .ge. 1) then
     do i=1,npft
        call libxml2f90__ll_selecttag('DOWN','pft',i) !select pft

        call getConfigINT   ('num','pft',i,myPFT,texist)
        if(myPFT .le. n_pft .and. myPFT .ge. 1) then
           call getConfigINT('include_pft','pft',i,ival,texist)
           if(texist) then
              include_pft(myPFT) = ival == 1
!           else
!              include_pft(myPFT) = 1  !! if a PFT is defined, assume it's meant to be included
           endif
           call getConfigINT('include_pft_ag','pft',i,ival,texist)
           if(texist) include_pft_ag(myPFT) = ival == 1
           call getConfigSTRING('name','pft',i,cval,texist)
           call getConfigREAL  ('SLA','pft',i,rval,texist)
           if(texist) SLA(myPFT) = real(rval)
           call getConfigREAL  ('b1Bl','pft',i,rval,texist)
           if(texist) b1Bl(myPFT) = real(rval)
           call getConfigREAL  ('b2Bl','pft',i,rval,texist)
           if(texist) b2Bl(myPFT) = real(rval)
           call getConfigREAL  ('b1Bs','pft',i,rval,texist)
           if (texist) then
              b1Bs_small(myPFT) = real(rval)
              b1Bs_big(myPFT)   = real(rval)
           end if
           call getConfigREAL  ('b2Bs','pft',i,rval,texist)
           if (texist) then
              b2Bs_small(myPFT) = real(rval)
              b2Bs_big  (myPFT) = real(rval)
           end if
           call getConfigREAL  ('b1Ht','pft',i,rval,texist)
           if(texist) b1Ht(myPFT) = real(rval)
           call getConfigREAL  ('b2Ht','pft',i,rval,texist)
           if(texist) b2Ht(myPFT) = real(rval)
           call getConfigREAL  ('Vm0','pft',i,rval,texist)
           if(texist) Vm0(myPFT) = real(rval)
           call getConfigINT   ('phenology','pft',i,ival,texist)
           if(texist) phenology(myPFT) = ival
           call getConfigREAL  ('q','pft',i,rval,texist)
           if(texist) q(myPFT) = real(rval)
           call getConfigREAL  ('clumping','pft',i,rval,texist)
           if(texist) clumping_factor(myPFT) = rval
!!$           call getConfigREAL  ('beta0','pft',i,rval,texist)
!!$           if(texist) beta(myPFT,0) = real(rval)
!!$           call getConfigREAL  ('beta1','pft',i,rval,texist)
!!$           if(texist) beta(myPFT,1) = real(rval)
!!$           call getConfigREAL  ('beta2','pft',i,rval,texist)
!!$           if(texist) beta(myPFT,2) = real(rval)
!!$           call getConfigREAL  ('beta3','pft',i,rval,texist)
!!$           if(texist) beta(myPFT,3) = real(rval)
!!$           call getConfigREAL  ('beta4','pft',i,rval,texist)
!!$           if(texist) beta(myPFT,4) = real(rval)
           call getConfigREAL  ('leaf_width','pft',i,rval,texist)
           if(texist) leaf_width(myPFT) = real(rval)
           call getConfigREAL  ('hgt_min','pft',i,rval,texist)
           if(texist) hgt_min(myPFT) = real(rval)
           call getConfigREAL  ('plant_min_temp','pft',i,rval,texist)
           if(texist) plant_min_temp(myPFT) = real(rval)
           call getConfigREAL  ('mort3','pft',i,rval,texist)
           if(texist) mort3(myPFT) = real(rval)
           call getConfigREAL  ('nonlocal_dispersal','pft',i,rval,texist)
           if(texist) nonlocal_dispersal(myPFT) = real(rval)
           call getConfigREAL  ('seed_rain','pft',i,rval,texist)
           if(texist) seed_rain(myPFT) = real(rval)
!!$           call getConfigREAL  ('init_dens','pft',i,rval,texist)
!!$           if(texist) init_dens(myPFT) = real(rval)
           call getConfigREAL  ('stomatal_slope','pft',i,rval,texist)
           if(texist) stomatal_slope(myPFT) = real(rval)
           call getConfigREAL  ('growth_resp_factor','pft',i,rval,texist)
           if(texist) growth_resp_factor(myPFT) = real(rval)
           call getConfigREAL  ('r_fract','pft',i,rval,texist)
           if(texist) r_fract(myPFT) = real(rval)
!!$           call getConfigREAL  ('c_fract','pft',i,rval,texist)
!!$           if(texist) c_fract(myPFT) = real(rval)
           call getConfigREAL  ('repro_min_h','pft',i,rval,texist)
           if(texist) repro_min_h(myPFT) = real(rval)
           call getConfigREAL  ('treefall_gt','pft',i,rval,texist)
           if(texist) treefall_s_gtht(myPFT) = real(rval)
           call getConfigREAL  ('treefall_lt','pft',i,rval,texist)
           if(texist) treefall_s_ltht(myPFT) = real(rval)
           call getConfigREAL  ('dark_respiration_factor','pft',i,rval,texist)
           if(texist) dark_respiration_factor(myPFT) = real(rval)
           call getConfigREAL  ('qsw','pft',i,rval,texist)
           if(texist) qsw(myPFT) = real(rval)
           call getConfigREAL  ('sapwood_ratio','pft',i,rval,texist)
           if(texist) sapwood_ratio(myPFT) = real(rval)
           call getConfigREAL  ('c2n_leaf','pft',i,rval,texist)
           if(texist) c2n_leaf(myPFT) = real(rval)
           call getConfigREAL  ('c2n_stem','pft',i,rval,texist)
           if(texist) c2n_stem(myPFT) = real(rval)
           call getConfigREAL  ('c2n_recruit','pft',i,rval,texist)
           if(texist) c2n_recruit(myPFT) = real(rval)
!!$           call getConfigREAL  ('c2n_storage','pft',i,rval,texist)
!!$           if(texist) c2n_storage(myPFT) = real(rval)
           call getConfigREAL  ('max_dbh','pft',i,rval,texist)
           if(texist) max_dbh(myPFT) = real(rval)
           call getConfigREAL  ('rho','pft',i,rval,texist)
           if(texist) rho(myPFT) = real(rval)
           call getConfigREAL  ('D0','pft',i,rval,texist)
           if(texist) D0(myPFT) = real(rval)
           call getConfigREAL  ('mort1','pft',i,rval,texist)
           if(texist) mort1(myPFT) = real(rval)
           call getConfigREAL  ('mort2','pft',i,rval,texist)
           if(texist) mort2(myPFT) = real(rval)

           call getConfigREAL  ('Vm_low_temp','pft',i,rval,texist)
           if(texist) Vm_low_temp(myPFT) = real(rval)
           call getConfigREAL  ('cuticular_cond','pft',i,rval,texist)
           if(texist) cuticular_cond(myPFT) = real(rval)
           call getConfigREAL  ('quantum_efficiency','pft',i,rval,texist)
           if(texist) quantum_efficiency(myPFT) = real(rval)
           call getConfigINT  ('photosyn_pathway','pft',i,ival,texist)
           if(texist) photosyn_pathway(myPFT) = ival
           call getConfigREAL  ('leaf_turnover_rate','pft',i,rval,texist)
           if(texist) leaf_turnover_rate(myPFT) = real(rval)
           call getConfigREAL  ('root_turnover_rate','pft',i,rval,texist)
           if(texist) root_turnover_rate(myPFT) = real(rval)
           call getConfigREAL  ('storage_turnover_rate','pft',i,rval,texist)
           if(texist) storage_turnover_rate(myPFT) = real(rval)
           call getConfigREAL  ('root_respiration_factor','pft',i,rval,texist)
           if(texist) root_respiration_factor(myPFT) = real(rval)
           call getConfigREAL  ('seedling_mortality','pft',i,rval,texist)
           if(texist) seedling_mortality(myPFT) = real(rval)
           call getConfigREAL  ('water_conductance','pft',i,rval,texist)
           if(texist) water_conductance(myPFT) = real(rval)

!!! PFT VARIABLES THAT ARE ACTUALLY IN CANOPY_RADIATION
           call getConfigREAL  ('leaf_scatter_vis','pft',i,rval,texist)
           if(texist) leaf_scatter_vis(myPFT) = real(rval)
           call getConfigREAL  ('diffuse_backscatter_vis','pft',i,rval,texist)
           if(texist) leaf_backscatter_vis(myPFT) = real(rval)
           call getConfigREAL  ('emis_v','pft',i,rval,texist)
           if(texist) leaf_emis(myPFT) = rval
           

!!! PFT VARIABLES THAT ARE ACTUALLY IN DECOMP
           call getConfigREAL  ('f_labile','pft',i,rval,texist)
           if(texist) f_labile(myPFT) = real(rval)

           print*,i,myPFT,trim(cval),SLA(myPFT)
        else
           print*,"INVALID PFT READ FROM CONFIG FILE ::", myPFT
        endif
        call libxml2f90__ll_selecttag('UP','config',1) !move back up to top level
     enddo
  endif

  !********* READ PFT CONSTANTS (PFTCONST) PARMS
  call libxml2f90__ll_selectlist(TRIM(FILENAME))       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','pftconst',ntag)    !get number of pft tags
  print*,"PFTCONST READ FROM FILE ::",ntag
  if(ntag .ge. 1) then
     do i=1,ntag
        call libxml2f90__ll_selecttag('DOWN','pftconst',i)
        
        call getConfigREAL  ('c2n_slow','pftconst',i,rval,texist)
        if(texist) c2n_slow = real(rval)
        call getConfigREAL  ('c2n_structural','pftconst',i,rval,texist)
        if(texist) c2n_structural = real(rval)
!        call getConfigREAL  ('c2n_stem','pftconst',i,rval,texist)
!        if(texist) c2n_stem = real(rval)
        call getConfigREAL  ('l2n_stem','pftconst',i,rval,texist)
        if(texist) l2n_stem = real(rval)
        call getConfigREAL  ('C2B','pftconst',i,rval,texist)
        if(texist) C2B = real(rval)
        call getConfigREAL  ('agf_bs','pftconst',i,rval,texist)
        if(texist) agf_bs = real(rval)
        call getConfigREAL  ('frost_mort','pftconst',i,rval,texist)
        if(texist) frost_mort = real(rval)
        
        call libxml2f90__ll_selecttag('UP','config',1) !move back up to top level
     enddo
  endif
  

  !********* READ HYDROLOGY PARMS
  call libxml2f90__ll_selectlist(TRIM(FILENAME))       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','hydro',ntag)    !get number of pft tags
  print*,"HYDROLOGY READ FROM FILE ::",ntag
  if(ntag .ge. 1) then
     do i=1,ntag
        call libxml2f90__ll_selecttag('DOWN','hydro',i)

         call getConfigINT  ('useTOPMODEL','hydro',i,ival,texist)
         if(texist) useTOPMODEL = ival
         call getConfigINT  ('useRUNOFF','hydro',i,ival,texist)
         if(texist) useRUNOFF = ival
         call getConfigINT  ('outputPeriod','hydro',i,ival,texist)
         if(texist) HydroOutputPeriod = ival
         call getConfigREAL  ('MoistRateTuning','hydro',i,rval,texist)
         if(texist) MoistRateTuning = real(rval)
         call getConfigREAL  ('MoistSatThresh','hydro',i,rval,texist)
         if(texist) MoistSatThresh = real(rval)
         call getConfigREAL  ('MoistdWT','hydro',i,rval,texist)
         if(texist) Moist_dWT = real(rval)
         call getConfigREAL  ('FracLiqRunoff','hydro',i,rval,texist)
         if(texist) FracLiqRunoff = real(rval)
         call getConfigREAL  ('runoff_vmax','hydro',i,rval,texist)
         if(texist) runoff_vmax = real(rval)
         call getConfigREAL  ('GrassLAImax','hydro',i,rval,texist)
         if(texist) GrassLAImax = real(rval)
!! redundant with SOIL runoff_time 
         call getConfigREAL  ('inverse_runoff_time','hydro',i,rval,texist)
         if(texist) then
            print*,"hydro::inverse_runoff_time was found to be redundant with soil::runoff_time"
            print*,"Please update your xml config file"
            stop
         endif
 
         call libxml2f90__ll_selecttag('UP','config',1) !move back up to top level

      enddo
  endif

  !********* READ CLIMATE SCENARIO PARMS
  call libxml2f90__ll_selectlist(TRIM(FILENAME))       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','scenario',ntag)    !get number of pft tags
  print*,"SCENARIO READ FROM FILE ::",ntag
  if(ntag .ge. 1) then
     do i=1,ntag
        call libxml2f90__ll_selecttag('DOWN','scenario',i)
        
        call getConfigREAL  ('atm_tmp_intercept','scenario',i,rval,texist)
        if(texist) atm_tmp_intercept = real(rval)
        call getConfigREAL  ('atm_tmp_slope','scenario',i,rval,texist)
        if(texist) atm_tmp_slope = real(rval)
        call getConfigREAL  ('prec_intercept','scenario',i,rval,texist)
        if(texist) prec_intercept = real(rval)
        call getConfigREAL  ('prec_slope','scenario',i,rval,texist)
        if(texist) prec_slope = real(rval)
        call getConfigINT  ('humid_scenario','scenario',i,ival,texist)
        if(texist)  humid_scenario = real(ival)
        
        call libxml2f90__ll_selecttag('UP','config',1) !move back up to top level
     enddo
  endif
   

  !********* READ LAPSE PARMS
  call libxml2f90__ll_selectlist(TRIM(FILENAME))       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','lapse',ntag)    !get number of pft tags
  print*,"LAPSE READ FROM FILE ::",ntag
  if(ntag .ge. 1) then
     do i=1,ntag
        call libxml2f90__ll_selecttag('DOWN','lapse',i)

         call getConfigREAL  ('geoht','lapse',i,rval,texist)
         if(texist) lapse%geoht = real(rval)
         call getConfigREAL  ('vels','lapse',i,rval,texist)
         if(texist) lapse%vels = real(rval)
         call getConfigREAL  ('atm_tmp','lapse',i,rval,texist)
         if(texist) lapse%atm_tmp = real(rval)
         call getConfigREAL  ('rv','lapse',i,rval,texist)
         if(texist) lapse%atm_shv = real(rval)
         call getConfigREAL  ('prss','lapse',i,rval,texist)
         if(texist) lapse%prss = real(rval)
         call getConfigREAL  ('pcpg','lapse',i,rval,texist)
         if(texist) lapse%pcpg = real(rval)
         call getConfigREAL  ('atm_co2','lapse',i,rval,texist)
         if(texist) lapse%atm_co2 = real(rval)
         call getConfigREAL  ('rlong','lapse',i,rval,texist)
         if(texist) lapse%rlong = real(rval)
         call getConfigREAL  ('par_diffuse','lapse',i,rval,texist)
         if(texist) lapse%par_diffuse = real(rval)
         call getConfigREAL  ('par_beam','lapse',i,rval,texist)
         if(texist) lapse%par_beam = real(rval)
         call getConfigREAL  ('nir_diffuse','lapse',i,rval,texist)
         if(texist) lapse%nir_diffuse = real(rval)
         call getConfigREAL  ('nir_beam','lapse',i,rval,texist)
         if(texist) lapse%nir_beam = real(rval)
         call getConfigREAL  ('pptnorm','lapse',i,rval,texist)
         if(texist) lapse%pptnorm = real(rval)


         call libxml2f90__ll_selecttag('UP','config',1) !move back up to top level
      enddo
  endif

  !********  CANOPY_RADIATION
  call libxml2f90__ll_selectlist(TRIM(FILENAME))       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','radiation',ntag)    !get number of pft tags
  print*,"CANOPY RADIATION READ FROM FILE ::",ntag
  if(ntag .ge. 1) then
     do i=1,ntag

        call libxml2f90__ll_selecttag('DOWN','radiation',i)
        
        call getConfigREAL  ('rlong_min','radiation',i,rval,texist)
        if(texist) rlong_min = real(rval) 
        call getConfigREAL  ('veg_temp_min','radiation',i,rval,texist)
        if(texist) rk4min_veg_temp = rval ! This is double precision. 
        call getConfigREAL  ('mubar','radiation',i,rval,texist)
        if(texist) mubar = rval
        call getConfigREAL  ('visible_fraction','radiation',i,rval,texist)
        if(texist) visible_fraction = real(rval)
        call getConfigREAL  ('visible_fraction_dir','radiation',i,rval,texist)
        if(texist) visible_fraction_dir = real(rval)
        call getConfigREAL  ('visible_fraction_dif','radiation',i,rval,texist)
        if(texist) visible_fraction_dif = real(rval)
        call getConfigREAL  ('leaf_scatter_nir','radiation',i,rval,texist)
        if(texist) leaf_scatter_nir = real(rval)
        call getConfigREAL  ('leaf_reflect_nir','radiation',i,rval,texist)
        if(texist) leaf_reflect_nir = real(rval)
        call getConfigREAL  ('leaf_trans_nir','radiation',i,rval,texist)
        if(texist) leaf_trans_nir = real(rval)
        call getConfigREAL  ('diffuse_backscatter_vis','radiation',i,rval,texist)
        if(texist) leaf_backscatter_vis = real(rval)
        call getConfigREAL  ('diffuse_backscatter_nir','radiation',i,rval,texist)
        if(texist) leaf_backscatter_nir = real(rval)
        
        
        call libxml2f90__ll_selecttag('UP','config',1) !move back up to top level
     enddo
  endif


  !********* SOILS
  call libxml2f90__ll_selectlist(TRIM(FILENAME))       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','soil',ntag)    !get number of pft tags
  print*,"SOIL READ FROM FILE ::",ntag
  if(ntag .ge. 1) then
     do i=1,ntag

        call libxml2f90__ll_selecttag('DOWN','soil',i)
        
        call getConfigREAL  ('water_stab_thresh','soil',i,rval,texist)
        if(texist) water_stab_thresh = real(rval)
        call getConfigINT  ('infiltration_method','soil',i,ival,texist)
         if(texist) infiltration_method = ival
         call getConfigREAL  ('dewmax','soil',i,rval,texist)
         if(texist) dewmax = real(rval)

         call getConfigSTRING  ('vegetation_database','soil',i,cval,texist)
         if (texist) then
            do ng=1,ngrids
               veg_database(ng) = trim(cval)
            end do
         end if

         call getConfigSTRING  ('soil_database','soil',i,cval,texist)
         if (texist) then
            do ng=1,ngrids
               soil_database(ng) = trim(cval)
            end do
         end if

         call getConfigSTRING  ('soilstate_db','soil',i,cval,texist)
         if (texist) soilstate_db = trim(cval)

         call getConfigSTRING  ('soildepth_db','soil',i,cval,texist)
         if (texist) soildepth_db = trim(cval)

         call getConfigINT  ('isoilstateinit','soil',i,ival,texist)
         if(texist) isoilstateinit = ival
         call getConfigINT  ('isoildepthflg','soil',i,ival,texist)
         if(texist) isoildepthflg = ival

         !!CHECK FOR CONFLICT WITH INVERSE_RUNOFF_TIME
         call getConfigREAL  ('runoff_time','soil',i,rval,texist) 
         if(texist) runoff_time = real(rval)
      
        call libxml2f90__ll_selecttag('UP','config',1) !move back up to top level
     enddo
  endif
  


  !********* DECOMPOSITION
  call libxml2f90__ll_selectlist(TRIM(FILENAME))       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','decomposition',ntag)    !get number of pft tags
  print*,"DECOMPOSITION READ FROM FILE ::",ntag
  if(ntag .ge. 1) then
     do i=1,ntag
        
        call libxml2f90__ll_selecttag('DOWN','decomposition',i)
        
        call getConfigREAL  ('cwd_frac','decomposition',i,rval,texist)
        if(texist) cwd_frac = real(rval)

        call getConfigREAL  ('resp_opt_water','decomposition',i,rval,texist)
        if(texist)  resp_opt_water = real(rval)
        call getConfigREAL  ('resp_water_below_opt','decomposition',i,rval,texist)
        if(texist) resp_water_below_opt = real(rval)
        call getConfigREAL  ('resp_water_above_opt','decomposition',i,rval,texist)
        if(texist) resp_water_above_opt = real(rval)
        call getConfigREAL  ('resp_temperature_increase','decomposition',i,rval,texist)
        if(texist) resp_temperature_increase = real(rval)
        call getConfigREAL  ('N_immobil_supply_scale','decomposition',i,rval,texist)
        if(texist) N_immobil_supply_scale = real(rval)
        call getConfigREAL  ('r_fsc','decomposition',i,rval,texist)
        if(texist)  r_fsc = real(rval)
        call getConfigREAL  ('r_stsc','decomposition',i,rval,texist)
        if(texist)  r_stsc = real(rval)
        call getConfigREAL  ('r_ssc','decomposition',i,rval,texist)
        if(texist)  r_ssc = real(rval)
        call getConfigREAL  ('K1','decomposition',i,rval,texist)
        if(texist)  K1 = real(rval)
        call getConfigREAL  ('K2','decomposition',i,rval,texist)
        if(texist)  K2 = real(rval)
        call getConfigREAL  ('K3','decomposition',i,rval,texist)
        if(texist)  K3 = real(rval)
        call getConfigINT   ('N_decomp_lim','decomposition',i,ival,texist)
        if(texist)  N_decomp_lim= ival

        
        call libxml2f90__ll_selecttag('UP','config',1) !move back up to top level
     enddo
  endif


  !******** FUSION/FISSION
  call libxml2f90__ll_selectlist(TRIM(FILENAME))       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','fusefiss',ntag)    !get number of pft tags
  print*,"FUSEFISS READ FROM FILE ::",ntag
  if(ntag .ge. 1) then
     do i=1,ntag
        
        call libxml2f90__ll_selecttag('DOWN','fusefiss',i)
        
        call getConfigREAL  ('min_recruit_size','fusefiss',i,rval,texist)
        if(texist) min_recruit_size = real(rval)
        call getConfigREAL  ('fusetol_h','fusefiss',i,rval,texist)
        if(texist) fusetol_h = real(rval)
        call getConfigREAL  ('lai_fuse_tol','fusefiss',i,rval,texist)
        if(texist) lai_fuse_tol = real(rval)
        call getConfigREAL  ('lai_tol','fusefiss',i,rval,texist)
        if(texist) lai_tol = real(rval)
        call getConfigREAL  ('ntol','fusefiss',i,rval,texist)
        
        call libxml2f90__ll_selecttag('UP','config',1) !move back up to top level
     enddo
  endif
    

  !*********  DISTURBANCE
  call libxml2f90__ll_selectlist(TRIM(FILENAME))       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','disturbance',ntag)    !get number of pft tags
  print*,"DISTURBANCE READ FROM FILE ::",ntag
  if(ntag .ge. 1) then
     do i=1,ntag

        call libxml2f90__ll_selecttag('DOWN','disturbance',i)
        
        !! GENERAL
        call getConfigREAL  ('min_new_patch_area','disturbance',i,rval,texist)
        if(texist) min_new_patch_area = real(rval)
        call getConfigINT  ('include_fire','disturbance',i,ival,texist)
        if(texist) include_fire = ival
        call getConfigINT  ('ianth_disturb','disturbance',i,ival,texist)
        if(texist) ianth_disturb = ival
 
        !! TREEFALL
        call getConfigREAL  ('treefall_disturbance_rate','disturbance',i,rval,texist)
        if(texist) treefall_disturbance_rate = real(rval)
        
        call getConfigREAL  ('Time2Canopy','disturbance',i,rval,texist)
        if(texist) Time2Canopy = real(rval)
        
        call getConfigREAL  ('treefall_hite_threshold','disturbance',i,rval,texist)
        if(texist) treefall_hite_threshold = real(rval)

        !! FORESTRY
        call getConfigINT  ('plantation_year','disturbance',i,ival,texist)
        if(texist) plantation_year = ival
        call getConfigREAL  ('plantation_rotation','disturbance',i,rval,texist)
        if(texist) plantation_rotation = real(rval)
        call getConfigREAL  ('mature_harvest_age','disturbance',i,rval,texist)
        if(texist) mature_harvest_age = real(rval)
        !! Possibly DEPRECATED 
        call getConfigINT  ('forestry_on','disturbance',i,ival,texist)
        if(texist) forestry_on = ival
        call getConfigINT  ('agriculture_on','disturbance',i,ival,texist)
        if(texist) agriculture_on = ival
        
        !! FIRE
        call getConfigREAL  ('fire_dryness_threshold','disturbance',i,rval,texist)
        if(texist) fire_dryness_threshold = real(rval)
        call getConfigREAL  ('fire_parameter','disturbance',i,rval,texist)
        if(texist) fire_parameter = real(rval)

        call libxml2f90__ll_selecttag('UP','config',1) !move back up to top level
     enddo
  endif
  

  !********** PHENOLOGY
  call libxml2f90__ll_selectlist(TRIM(FILENAME))       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','phenology',ntag)    !get number of pft tags
  print*,"PHENOLOGY READ FROM FILE ::",ntag
  if(ntag .ge. 1) then
     do i=1,ntag
        
        call libxml2f90__ll_selecttag('DOWN','phenology',i)
        
        call getConfigREAL  ('retained_carbon_fraction','phenology',i,rval,texist)
        if(texist)  retained_carbon_fraction = real(rval)
        call getConfigREAL  ('theta_crit','phenology',i,rval,texist)
        if(texist)  theta_crit= real(rval)
        call getConfigREAL  ('dl_tr','phenology',i,rval,texist)
        if(texist)  dl_tr = real(rval)
        call getConfigREAL  ('st_tr1','phenology',i,rval,texist)
        if(texist)  st_tr1 = real(rval)
        call getConfigREAL  ('st_tr2','phenology',i,rval,texist)
        if(texist)  st_tr2 = real(rval)
        call getConfigREAL  ('phen_a','phenology',i,rval,texist)
        if(texist)  phen_a = real(rval)
        call getConfigREAL  ('phen_b','phenology',i,rval,texist)
        if(texist)  phen_b = real(rval)
        call getConfigREAL  ('phen_c','phenology',i,rval,texist)
        if(texist)  phen_c = real(rval)
        call getConfigINT  ('iphen_scheme','phenology',i,ival,texist)
        if(texist) iphen_scheme = ival

        call getConfigINT  ('iphenys1','phenology',i,ival,texist)
        if(texist) iphenys1 = ival
        call getConfigINT  ('iphenysf','phenology',i,ival,texist)
        if(texist) iphenysf = ival
        call getConfigINT  ('iphenyf1','phenology',i,ival,texist)
        if(texist) iphenyf1 = ival
        call getConfigINT  ('iphenyff','phenology',i,ival,texist)
        if(texist) iphenyff = ival

        
        call libxml2f90__ll_selecttag('UP','config',1) !move back up to top level
     enddo
  endif


  !********** PHYSIOLOGY
  call libxml2f90__ll_selectlist(TRIM(FILENAME))       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','physiology',ntag)    !get number of pft tags
  print*,"PHYSIOLOGY READ FROM FILE ::",ntag
  if(ntag .ge. 1) then
     do i=1,ntag
        
        call libxml2f90__ll_selecttag('DOWN','physiology',i)

        call getConfigINT  ('istoma_scheme','physiology',i,ival,texist)
        if(texist) istoma_scheme = ival
        call getConfigINT  ('n_plant_lim','physiology',i,ival,texist)
        if(texist) n_plant_lim = ival
        
        call libxml2f90__ll_selecttag('UP','config',1) !move back up to top level
     enddo
  endif


  !********** INITIAL CONDITIONS
  call libxml2f90__ll_selectlist(TRIM(FILENAME))       
  call libxml2f90__ll_selecttag('ACT','config',1) !select upper level tag
  call libxml2f90__ll_exist('DOWN','initcond',ntag)    !get number of pft tags
  init_fsc = -1.0
  init_stsc = -1.0
  init_ssc = -1.0
  init_stsl = -1.0
  init_fsn = -1.0
  init_msn = -1.0
  
  print*,"INITCOND READ FROM FILE ::",ntag
  if(ntag .ge. 1) then
     do i=1,ntag
        
        call libxml2f90__ll_selecttag('DOWN','initcond',i)

        call getConfigREAL  ('fsc','initcond',i,rval,texist)
        if(texist) init_fsc = rval
        call getConfigREAL  ('stsc','initcond',i,rval,texist)
        if(texist) init_stsc = rval
        call getConfigREAL  ('ssc','initcond',i,rval,texist)
        if(texist) init_ssc = rval
        call getConfigREAL  ('stsl','initcond',i,rval,texist)
        if(texist) init_stsl = rval
        call getConfigREAL  ('fsn','initcond',i,rval,texist)
        if(texist) init_fsn = rval
        call getConfigREAL  ('msn','initcond',i,rval,texist)
        if(texist) init_msn = rval
        

        
        call libxml2f90__ll_selecttag('UP','config',1) !move back up to top level
     enddo
  endif
  


  !********* LOST PARAMETERS
!!$         call getConfigINT  ('estimator','numerics',i,ival,texist)
!!$         if(texist) data%estimator = ival
!!$         call getConfigREAL  ('Vm_amp','numerics',i,rval,texist)
!!$         if(texist) data%Vm_amp = real(rval)
!!$         call getConfigINT  ('n_poi','numerics',i,ival,texist)
!!$         if(texist) data%n_poi = ival
!!$         call getConfigINT  ('max_cohorts','numerics',i,ival,texist)
!!$         if(texist) data%max_cohorts = ival
!!$         call getConfigREAL  ('dbhmax','numerics',i,rval,texist)
!!$         if(texist) data%dbhmax = real(rval)
!!$         call getConfigREAL  ('f_area','community',i,rval,texist)
!!$         if(texist) data%f_area = real(rval)
!!$         call getConfigINT  ('forecast_beg_year','community',i,ival,texist)
!!$         if(texist) data%forecast_beg_year = ival
!!$         call getConfigINT  ('project_landuse_rates','community',i,ival,texist)
!!$         if(texist) data%project_landuse_rates = ival
!!$         call getConfigREAL  ('k_fw','community',i,rval,texist)
!!$         if(texist) data%k_fw = real(rval)
!!$         call getConfigREAL  ('fp1','community',i,rval,texist)
!!$         if(texist) data%fp1 = real(rval)
!!$         call getConfigREAL  ('btol','community',i,rval,texist)
!!$         if(texist) data%btol = real(rval)
!!$         call getConfigREAL  (' theta_crit_lo','community',i,rval,texist)
!!$         if(texist) data%theta_crit_lo = real(rval)
!!$         call getConfigREAL  (' theta_crit_hi','community',i,rval,texist)
!!$         if(texist) data%theta_crit_hi= real(rval)
!!$         call getConfigREAL  ('gee_delay_time','community',i,rval,texist)
!!$         if(texist) data%gee_delay_time= real(rval)
!!$         call getConfigREAL  ('k_nlr','ecosystem',i,rval,texist)
!!$         if(texist) data%k_nlr = real(rval)
!!$         call getConfigREAL  ('nlr_low_temp','ecosystem',i,rval,texist)
!!$         if(texist) data%nlr_low_temp = real(rval)
!!$         call getConfigREAL  ('nlr_slope','ecosystem',i,rval,texist)
!!$         if(texist) data%nlr_slope = real(rval)
!!$         call getConfigINT  ('winter_decay','ecosystem',i,ival,texist)
!!$         if(texist) data%winter_decay = ival
!!$         call getConfigREAL  ('resp_opt_temp','ecosystem',i,rval,texist)
!!$         if(texist) data%resp_opt_temp = real(rval)
!!$         call getConfigREAL  ('resp_tshr','ecosystem',i,rval,texist)
!!$         if(texist) data%resp_tshr = real(rval)
!!$         call getConfigREAL  ('resp_tshl','ecosystem',i,rval,texist)
!!$         if(texist) data%resp_tshl = real(rval)
!!$         call getConfigREAL  ('resp_max_temp','ecosystem',i,rval,texist)
!!$         if(texist) data%resp_max_temp = real(rval)
!!$         call getConfigREAL  ('resp_opt_water','ecosystem',i,rval,texist)
!!$         if(texist) data%resp_opt_water = real(rval)
!!$         call getConfigREAL  ('resp_water_decay','ecosystem',i,rval,texist)
!!$         if(texist) data%resp_water_decay = real(rval)
!!$         call getConfigREAL  ('resp_water_growth','ecosystem',i,rval,texist)
!!$         if(texist) data%resp_water_growth = real(rval)
!!$         call getConfigREAL  ('nitrogen2','ecosystem',i,rval,texist)
!!$         if(texist) data%nitrogen2 = real(rval)
!!$         call getConfigREAL  ('nitrogen1','ecosystem',i,rval,texist)
!!$         if(texist) data%nitrogen1 = real(rval)
!!$
!!$         call getConfigINT  ('phenol_spring_firstyear','settings',i,ival,texist)
!!$         if(texist) data%phenol_spring_firstyear = ival
!!$         call getConfigINT  ('phenol_spring_lastyear','settings',i,ival,texist)
!!$         if(texist) data%phenol_spring_lastyear = ival
!!$         call getConfigINT  ('phenol_fall_firstyear','settings',i,ival,texist)
!!$         if(texist) data%phenol_fall_firstyear = ival
!!$         call getConfigINT  ('phenol_fall_lastyear','settings',i,ival,texist)
!!$         if(texist) data%phenol_fall_lastyear = ival
!!$
!!$         call getConfigINT  ('landuse_flag','settings',i,ival,texist)
!!$         if(texist) data%landuse_flag = ival
!!$         call getConfigINT  ('physiology_flag','settings',i,ival,texist)
!!$         if(texist) data%physiology_flag = ival
!!$         call getConfigINT  ('phenology_flag','settings',i,ival,texist)
!!$         if(texist) data%phenology_flag = ival
!!$         call getConfigINT  ('mortality_flag','settings',i,ival,texist)
!!$         if(texist) data%mortality_flag = ival
!!$         call getConfigINT  ('plant_respiration_flag','settings',i,ival,texist)
!!$         if(texist) data%plresp = ival
!!$         call getConfigSTRING  ('meteorology_filepath','settings',i,cval,texist)
!!$         if(texist) data%meteorology_filepath = trim(cval)
!!$         call getConfigSTRING  ('optimizer_filepath','settings',i,cval,texist)
!!$         if(texist) data%opt_inputs_fn = trim(cval)
!!$         call getConfigREAL  ('grid_resolution','settings',i,rval,texist)
!!$         if(texist) data%grid_resolution = real(rval)


       

end subroutine read_ed_xml_config

subroutine getConfigINT(ctag,ptag,pwhich,val,texist)
  character(*) :: ctag,ptag
  integer(4) :: val,j,pwhich
  logical(4) :: texist
!! Arguments:  same for REAL and STRING versions
!!  ctag = child tag (string)
!!  ptag = parent tag (string)
!!  pwhich = which node with matching parent tag (integer)
!!  val = value extracted from file
!!  texist = did the value exist? (logical)

  j=0
  texist = .false.
  
  call libxml2f90__ll_exist('DOWN',ctag,j)
  if(j .gt. 0) then
     call libxml2f90__ll_selecttag('DOWN',ctag,1)
     call libxml2f90__existid(ctag,texist)
     if(texist) call libxml2f90__ll_geti4(ctag,1,val)
     call libxml2f90__ll_selecttag('UP',ptag,pwhich)
  endif

end subroutine getConfigINT

subroutine getConfigSTRING(ctag,ptag,pwhich,val,texist)
  character(*) :: ctag,ptag,val
  integer(4) :: j,len,pwhich
  logical(4) :: texist
  j=0
  texist = .false.
  
  call libxml2f90__ll_exist('DOWN',ctag,j)
  if(j .gt. 0) then
     call libxml2f90__ll_selecttag('DOWN',ctag,1)
     call libxml2f90__existid(ctag,texist)
     if(texist) then 
        call libxml2f90__ll_getsize(ctag,len)
        call libxml2f90__ll_getch(ctag,len,val)
        val = val(1:len)
     endif
     call libxml2f90__ll_selecttag('UP',ptag,pwhich)
  endif

end subroutine getConfigSTRING

subroutine getConfigREAL(ctag,ptag,pwhich,val,texist)
  character(*) :: ctag,ptag
  integer(4) :: j,pwhich
  real(8)    :: val
  logical(4) :: texist
  j=0
  texist = .false.
  
  call libxml2f90__ll_exist('DOWN',ctag,j)
  if(j .gt. 0) then
     call libxml2f90__ll_selecttag('DOWN',ctag,1)
     call libxml2f90__existid(ctag,texist)
     if(texist) call libxml2f90__ll_getr8(ctag,1,val)
     call libxml2f90__ll_selecttag('UP',ptag,pwhich)
  endif

end subroutine getConfigREAL

subroutine write_ed_xml_config
!!produce a record of parameters used in a run
  use ed_max_dims, only: n_pft
  use pft_coms
  use canopy_radiation_coms
  use decomp_coms
  use ed_misc_coms, only: sfilout
  integer :: ival
  character(512) :: xfilout 

  write(xfilout,"(a)") trim(sfilout)//".xml"

  !! construct list
  call libxml2f90_ll_add_list("OUTCONFIG")
  call libxml2f90_ll_opentag("config")

  !************   PFT  *****************
  do i=1,n_pft
     if(include_pft(i)) then !! loop over included PFT's

        call libxml2f90_ll_opentag("pft")
     
        call putConfigINT("num",i)
        if (include_pft_ag(i)) then
           ival = 1
        else
           ival = 0
        end if
        call putConfigINT("include_pft_ag",ival)
        if (include_pft(i)) then
           ival = 1
        else
           ival = 0
        end if
        call putConfigINT("include_pft",ival)
        
        call putConfigREAL("SLA",SLA(i))
        call putConfigREAL("b1Bl",b1Bl(i))
        call putConfigREAL("b2Bl",b2Bl(i))
        call putConfigREAL("b1Bs",b1Bs_small(i))
        b1Bs_big(i) = b1Bs_big(i)
        call putConfigREAL("b2Bs",b2Bs_small(i))
        b2Bs_big(i) = b2Bs_big(i)
        call putConfigREAL("b1Ht",b1Ht(i))
        call putConfigREAL("b2Ht",b2Ht(i))
        call putConfigREAL("Vm0",Vm0(i))
        call putConfigINT("phenology",phenology(i))
        call putConfigREAL("q",q(i))
        call putConfigREAL8("clumping",clumping_factor(i))
        call putConfigREAL("leaf_width",leaf_width(i))
        call putConfigREAL("hgt_min",hgt_min(i))
        call putConfigREAL("plant_min_temp",plant_min_temp(i))
        call putConfigREAL("mort3",mort3(i))
        call putConfigREAL("nonlocal_dispersal",nonlocal_dispersal(i))
        call putConfigREAL("seed_rain",seed_rain(i))
        call putConfigREAL("stomatal_slope",stomatal_slope(i))
        call putConfigREAL("growth_resp_factor",growth_resp_factor(i))
        call putConfigREAL("r_fract",r_fract(i))
        call putConfigREAL("repro_min_h",repro_min_h(i))
        call putConfigREAL("treefall_gt",treefall_s_gtht(i))
        call putConfigREAL("treefall_lt",treefall_s_ltht(i))
        call putConfigREAL("dark_respiration_factor",dark_respiration_factor(i))
        call putConfigREAL("qsw",qsw(i))
        call putConfigREAL("c2n_leaf",c2n_leaf(i))
        call putConfigREAL("c2n_recruit",c2n_recruit(i))
        call putConfigREAL("max_dbh",max_dbh(i))
        call putConfigREAL("rho",rho(i))
        call putConfigREAL("D0",D0(i))
        call putConfigREAL("mort1",mort1(i))
        call putConfigREAL("mort2",mort2(i))
        call putConfigREAL("Vm_low_temp",Vm_low_temp(i))
        call putConfigREAL("cuticular_cond",cuticular_cond(i))
        call putConfigREAL("quantum_efficiency",quantum_efficiency(i))
        call putConfigINT("photosyn_pathway",photosyn_pathway(i))
        call putConfigREAL("leaf_turnover_rate",leaf_turnover_rate(i))
        call putConfigREAL("root_turnover_rate",root_turnover_rate(i))
        call putConfigREAL("storage_turnover_rate",storage_turnover_rate(i))
        call putConfigREAL("root_respiration_factor",root_respiration_factor(i))
        call putConfigREAL("seedling_mortality",seedling_mortality(i))
        call putConfigREAL("water_conductance",water_conductance(i))
        call putConfigREAL8("leaf_scatter_vis",leaf_scatter_vis(i))
        call putConfigREAL8("diffuse_backscatter_vis",leaf_backscatter_vis(i))
        call putConfigREAL8("emis_v",leaf_emis(i))
        call putConfigREAL("f_labile",f_labile(i))
        call libxml2f90_ll_closetag("pft")
     endif
  end do
 
  call libxml2f90_ll_closetag("config")
  !! write list
  print*,"Wrote Config Record:", trim(xfilout)
  open(12,file=trim(xfilout),form='formatted',status='replace')
  call libxml2f90__ll_report("OUTCONFIG",12,.false.)
  close(12)

end subroutine write_ed_xml_config


subroutine putConfigSTRING(tag,value)
  character(*),intent(in) :: tag 
  character(*),intent(in) :: value
  integer :: lenval 
  lenval = len(value)
  call libxml2f90_ll_opentag(tag)
  call libxml2f90_ll_addid(trim(tag),lenval,trim(value))
  call libxml2f90_ll_closetag(tag)
end subroutine putConfigSTRING

subroutine putConfigINT(tag,ivalue)
  use ed_max_dims, only : str_len
  character(*),intent(in) :: tag 
  integer,intent(in) :: ivalue
  character(str_len) :: value
  integer :: lenval 
  write(value,"(i11.1)") ivalue
  lenval = len(trim(value))
  call libxml2f90_ll_opentag(tag)
  call libxml2f90_ll_addid(trim(tag),lenval,trim(value))
  call libxml2f90_ll_closetag(tag)
end subroutine putConfigINT

subroutine putConfigREAL(tag,rvalue)
  use ed_max_dims, only : str_len
  character(*),intent(in) :: tag 
  real,intent(in) :: rvalue
  character(str_len) :: value
  integer :: lenval 
  write(value,"(f20.10)") rvalue
  lenval = len(trim(value))
  call libxml2f90_ll_opentag(tag)
  call libxml2f90_ll_addid(trim(tag),lenval,trim(value))
  call libxml2f90_ll_closetag(tag)
end subroutine putConfigREAL

subroutine putConfigREAL8(tag,rvalue)
  use ed_max_dims, only : str_len
  character(*),intent(in) :: tag 
  real(kind=8),intent(in) :: rvalue
  character(str_len) :: value
  integer :: lenval 
  write(value,"(f40.20)") rvalue
  lenval = len(trim(value))
  call libxml2f90_ll_opentag(tag)
  call libxml2f90_ll_addid(trim(tag),lenval,trim(value))
  call libxml2f90_ll_closetag(tag)
end subroutine putConfigREAL8

