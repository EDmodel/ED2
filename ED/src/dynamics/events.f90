!!! Module for handling discrete events
!! evaluated daily

subroutine read_events_xml(filename,success)

  character*(*),intent(in) :: filename
  logical                  :: iamhere
  logical, intent(out)     :: success 

  success = .false.

  if (filename /= '') then
     inquire (file=filename,exist=iamhere)
     if (iamhere) then
        print*,"LOADING EVENT FILE ::",trim(filename)
        call libxml2f90__setformat(3) !set to pure XML but with blank removal
        call libxml2f90__readin_file(trim(filename),'EVENTS')
        success = .true.

        call libxml2f90__ll_selectlist('EVENTS')
        call libxml2f90__ll_selecttag('ACT','eventlist',1) !select upper level tag  
        call libxml2f90__ll_exist('DOWN','event',nevent)  
        print*,"NUMBER OF EVENTS FOUND", nevent

     end if
  endif

end subroutine read_events_xml


subroutine prescribed_event(year,doy)

  use misc_coms, only: event_file 

  integer, intent(in) :: year
  integer, intent(in) :: doy !! day of year
  logical,save :: first = .true.
  integer,save :: next_year = 100000
  integer,save :: next_doy = 367
  integer :: curr_year, curr_doy
  integer(4) :: i,j,nevent,nrep
  integer(4),dimension(20) ::ival,pft
  real(8),dimension(20):: rval
  logical(4) :: texist = .false.
  


  if (event_file /= '') then

     if(first) then  !! is first pass, init vars

        call read_events_xml(trim(event_file),first)

        if(first) then !! if read was successful

           !! access list and find next event
           call libxml2f90__ll_selectlist('EVENTS')  
           call libxml2f90__ll_selecttag('ACT','eventlist',1) !select upper level tag
           call libxml2f90__ll_exist('DOWN','event',nevent)           
           print*,"NUMBER OF EVENTS FOUND", nevent
           if(nevent .ge. 1) then
              do i=1,nevent
                 call libxml2f90__ll_selecttag('DOWN','event',i) !select event
                 call getConfigINT   ('year','event',i,curr_year,texist)
                 call getConfigINT   ('doy','event',i,curr_doy,texist)
                 print*,curr_year,curr_doy,texist,next_year,next_doy
                 !! see if the event bumps the closest event
                 if((curr_year == year .and. curr_doy .ge. doy) .or. (curr_year .gt. year)) then
                    if((curr_year .lt. next_year) .or. (curr_year == next_year .and. curr_doy .lt. next_doy)) then
                       next_year = curr_year
                       next_doy  = curr_doy
                    end if
                 end if
                 call libxml2f90__ll_selecttag('UP','eventlist',1) !move back up to top level
              enddo
              print*,"FIRST EVENT",next_doy,next_year
           end if
              
        endif

        first = .false.
     end if

     if(next_year == year .and. next_doy == doy) then

        !! access list and find todays events
        call libxml2f90__ll_selectlist('EVENTS')  
        call libxml2f90__ll_selecttag('ACT','eventlist',1) !select upper level tag
        call libxml2f90__ll_exist('DOWN','event',nevent)
 
!        print*,"EVENT",nevent
        findevent: do i=1,nevent !! Find Event in list
           call libxml2f90__ll_selecttag('DOWN','event',i) !select event
           call getConfigINT   ('year','event',i,curr_year,texist)
           call getConfigINT   ('doy','event',i,curr_doy,texist)
           !print*,curr_year,curr_doy,texist,next_year,next_doy
           !! see if the event bumps the closest event
           if((curr_year == year .and. curr_doy == doy)) then
              print*,"EVENT",i,"of",nevent
              
              !!! PROCESS AGRICULTURAL EVENTS

              !! PLANTING
              call libxml2f90__ll_exist('DOWN','plant',nrep)
              if(nrep .gt. 0) then
                 do j = 1,nrep
                    call libxml2f90__ll_selecttag('DOWN','plant',j)
                    call getConfigINT   ('pft','plant',j,pft(1),texist)
                    if(.not.texist) then
                       PRINT*,"Planting prescribed but PFT not specified"
                       stop
                    endif
                    call getConfigREAL   ('density','plant',j,rval(1),texist)
                    if(.not.texist) rval(1) = 1.0
                    call event_planting(pft(1),rval(1))
                    
                    call libxml2f90__ll_selecttag('UP','event',i)
                 end do
              end if

              !!! HARVESTING
              call libxml2f90__ll_exist('DOWN','harvest',nrep)
              if(nrep .gt. 0) then
                 print*,"HARVEST FESTIVAL"                 
                 ! agb_frac
                 ! bgb_frac 
                 do j = 1,nrep
                    call libxml2f90__ll_selecttag('DOWN','harvest',j)
                    !! Above-ground coarse biomass fraction removed
                    call getConfigREAL  ('agb_frac','harvest',j,rval(1),texist)
                    if(.not.texist) rval(1) = 1.0
                    !! Below-ground biomass fraction removed
                    call getConfigREAL   ('bgb_frac','harvest',j,rval(2),texist)
                    if(.not.texist) rval(2) = 0.0
                    !! Foliage fraction removed
                    call getConfigREAL   ('fol_frac','harvest',j,rval(3),texist)
                    if(.not.texist) rval(3) = 1.0
                    call event_harvest(rval(1),rval(2),rval(3))
                    
                    call libxml2f90__ll_selecttag('UP','event',i)
                 end do
              end if

              call libxml2f90__ll_exist('DOWN','fertilize',nrep)
              if(nrep .gt. 0) print*,"I say HABER, you say PROCESS"
              ! nitrogen

              call libxml2f90__ll_exist('DOWN','till',nrep)
              if(nrep .gt. 0) print*,"TILL"
              ! depth

              call libxml2f90__ll_exist('DOWN','pesticide',nrep)
              if(nrep .gt. 0) print*,"PESTICIDE"
              
              call libxml2f90__ll_exist('DOWN','herbicide',nrep)
              if(nrep .gt. 0) print*,"HERBICIDE"
              !pft !! spp removed

              call libxml2f90__ll_exist('DOWN','grazer',nrep)
              if(nrep .gt. 0) print*,"grazer"
              !pft !! spp removed
              !fol_frac



              !!! PROCESS FORESTRY MANAGEMENT
              call libxml2f90__ll_exist('DOWN','thin',nrep)
              if(nrep .gt. 0) print*,"There's too many trees around here"
              !pft, default = all
              
              !!! PROCESS NATURAL DISTURBANCES
              call libxml2f90__ll_exist('DOWN','wind',nrep)
              if(nrep .gt. 0) print*,"If the wind blows, your gone"
              !fol_frac
              !agb_frac

              call libxml2f90__ll_exist('DOWN','fire',nrep)
              if(nrep .gt. 0) print*,"FIRE, FIRE, FIRE"

              call libxml2f90__ll_exist('DOWN','pathogen',nrep)
              if(nrep .gt. 0) print*,"How can you not love the Beatles"

              call libxml2f90__ll_exist('DOWN','invasive',nrep)
              if(nrep .gt. 0) print*,"introduced species"

              
              call libxml2f90__ll_selecttag('UP','eventlist',1) !move back up to top level
              exit findevent
           endif
           call libxml2f90__ll_selecttag('UP','eventlist',1) !move back up to top level
        enddo findevent


        !! find next event
        next_year = 100000000
        next_doy  = 365
        do i=1,nevent
           call libxml2f90__ll_selecttag('DOWN','event',i) !select event
           call getConfigINT   ('year','event',i,curr_year,texist)
           call getConfigINT   ('doy','event',i,curr_doy,texist)
           
           !! see if the event bumps the closest event
           if((curr_year == year .and. curr_doy .gt. doy) .or. (curr_year .gt. year)) then
              if((curr_year .lt. next_year) .or. (curr_year == next_year .and. curr_doy .lt. next_doy)) then
                 next_year = curr_year
                 next_doy  = curr_doy
              end if
           end if
           call libxml2f90__ll_selecttag('UP','eventlist',1) !move back up to top level
        enddo
        print*,"NEXT EVENT",next_doy,next_year


     endif
  end if

end subroutine prescribed_event


subroutine event_harvest(agb_frac8,bgb_frac8,fol_frac8)
  use grid_coms, only : ngrids
  use ed_state_vars,only: edgrid_g, &
       edtype,polygontype,sitetype, &
       patchtype,allocate_patchtype,copy_patchtype,deallocate_patchtype 
  use pft_coms, only:sla,qsw,q,hgt_min
  use misc_coms, only: integration_scheme
  use disturbance_utils_ar,only: plant_patch_ar

  real(8),intent(in) :: agb_frac8
  real(8),intent(in) :: bgb_frac8
  real(8),intent(in) :: fol_frac8
  real :: agb,bgb
  integer :: ifm,ipy,isi,ipa
  type(edtype), pointer :: cgrid
  type(polygontype), pointer :: cpoly
  type(sitetype),pointer :: csite
  type(patchtype),pointer :: cpatch

  agb = real(agb_frac8)
  bgb = real(bgb_frac8)
  fol = real(fol_frac8)

  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"         HARVESTING "
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,""
  print*,agb,bgb

  do ifm = 1,ngrids
     cgrid => edgrid_g(ifm) 
     do ipy = 1,cgrid%npolygons
        
        cpoly => cgrid%polygon(ipy)
        
        do isi = 1,cpoly%nsites
           
           csite => cpoly%site(isi)

           do ipa=1,csite%npatches
              
              cpatch => csite%patch(ipa)

              !! remove agb + bgb + fol

              !! move residual frac to debris/litter pools

              !! reset patch and cohort vars

!              call plant_patch_ar(csite,ipa,pft,density,planting_ht,cpoly%lsl(isi))            
!              call update_patch_derived_props_ar(csite, cpoly%lsl(isi), cpoly%met(isi)%rhos,ipa)
!              call new_patch_sfc_props_ar(csite, ipa, cpoly%met(isi)%rhos)
           enddo

           ! Update site properties. ## THINK ABOUT WHAT TO SET FLAG##########
!           call update_site_derived_props_ar(cpoly,1,isi)           

        enddo

     end do
  end do

  ! Re-allocate integration buffer
  if(integration_scheme == 1) call initialize_rk4patches_ar(0)

end subroutine event_harvest


subroutine event_planting(pft,density8)
  use grid_coms, only : ngrids
  use ed_state_vars,only: edgrid_g, &
       edtype,polygontype,sitetype, &
       patchtype,allocate_patchtype,copy_patchtype,deallocate_patchtype 
  use pft_coms, only:sla,qsw,q,hgt_min
  use misc_coms, only: integration_scheme
  use disturbance_utils_ar,only: plant_patch_ar

  integer(4),intent(in) :: pft
  real(8),intent(in) :: density8
  real :: density
  integer :: ifm,ipy,isi,ipa
  type(edtype), pointer :: cgrid
  type(polygontype), pointer :: cpoly
  type(sitetype),pointer :: csite
  real, parameter :: planting_ht = 0.1

  density = real(density8)

  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"         PLANTING "
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,""
  print*,pft,density

  do ifm = 1,ngrids
     cgrid => edgrid_g(ifm) 
     do ipy = 1,cgrid%npolygons
        
        cpoly => cgrid%polygon(ipy)
        
        do isi = 1,cpoly%nsites
           
           csite => cpoly%site(isi)

           do ipa=1,csite%npatches
              
              call plant_patch_ar(csite,ipa,pft,density,planting_ht,cpoly%lsl(isi))            
              call update_patch_derived_props_ar(csite, cpoly%lsl(isi), cpoly%met(isi)%rhos,ipa)
              call new_patch_sfc_props_ar(csite, ipa, cpoly%met(isi)%rhos)
           enddo

           ! Update site properties. ## THINK ABOUT WHAT TO SET FLAG##########
           call update_site_derived_props_ar(cpoly,1,isi)           

        enddo

     end do
  end do

  ! Re-allocate integration buffer
  if(integration_scheme == 1) call initialize_rk4patches_ar(0)

end subroutine event_planting

!!! scrap code from writing event_planting
!!$  !!resize
!!$              ico = cpatch%ncohorts + 1
!!$print*,"Resize",ico
!!$              call allocate_patchtype(temppatch,cpatch%ncohorts)
!!$              call copy_patchtype(cpatch,temppatch,1,cpatch%ncohorts,1,cpatch%ncohorts)
!!$              call deallocate_patchtype(cpatch)
!!$              call allocate_patchtype(cpatch,ico)
!!$              call copy_patchtype(temppatch,cpatch,1,temppatch%ncohorts,1,temppatch%ncohorts)
!!$print*,"setvals"
!!$              !!enter values for new cohort
!!$              cpatch%pft(ico)     = pft
!!$              cpatch%nplant(ico)  = density
!!$              cpatch%hite(ico)    = hgt_min(pft)
!!$              cpatch%dbh(ico)     = h2dbh(hgt_min(pft),pft)
!!$              cpatch%bdead(ico)   = dbh2bd(cpatch%dbh(ico),&
!!$                   & hgt_min(pft),pft)
!!$              cpatch%bleaf(ico)   = dbh2bl(cpatch%dbh(ico),pft)
!!$print*,cpatch%hite(ico),cpatch%dbh(ico),cpatch%bdead(ico),cpatch%bleaf(ico)
!!$              cpatch%phenology_status(ico) = 0
!!$              cpatch%balive(ico)  = cpatch%bleaf(ico)* &
!!$                   & (1.0 + q(pft) + qsw(pft) * cpatch%hite(ico))
!!$              cpatch%lai(ico)     = cpatch%bleaf(ico) * density * sla(pft)
!!$              cpatch%bstorage(ico)  = 0.0
!!$              cpatch%veg_temp(ico)  = csite%can_temp(ipa)
!!$              cpatch%veg_water(ico) = 0.0
!!$print*,"vegtemp",cpatch%balive(ico),cpatch%veg_temp(ico)
!!$              !!update site lai
!!$              csite%lai(ipa) = sum(cpatch%lai(1:ico))
!!$print*,csite%lai(ipa),cpatch%lai
!!$              hcapveg = hcapveg_ref * max(cpatch%hite(1),cpatch%hite(ico),heathite_min) * cpatch%lai(ico)/csite%lai(ipa)
!!$
!!$print*,hcapveg
!!$              cpatch%veg_energy(ico) =  hcapveg * (cpatch%veg_temp(ico)-t3ple)
!!$print*,cpatch%veg_energy(ico)
!!$
!!$              !! count as new recruitment?
!!$              cpatch%new_recruit_flag(ico) = 1 
!!$print*,"init_cohorts"
!!$              call init_ed_cohort_vars_array(cpatch,ico,cpoly%lsl(isi))
!!$              
!!$              csite%cohort_count(ipa) = csite%cohort_count(ipa) + 1
!!$              cpatch%ncohorts = ico
!!$              csite%paco_n(ipa) = ico
!!$
!!$              ! Sort the cohorts so that the new cohort is at the correct height bin
!!$
!!$              call sort_cohorts_ar(cpatch)           
!!$   
!!$!              call update_patch_derived_props_ar(csite, cpoly%lsl(isi), cpoly%met(isi)%rhos, ipa)
!!$          
