!!! Module for handling discrete events
!! evaluated daily

subroutine read_events_xml(filename,success)
  implicit none
  character*(*),intent(in) :: filename
  logical                  :: iamhere
  logical, intent(out)     :: success
  integer                  :: nevent
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
  use ed_misc_coms, only: event_file
  implicit none

  integer, intent(in) :: year
  integer, intent(in) :: doy !! day of year
  logical,save :: first = .true.
  integer,save :: next_year = 100000
  integer,save :: next_doy = 367
  integer :: curr_year, curr_doy
  integer(kind=4) :: i,j,nevent,nrep
  integer(kind=4),dimension(20) ::pft
  real(kind=8),dimension(20):: rval
  logical(kind=4) :: texist = .false.



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
                    if(.not.texist) rval(1) = 1.0d+0
                    !! Below-ground biomass fraction removed
                    call getConfigREAL   ('bgb_frac','harvest',j,rval(2),texist)
                    if(.not.texist) rval(2) = 1.0d+0
                    !! Foliage fraction removed
                    call getConfigREAL   ('fol_frac','harvest',j,rval(3),texist)
                    if(.not.texist) rval(3) = 1.0d+0
                    !! storage fraction removed
                    call getConfigREAL   ('stor_frac','harvest',j,rval(4),texist)
                    if(.not.texist) rval(4) = 1.0d+0

                    call event_harvest(rval(1),rval(2),rval(3),rval(4))

                    call libxml2f90__ll_selecttag('UP','event',i)
                 end do
              end if

              call libxml2f90__ll_exist('DOWN','fertilize',nrep)
              if(nrep .gt. 0) then
                 print*,"I say HABER, you say PROCESS"
                 do j = 1,nrep

                    call libxml2f90__ll_selecttag('DOWN','fertilize',j)
                    rval = 0.0d+0
                    ! nitrogen
                    call getConfigREAL  ('NH4','fertilize',j,rval(1),texist)
                    if(.not.texist) rval(1) = 0.0d+0
                    call getConfigREAL  ('N03','fertilize',j,rval(2),texist)
                    if(.not.texist) rval(2) = 0.0d+0
                    ! Phosphorus
                    call getConfigREAL  ('P','fertilize',j,rval(3),texist)
                    if(.not.texist) rval(3) = 0.0d+0
                    ! Potassium
                    call getConfigREAL  ('K','fertilize',j,rval(4),texist)
                    if(.not.texist) rval(4) = 0.0d+0
                    ! Calcium
                    call getConfigREAL  ('Ca','fertilize',j,rval(5),texist)
                    if(.not.texist) rval(5) = 0.0d+0

                    call event_fertilize(rval(1:5))

                    call libxml2f90__ll_selecttag('UP','event',i)
                 end do
              endif

              call libxml2f90__ll_exist('DOWN','till',nrep)
              if(nrep .gt. 0) then
                 print*,"TILL"
                 do j = 1,nrep
                    ! depth (m)
                    call libxml2f90__ll_selecttag('DOWN','till',j)

                    call getConfigREAL  ('depth','till',j,rval(1),texist)
                    if(.not.texist) rval(1) = 0.15d+0 !assume a default of 15cm

                    call event_till(rval(1))

                    call libxml2f90__ll_selecttag('UP','event',i)

                 enddo
              end if

              call libxml2f90__ll_exist('DOWN','irrigate',nrep)
              if(nrep .gt. 0) then
                 print*,"IRRIGATE"
                 do j = 1,nrep

                    call getConfigREAL ('irrigate','event',j,rval(1),texist)
                    if(.not.texist) rval(1) = 0.0d+0

                    call event_irrigate(rval(1))

                 enddo
              endif

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
              if(nrep .gt. 0) then
                 print*,"FIRE, FIRE, FIRE"
                 do j = 1,nrep
                    call libxml2f90__ll_selecttag('DOWN','fire',j)

                    call event_fire()

                    call libxml2f90__ll_selecttag('UP','event',i)
                 end do
              end if

              call libxml2f90__ll_exist('DOWN','pathogen',nrep)
              if(nrep .gt. 0) print*,"How can you not love the Beatles"

              call libxml2f90__ll_exist('DOWN','invasive',nrep)
              if(nrep .gt. 0) print*,"introduced species"

              !! PLANTING (assumed to occur AFTER the other events)
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
                    if(.not.texist) rval(1) = 1.0d+0
                    call event_planting(pft(1),rval(1))

                    call libxml2f90__ll_selecttag('UP','event',i)
                 end do
              end if



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


subroutine event_harvest(agb_frac8,bgb_frac8,fol_frac8,stor_frac8)
  use stable_cohorts, only : is_resolvable
  use update_derived_utils, only : update_patch_derived_props, update_site_derived_props
  use grid_coms, only : ngrids
  use ed_state_vars,only: edgrid_g, &
       edtype,polygontype,sitetype, &
       patchtype,allocate_patchtype,copy_patchtype,deallocate_patchtype 
  use met_driver_coms, only : met_driv_state             ! ! structure
  use pft_coms, only:qsw,qbark,q,hgt_min, agf_bs, is_grass
  use ed_therm_lib, only: calc_veg_hcap,update_veg_energy_cweh
  use fuse_fiss_utils, only: terminate_cohorts
  use allometry, only : bd2dbh, dbh2h, bl2dbh, bl2h, h2dbh, area_indices &
                      , ed_biomass,size2bt,size2xb,ed_balive
  use consts_coms, only : pio4
  use ed_misc_coms     , only : igrass,frqsumi               ! ! intent(in)
  use plant_hydro,     only : rwc2tw,twi2twe                   ! ! sub-routine
  implicit none
  real(kind=8),intent(in) :: agb_frac8
  real(kind=8),intent(in) :: bgb_frac8
  real(kind=8),intent(in) :: fol_frac8
  real(kind=8),intent(in) :: stor_frac8
  real :: ialloc,bdeada_new,bdeadb_new,bsapa_new,bsapb_new,bbarka_new,bbarkb_new
  real :: bleaf_new,bfr_new,bstore_new
  real :: agb_frac,bgb_frac,fol_frac,stor_frac
  real :: old_leaf_hcap
  real :: old_wood_hcap
  real :: old_leaf_water
  real :: old_wood_water
  real :: old_leaf_water_im2
  real :: old_wood_water_im2
  real :: elim_nplant
  real :: elim_lai
  integer :: ifm,ipy,isi,ipa,ico,pft
  type(edtype), pointer :: cgrid
  type(polygontype), pointer :: cpoly
  type(met_driv_state), pointer :: cmet
  type(sitetype),pointer :: csite
  type(patchtype),pointer :: cpatch

  agb_frac  = real(agb_frac8)
  bgb_frac  = real(bgb_frac8)
  fol_frac  = real(fol_frac8)
  stor_frac = real(stor_frac8)

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
  print*,agb_frac,bgb_frac,fol_frac

  do ifm = 1,ngrids
     cgrid => edgrid_g(ifm)
     do ipy = 1,cgrid%npolygons

        cpoly => cgrid%polygon(ipy)

        do isi = 1,cpoly%nsites

           csite => cpoly%site(isi)
           cmet  => cpoly%met(isi)

           do ipa=1,csite%npatches

              cpatch => csite%patch(ipa)

              if(cpatch%ncohorts > 0) then

              do ico=1,cpatch%ncohorts

                 pft = cpatch%pft(ico)
                 !! calc new pool sizes

                 ialloc     =  1.0 / (1.0 + q(pft) + (qsw(pft)+qbark(pft))*cpatch%hite(ico))
                 bdeada_new = cpatch%bdeada(ico) * ( 1.0 - agb_frac)
                 bdeada_new = max(0.0, bdeada_new)
                 bdeadb_new = cpatch%bdeadb(ico) * ( 1.0 - bgb_frac)
                 bdeadb_new = max(0.0, bdeadb_new)
                 bsapa_new  = cpatch%balive(ico) * qsw(pft) * cpatch%hite(ico) * ialloc    &
                            * agf_bs(pft) * (1.0-agb_frac)
                 bsapa_new = max(0.0, bsapa_new)
                 bsapb_new  = cpatch%balive(ico) * qsw(pft) * cpatch%hite(ico) * ialloc    &
                            * (1.0 - agf_bs(pft)) * (1.0-bgb_frac)
                 bsapb_new = max(0.0, bsapb_new)
                 bbarka_new = cpatch%balive(ico) * qbark(pft) * cpatch%hite(ico) * ialloc  &
                            * agf_bs(pft) * (1.0-agb_frac)
                 bbarka_new = max(0.0, bbarka_new)
                 bbarkb_new = cpatch%balive(ico) * qbark(pft) * cpatch%hite(ico) * ialloc  &
                            * (1.0 - agf_bs(pft)) * (1.0-bgb_frac)
                 bbarkb_new = max(0.0, bbarkb_new)

                 bstore_new = cpatch%bstorage(ico) * (1.0-stor_frac)
                 bstore_new = max(0.0, bstore_new)
                 bleaf_new  = cpatch%balive(ico)   * ialloc *(1.0-fol_frac)
                 bfr_new    = cpatch%balive(ico)   * q(pft) * ialloc * (1.0-bgb_frac)
                 bfr_new    = max(0.0, bfr_new)

                 !! move residual frac to debris/litter pools
                   !! For now assume 100% removal [[needs to be updated]]

                 !! Adjust C balance to reflect new values
                 !! See https://github.com/EDmodel/ED2/issues/300
                 csite%cbudget_loss2yield(ipa) = csite%cbudget_loss2yield(ipa) &
                      + frqsumi * cpatch%nplant(ico) * &
                      ( (cpatch%broot(ico) - bfr_new) &
                      + (cpatch%bsapwooda(ico) - bsapa_new) &
                      + (cpatch%bsapwoodb(ico) - bsapb_new) &
                      + (cpatch%bbarka(ico) - bbarka_new)  &
                      + (cpatch%bbarkb(ico) - bbarkb_new) &
                      + (cpatch%bdeada(ico) - bdeada_new) &
                      + (cpatch%bdeadb(ico) - bdeadb_new)   &
                      + (cpatch%bstorage(ico) - bstore_new) )

                 !! update biomass pools
                 !! [[this needs to be more sophisticated]]

                 cpatch%broot(ico)     = bfr_new
                 cpatch%bsapwooda(ico) = bsapa_new
                 cpatch%bsapwoodb(ico) = bsapb_new
                 cpatch%bbarka(ico)    = bbarka_new
                 cpatch%bbarkb(ico)    = bbarkb_new
                 cpatch%bdeada(ico)    = bdeada_new
                 cpatch%bdeadb(ico)    = bdeadb_new
                 cpatch%bstorage(ico)  = bstore_new

                 cpatch%balive(ico)    = ed_balive(cpatch,ico)
                 
                 if(bleaf_new .le. tiny(1.0)) then
                    cpatch%phenology_status(ico) = -2
                    cpatch%elongf(ico)           = 0.0
                    !----- No leaves, then set it to zero. --------------------------------!
                    cpatch%bleaf(ico)      = 0.0
                 else
                    cpatch%phenology_status(ico) = 1
                    !! Again, adjust C balance accordingly.
                    csite%cbudget_loss2yield(ipa) = csite%cbudget_loss2yield(ipa) &
                         + (cpatch%bleaf(ico) - bleaf_new) * cpatch%nplant(ico) * frqsumi
                    cpatch%bleaf(ico)   = bleaf_new
                 end if

                 if((cpatch%bdeada(ico) + cpatch%bdeadb(ico)) > tiny(1.0)) then
                    if(is_grass(cpatch%pft(ico)) .and. igrass==1) then
                       cpatch%hite(ico) = max( hgt_min(pft),                               &
                                               bl2h(cpatch%bleaf(ico),cpatch%sla(ico),pft))
                       cpatch%dbh (ico) = h2dbh(cpatch%hite(ico),pft)
                    else
                       cpatch%dbh (ico) = bd2dbh(cpatch%pft(ico), cpatch%bdeada(ico)       &
                                                ,cpatch%bdeadb(ico))
                       cpatch%hite(ico) = dbh2h (cpatch%pft(ico), cpatch%dbh(ico))
                    end if
                 else
                    cpatch%dbh(ico)  = 0.0
                    cpatch%hite(ico) = 0.0
                 end if

                 !----- Update LAI, WAI, and CAI ------------------------------------------!
                 call area_indices(cpatch, ico)

                 !----- Update basal area and above-ground biomass. -----------------------!
                 cpatch%basarea(ico) = pio4 * cpatch%dbh(ico) * cpatch%dbh(ico)
                 cpatch%btimber(ico) = size2bt(cpatch%dbh(ico),cpatch%hite(ico)            &
                                              ,cpatch%bdeada(ico),cpatch%bsapwooda(ico)    &
                                              ,cpatch%bbarka(ico),cpatch%pft(ico))
                 cpatch%thbark (ico) = size2xb(cpatch%dbh(ico),cpatch%hite(ico)            &
                                              ,cpatch%bbarka(ico),cpatch%bbarkb(ico)       &
                                              ,cpatch%sla(ico),cpatch%pft(ico))
                 cpatch%agb(ico)     = ed_biomass(cpatch, ico)
                 !-------------------------------------------------------------------------!
                 !    Here we are leaving all water in the branches and twigs... Do not    !
                 ! worry, if there is any, it will go down through shedding the next       !
                 ! step.                                                                   !
                 !-------------------------------------------------------------------------!
                 old_leaf_hcap      = cpatch%leaf_hcap     (ico)
                 old_wood_hcap      = cpatch%wood_hcap     (ico)
                 old_leaf_water     = cpatch%leaf_water    (ico)
                 old_wood_water     = cpatch%wood_water    (ico)
                 old_leaf_water_im2 = cpatch%leaf_water_im2(ico)
                 old_wood_water_im2 = cpatch%wood_water_im2(ico)
                 call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdeada(ico)                   &
                                   ,cpatch%bsapwooda(ico),cpatch%bbarka(ico)               &
                                   ,cpatch%nplant(ico),cpatch%pft(ico)                     &
                                   ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico))
                 call rwc2tw(cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico),cpatch%bleaf(ico)   &
                            ,cpatch%bsapwooda(ico),cpatch%bsapwoodb(ico)                   &
                            ,cpatch%bdeada(ico),cpatch%bdeadb(ico),cpatch%broot(ico)       &
                            ,cpatch%dbh(ico),cpatch%pft(ico),cpatch%leaf_water_int(ico)    &
                            ,cpatch%wood_water_int(ico))
                 call twi2twe(cpatch%leaf_water_int(ico),cpatch%wood_water_int(ico)        &
                             ,cpatch%nplant(ico),cpatch%leaf_water_im2(ico)                &
                             ,cpatch%wood_water_im2(ico))
                 call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap     &
                                            ,old_leaf_water,old_wood_water                 &
                                            ,old_leaf_water_im2,old_wood_water_im2         &
                                            ,.true.,.false.)

                 !----- Update flags telling whether leaves and branches can be solved. ---!
                 call is_resolvable(csite,ipa,ico,.false.,.false.,'event_harvest')

              enddo

              !! remove small cohorts
              call terminate_cohorts(csite,ipa,cmet,.false.,elim_nplant,elim_lai)

              call update_patch_derived_props(csite,ipa,.true.)
            end if  !! check to make sure there ARE cohorts

           enddo
           ! Update site properties. ## THINK ABOUT WHAT TO SET FLAG##########
           call update_site_derived_props(cpoly,0,isi)
        end do
     end do

  end do

  print*,"done HARVEST"

end subroutine event_harvest


subroutine event_planting(pft,density8)
  use rk4_integ_utils, only : initialize_rk4patches
  use update_derived_utils, only : update_patch_thermo_props  &
                                 , update_patch_thermo_fmean  &
                                 , update_patch_derived_props &
                                 , update_site_derived_props
  use grid_coms, only : ngrids,nzg,nzs
  use ed_state_vars,only: edgrid_g, &
       edtype,polygontype,sitetype, &
       patchtype,allocate_patchtype,copy_patchtype,deallocate_patchtype, &
       filltab_alltypes
  use disturbance      , only : plant_patch
  use ed_type_init     , only : new_patch_sfc_props
  implicit none
  integer(kind=4),intent(in) :: pft
  real(kind=8),intent(in) :: density8
  real :: density
  integer :: ifm,ipy,isi,ipa
  type(edtype), pointer :: cgrid
  type(polygontype), pointer :: cpoly
  type(sitetype),pointer :: csite
  real, parameter :: planting_ht = 1.0 !multiple of hgt_min

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
              call update_patch_thermo_props(csite,ipa,ipa,nzg,nzs,cpoly%lsl(isi) &
                                            ,cpoly%ntext_soil(:,isi))
              call update_patch_thermo_fmean(csite,ipa,ipa,nzg,cpoly%lsl(isi)     &
                                            ,cpoly%ntext_soil(:,isi))
              call plant_patch(csite,ipa,nzg,pft,density,cpoly%ntext_soil(:,isi)  &
                              ,planting_ht,cpoly%lsl(isi))
              call update_patch_derived_props(csite,ipa,.true.)
              call new_patch_sfc_props(csite, ipa,nzg,nzs,cpoly%ntext_soil(:,isi))

           enddo

           ! Update site properties. ## THINK ABOUT WHAT TO SET FLAG##########
           call update_site_derived_props(cpoly,1,isi)

        enddo

     end do
  end do

  ! Re-allocate integration buffer
  call initialize_rk4patches(.false.)

  ! Reset hdf vars since number of cohorts changed mid-month
  call filltab_alltypes

end subroutine event_planting

subroutine event_fertilize(rval8)
  use update_derived_utils, only : update_patch_derived_props, update_site_derived_props
  use grid_coms, only : ngrids
  use ed_state_vars,only: edgrid_g, &
       edtype,polygontype,sitetype, &
       patchtype,allocate_patchtype,copy_patchtype,deallocate_patchtype
  real(kind=8),intent(in),dimension(5) :: rval8

  real :: nh4,no3,p,k,ca
  integer :: ifm,ipy,isi,ipa
  type(edtype), pointer :: cgrid
  type(polygontype), pointer :: cpoly
  type(sitetype),pointer :: csite
  type(patchtype),pointer :: cpatch

  nh4 = real(rval8(1))
  no3 = real(rval8(2))
  p   = real(rval8(3))
  k   = real(rval8(4))
  ca  = real(rval8(5))

  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"         FERTILIZE          "
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,""
  print*,rval8

  do ifm = 1,ngrids
     cgrid => edgrid_g(ifm)
     do ipy = 1,cgrid%npolygons

        cpoly => cgrid%polygon(ipy)

        do isi = 1,cpoly%nsites

           csite => cpoly%site(isi)

           do ipa=1,csite%npatches

              cpatch => csite%patch(ipa)

              csite%mineralized_soil_N(ipa) = max(0.0,csite%mineralized_soil_N(ipa) + nh4 + no3)

              !! update patch properties
              call update_patch_derived_props(csite,ipa,.true.)
           enddo

           ! Update site properties. ## THINK ABOUT WHAT TO SET FLAG##########
           call update_site_derived_props(cpoly,0,isi)
        end do
     end do

  end do

end subroutine event_fertilize

subroutine event_irrigate(rval8)
  use grid_coms, only : ngrids,nzg
  use ed_state_vars,only: edgrid_g, &
       edtype,polygontype,sitetype, &
       patchtype
  use therm_lib, only: uint2tl,tl2uint
  use consts_coms, only : wdns,wdnsi,t3ple
  implicit none
  real(kind=8),intent(in) :: rval8

  real :: iwater,ienergy,fliq,soil_temp
  integer :: ifm,ipy,isi,ipa,k
  type(edtype), pointer :: cgrid
  type(polygontype), pointer :: cpoly
  type(sitetype),pointer :: csite
  type(patchtype),pointer :: cpatch

  iwater = real(rval8)

  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"         IRRIGATE           "
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,""
  print*,rval8

  if(iwater .lt. 0.0) then
     print*,"INVALID IRRIGATION VALUE"
     stop
  end if

  do ifm = 1,ngrids
     cgrid => edgrid_g(ifm)
     do ipy = 1,cgrid%npolygons

        cpoly => cgrid%polygon(ipy)

        do isi = 1,cpoly%nsites

           csite => cpoly%site(isi)

           do ipa=1,csite%npatches

              cpatch => csite%patch(ipa)

              !! note, assume irrigation water is at same temperature as soil
              if (csite%soil_tempk(nzg,ipa) == t3ple)then
                 fliq = 0.5
              elseif (csite%soil_tempk(nzg,ipa) > t3ple)then
                 fliq = 1.0
              else
                 fliq = 0.0
              end if
              ienergy = wdns * tl2uint(csite%soil_tempk(nzg,ipa),fliq)

              k = csite%nlev_sfcwater(ipa)
              if(k .eq. 0) then
                 csite%sfcwater_mass(1,ipa)    = iwater
                 csite%sfcwater_tempk(1,ipa)   = soil_temp
                 csite%sfcwater_energy(1,ipa)  = ienergy
                 csite%sfcwater_depth(1,ipa)   = iwater * wdnsi
                 csite%sfcwater_fracliq(1,ipa) = fliq
                 csite%nlev_sfcwater(ipa)      = 1
              else
                 csite%sfcwater_energy(k,ipa) = (csite%sfcwater_energy(k,ipa)*csite%sfcwater_mass(k,ipa) &
                      + ienergy*iwater)/(csite%sfcwater_mass(k,ipa) + iwater)
                 csite%sfcwater_mass(k,ipa)   = csite%sfcwater_mass(k,ipa) + iwater
                 csite%sfcwater_depth(k,ipa)  = csite%sfcwater_depth(k,ipa) + iwater*wdnsi
                 call uint2tl(csite%sfcwater_energy(k,ipa),csite%sfcwater_tempk(k,ipa),csite%sfcwater_fracliq(k,ipa))
              endif

              !! do we need to call infiltration?
              !! redistribute_snow takes a rk4 patch, not a normal patch

           enddo

        end do
     end do

  end do

end subroutine

subroutine event_fire()
  implicit none
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"         FIRE           "
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,""
  print*,"<<---- CURRENTLY UNIMPLEMENTED ---->>"


end subroutine event_fire

subroutine event_till(rval8)
  use update_derived_utils, only : update_patch_derived_props, update_site_derived_props
  use grid_coms, only : ngrids
  use ed_state_vars,only: edgrid_g, &
       edtype,polygontype,sitetype, &
       patchtype,allocate_patchtype,copy_patchtype,deallocate_patchtype
  use met_driver_coms, only : met_driv_state             ! ! structure
  use pft_coms, only: c2n_storage,c2n_leaf,c2n_stem,l2n_stem,f_labile_leaf,f_labile_stem,agf_bs
  use fuse_fiss_utils, only: terminate_cohorts
  use ed_therm_lib, only : calc_veg_hcap
  use therm_lib        , only : cmtl2uext
  use plant_hydro,     only : rwc2tw, twi2twe       ! ! sub-routine
  use consts_coms,     only : t3ple ! ! intent(in)
  implicit none
  real(kind=8),intent(in) :: rval8

  real :: depth
  integer :: ifm,ipy,isi,ipa,pft,ico
  real :: elim_nplant,elim_lai,a_bfast,b_bfast,a_bstruct,b_bstruct,a_bstorage,b_bstorage
  type(edtype), pointer :: cgrid
  type(polygontype), pointer :: cpoly
  type(sitetype),pointer :: csite
  type(met_driv_state), pointer :: cmet
  type(patchtype),pointer :: cpatch


  depth = real(rval8)


  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"         TILLING            "
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,"----------------------------"
  print*,""
  print*,rval8

  do ifm = 1,ngrids
     cgrid => edgrid_g(ifm)
     do ipy = 1,cgrid%npolygons

        cpoly => cgrid%polygon(ipy)

        do isi = 1,cpoly%nsites

           csite => cpoly%site(isi)
           cmet  => cpoly%met(isi)

           do ipa=1,csite%npatches

              cpatch => csite%patch(ipa)
              print*,cpatch%ncohorts
              if(cpatch%ncohorts > 0) then

              do ico=1,cpatch%ncohorts

                 pft = cpatch%pft(ico)
                 a_bfast  = f_labile_leaf(pft)           * cpatch%bleaf (ico)              &
                          + f_labile_stem(pft)                                             &
                          * ( cpatch%bsapwooda(ico)      + cpatch%bbarka(ico)              &
                            + cpatch%bdeada(ico)                             )
                 b_bfast  = f_labile_leaf(pft)           * cpatch%broot (ico)              &
                          + f_labile_stem(pft)                                             &
                          * ( cpatch%bsapwoodb(ico)      + cpatch%bbarkb(ico)              &
                            + cpatch%bdeadb(ico)                             )
                 a_bstruct  = (1.0 - f_labile_leaf(pft)) * cpatch%bleaf (ico)              &
                            + (1.0 - f_labile_stem(pft))                                   &
                            * ( cpatch%bsapwooda(ico)    + cpatch%bbarka(ico)              &
                              + cpatch%bdeada(ico)                           )
                 b_bstruct  = (1.0 - f_labile_leaf(pft)) * cpatch%broot (ico)              &
                            + (1.0 - f_labile_stem(pft))                                   &
                            * ( cpatch%bsapwoodb(ico)    + cpatch%bbarkb(ico)              &
                              + cpatch%bdeadb(ico)                           )
                 a_bstorage = agf_bs(pft) * cpatch%bstorage(ico)
                 b_bstorage = (1.0 - agf_bs(pft)) * cpatch%bstorage(ico)
                 !! MLO - I multiplied bfast/bstruct/bstorage by cpatch%nplant.  
                 !!       I never used this routine and I am not sure about this fix,
                 !!       but it seems to me that the units were inconsistent
                 !!       (fast_soil_C/struct_soil_C are in kgC/m2 whereas
                 !!       bfast/bstruct/bstorage are in kgC/plant)
                 
                 !! move biomass to debris/litter pools
                 csite%fast_grnd_C(ipa)       = csite%fast_grnd_C(ipa)                     &
                                              + cpatch%nplant(ico) * (a_bfast + a_bstorage)
                 csite%fast_soil_C(ipa)       = csite%fast_soil_C(ipa)                     &
                                              + cpatch%nplant(ico) * (b_bfast + b_bstorage)

                 csite%structural_grnd_C(ipa) = csite%structural_grnd_C(ipa)               &
                                              + cpatch%nplant(ico) * a_bstruct
                 csite%structural_soil_C(ipa) = csite%structural_soil_C(ipa)               &
                                              + cpatch%nplant(ico) * b_bstruct
                 csite%structural_grnd_L(ipa) = csite%structural_grnd_L(ipa)               &
                                              + cpatch%nplant(ico) * a_bstruct             &
                                              * l2n_stem / c2n_stem(pft)
                 csite%structural_soil_L(ipa) = csite%structural_soil_L(ipa)               &
                                              + cpatch%nplant(ico) * b_bstruct             &
                                              * l2n_stem / c2n_stem(pft)
                 csite%fast_grnd_N(ipa) = csite%fast_grnd_N(ipa)                           &
                                        + cpatch%nplant(ico) * ( a_bfast/c2n_leaf(pft)     &
                                                               + a_bstorage / c2n_storage )
                 csite%fast_soil_N(ipa) = csite%fast_soil_N(ipa)                           &
                                        + cpatch%nplant(ico) * ( b_bfast/c2n_leaf(pft)     &
                                                               + b_bstorage / c2n_storage )


                 !! update biomass pools
                 cpatch%balive(ico)           = 0.0
                 cpatch%broot(ico)            = 0.0
                 cpatch%bsapwooda(ico)        = 0.0
                 cpatch%bsapwoodb(ico)        = 0.0
                 cpatch%bbarka(ico)           = 0.0
                 cpatch%bbarkb(ico)           = 0.0
                 cpatch%bdeada(ico)           = 0.0
                 cpatch%bdeadb(ico)           = 0.0
                 cpatch%bstorage(ico)         = 0.0
                 cpatch%thbark(ico)           = 0.0
                 cpatch%nplant(ico)           = 0.0
                 cpatch%lai(ico)              = 0.0
                 cpatch%wai(ico)              = 0.0
                 cpatch%crown_area(ico)       = 0.0
                 cpatch%bleaf(ico)            = 0.0
                 cpatch%leaf_water(ico)       = 0.0
                 cpatch%leaf_temp(ico)        = csite%can_temp(ipa)
                 cpatch%leaf_temp_pv(ico)     = csite%can_temp(ipa)
                 cpatch%leaf_vpdef(ico)       = csite%can_vpdef(ipa)
                 cpatch%wood_water(ico)       = 0.0
                 cpatch%wood_temp(ico)        = csite%can_temp(ipa)
                 cpatch%wood_temp_pv(ico)     = csite%can_temp(ipa)


                 if (cpatch%leaf_temp(ico) == t3ple) then
                    cpatch%leaf_fliq(ico) = 0.5
                 elseif (cpatch%leaf_temp(ico) > t3ple) then
                    cpatch%leaf_fliq(ico) = 1.0
                 else
                    cpatch%leaf_fliq(ico) = 0.0
                 end if
                 if (cpatch%wood_temp(ico) == t3ple) then
                    cpatch%wood_fliq(ico) = 0.5
                 elseif (cpatch%wood_temp(ico) > t3ple) then
                    cpatch%wood_fliq(ico) = 1.0
                 else
                    cpatch%wood_fliq(ico) = 0.0
                 end if
                 call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdeada(ico)                   &
                                   ,cpatch%bsapwooda(ico),cpatch%bbarka(ico)               &
                                   ,cpatch%nplant(ico),cpatch%pft(ico)                     &
                                   ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico))
                 call rwc2tw(cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico)                     &
                            ,cpatch%bleaf(ico),cpatch%bsapwooda(ico),cpatch%bsapwoodb(ico) &
                            ,cpatch%bdeada(ico),cpatch%bdeadb(ico),cpatch%broot(ico)       &
                            ,cpatch%dbh(ico),cpatch%pft(ico),cpatch%leaf_water_int(ico)    &
                            ,cpatch%wood_water_int(ico))
                 call twi2twe(cpatch%leaf_water_int(ico),cpatch%wood_water_int(ico)        &
                             ,cpatch%nplant(ico),cpatch%leaf_water_im2(ico)                &
                             ,cpatch%wood_water_im2(ico))
                 cpatch%leaf_energy(ico) = cmtl2uext( cpatch%leaf_hcap     (ico)           &
                                                    , cpatch%leaf_water    (ico)           &
                                                    + cpatch%leaf_water_im2(ico)           &
                                                    , cpatch%leaf_temp     (ico)           &
                                                    , cpatch%leaf_fliq     (ico) )
                 cpatch%wood_energy(ico) = cmtl2uext( cpatch%wood_hcap     (ico)           &
                                                    , cpatch%wood_water    (ico)           &
                                                    + cpatch%wood_water_im2(ico)           &
                                                    , cpatch%wood_temp     (ico)           &
                                                    , cpatch%wood_fliq     (ico) )
                 cpatch%phenology_status(ico) = -2
                 cpatch%elongf(ico)           = 0.0

              enddo
              !! remove small cohorts
              call terminate_cohorts(csite,ipa,cmet,.false.,elim_nplant,elim_lai)

              !! update patch properties
              call update_patch_derived_props(csite,ipa,.true.)
              endif
           enddo
           ! Update site properties. ## THINK ABOUT WHAT TO SET FLAG##########
           call update_site_derived_props(cpoly,0,isi)
        end do
     end do

  end do

  print*,"end TILL"

end subroutine event_till


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
!!$              cpatch%bdead(ico)   = size2bd(cpatch%dbh(ico),cpatch%hite(ico),pft)
!!$              cpatch%bleaf(ico)   = size2bl(cpatch%dbh(ico),cpatch%hite(ico),cpatch%sla(ico),pft)
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
!!$              cpatch%leaf_energy(ico) =  hcapveg * (cpatch%leaf_temp(ico)-t3ple)
!!$print*,cpatch%veg_energy(ico)
!!$
!!$              !! count as new recruitment?
!!$              cpatch%new_recruit_flag(ico) = 1
!!$print*,"init_cohorts"
!!$              call init_ed_cohort_vars(cpatch,ico,cpoly%lsl(isi))
!!$
!!$              csite%cohort_count(ipa) = csite%cohort_count(ipa) + 1
!!$              cpatch%ncohorts = ico
!!$              csite%paco_n(ipa) = ico
!!$
!!$              ! Sort the cohorts so that the new cohort is at the correct height bin
!!$
!!$              call sort_cohorts(cpatch)
!!$
!!$!              call update_patch_derived_props(csite, ipa,.true.)
!!$
