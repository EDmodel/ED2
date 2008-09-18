subroutine read_ed1_history_file_array


  use max_dims, only: n_pft,huge_patch,huge_cohort,max_water,str_len
  use pft_coms, only: SLA, q, qsw, hgt_min, include_pft, include_pft_ag, phenology,pft_1st_check,include_these_pft
  use misc_coms, only: sfilin, ied_init_mode
  use mem_sites, only: grid_res,edres
  use consts_coms, only: pio180
  use ed_misc_coms, only: use_target_year, restart_target_year
  use ed_state_vars,only:polygontype,sitetype,patchtype,edtype, &
       edgrid_g,allocate_sitetype,allocate_patchtype
  use grid_coms,only:ngrids

  implicit none

  integer :: year
  real, external :: dbh2h

  type(edtype),pointer      :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer    :: csite
  type(patchtype),pointer   :: cpatch

  real :: h2dbh
  real :: dbh2bd
  real :: dbh2bl
  
  integer :: ii,pft
  character(len=str_len) :: ed_fname
  real :: flat,lat
  real :: flon,lon
  integer :: ilat,ilon
  
  logical, parameter :: renumber_pfts = .true.
  character(len=str_len) :: pss_name
  character(len=str_len) :: css_name
  character(len=str_len) :: site_name
  real :: dist
  real :: best_dist

  character(len=str_len) :: best_pss_name
  character(len=str_len) :: best_css_name
  logical :: restart_exist

  integer :: nsites
  integer :: igr,ipy,ipy2,isi,ipa,ico
  integer :: ip,ip2,ic,ic2

  character(len=str_len) :: cdum
  integer :: nwater
  real , dimension(max_water) :: depth
  integer, dimension(n_pft) :: include_pft_ep
  ! Patch and site variables from the restart

  integer,dimension(huge_patch) :: sitenum
  real ,dimension(huge_patch) :: time
  character(len=str_len) ,dimension(huge_patch) :: pname
  integer, dimension(huge_patch) :: trk
  real, dimension(huge_patch) :: age
  real, dimension(huge_patch) :: area
  real, dimension(huge_patch) :: fsc
  real, dimension(huge_patch) :: stsc
  real, dimension(huge_patch) :: stsl
  real, dimension(huge_patch) :: ssc
  real, dimension(huge_patch) :: psc
  real, dimension(huge_patch) :: msn
  real, dimension(huge_patch) :: fsn
  real, dimension(max_water,huge_patch) :: water

  !    All variables became double precision. ED-1 apparently had it, so some values are often smaller 
  ! than 1.18E-38, which causes FPE problems.
  real(kind=8)            :: dage,darea,dfsc,dstsc,dstsl,dssc,dpsc,dmsn,dfsn
  real(kind=8), dimension(max_water) :: dwater
  
  ! This is the smallest representable number in single precision
  real(kind=8), parameter :: snglmin=dble(tiny(1.))
  real, parameter :: min_area=epsilon(1.) ! Doesn't need to be the machine epsilon, chose just a small number.
  real(kind=8), parameter :: min_ok=dble(tiny(1.)/epsilon(1.)) ! Chose a number small enough, but with some
                                                               ! room for multiplying by a small area.
  ! Cohort variables from the restart

  integer :: ierr
  real, dimension(12,huge_cohort) :: cb
  real, dimension(12,huge_cohort) :: cb_max
  integer, dimension(huge_cohort) :: leaves_on
  real, dimension(huge_cohort) :: balive
  real, dimension(huge_cohort) :: avgRg
  real, dimension(huge_cohort) :: bdead
  real, dimension(huge_cohort) :: nplant
  real, dimension(huge_cohort) :: hite
  real, dimension(huge_cohort) :: dbh
  integer, dimension(huge_cohort) :: ipft
  character(len=64), dimension(huge_cohort) :: cname,cpname
  real ,dimension(huge_cohort) :: ctime
  logical, dimension(huge_cohort) :: add_this_cohort

  real :: area_tot
  real :: area_sum
  real :: patch_lai,poly_lai
  real :: site_lai
  integer :: ncohorts,npatchco
  integer :: npatches,nsitepat,npatch2
  integer, parameter :: harvard_override = 0
  logical :: site_match
  real :: dist_gc
  integer :: nw
  
  real, external :: ed_biomass
  
  
  ! Loop over all grids, polygons, and sites

  do igr = 1,ngrids

     cgrid => edgrid_g(igr)

     do ipy = 1,cgrid%npolygons
        
        cpoly => cgrid%polygon(ipy)

        ! Override of LSL, ntext_soil
        if(harvard_override == 1)then

           cpoly%lsl(1)   = 1
           cgrid%lsl(ipy) = 1
           cgrid%ntext_soil(1,ipy) = 2
           cgrid%ntext_soil(2,ipy) = 2
           cgrid%ntext_soil(3,ipy) = 3
           cgrid%ntext_soil(4,ipy) = 3
        endif

        ! =================================
        ! Part I: Find the restart files
        ! =================================

        call create_ed1_fname(cgrid%lat(ipy), edres, cgrid%lon(ipy),  &
             trim(sfilin), pss_name, css_name, site_name)
        
        inquire(file=trim(pss_name),exist=restart_exist)
        if(.not.restart_exist)then
           write(unit=*,fmt='(a)') ' + Your restart file does not exist:'
           write(unit=*,fmt='(a)') '    - Not found:'//trim(pss_name)
           ! Find nearest neighbor.
           best_dist = huge(6.)
           best_pss_name = 'null'

           ! IT IS POSSIBLE THAT NONE OF THE POLYGONS ON THIS TILE
           ! ARE CLOSE TO VALUES IN THE DATABASE
           do ipy2 = 1,cgrid%npolygons
              
              call create_ed1_fname(cgrid%lat(ipy2), edres, cgrid%lon(ipy2),  &
                   trim(sfilin), pss_name, css_name, site_name)
              inquire(file=trim(pss_name),exist=restart_exist)
              if(restart_exist)then
                 
                 dist = dist_gc(cgrid%lon(ipy),cgrid%lon(ipy2),cgrid%lat(ipy),cgrid%lat(ipy2))
                 if(dist < best_dist)then
                    best_dist = dist
                    best_pss_name = trim(pss_name)
                    best_css_name = trim(css_name)
                 endif
              endif
              
           enddo

           
           if(trim(best_pss_name) /= 'null')then
              pss_name = trim(best_pss_name)
              css_name = trim(best_css_name)
              write(unit=*,fmt='(a)') '    - Using:'//trim(pss_name)//' instead.'
           else

              ! LAST RESORT - USE SPARINGLY
              do ilat = 1,180
                 do ilon = 1,360
                    lat = -90.0 + real(ilat-1)
                    lon = -180.0 + real(ilon-1)
                    call create_ed1_fname(lat, edres,lon,  &
                         trim(sfilin), pss_name, css_name, site_name)
                    inquire(file=trim(pss_name),exist=restart_exist)
                    if(restart_exist)then
                       
                       dist = dist_gc(cgrid%lon(ipy),cgrid%lon(ipy2),cgrid%lat(ipy),cgrid%lat(ipy2))
                       if(dist < best_dist)then
                          best_dist = dist
                          best_pss_name = trim(pss_name)
                          best_css_name = trim(css_name)
                       endif
                    endif
                 enddo
              enddo
              if(trim(best_pss_name) /= 'null')then
                 pss_name = trim(best_pss_name)
                 css_name = trim(best_css_name)
                 write(unit=*,fmt='(a)') '    - Using:'//trim(pss_name)//' instead.'
              else
                 
                 call fatal_error('Cannot find a suitable restart file.' &
                      ,'read_ed1_history_file','ed_history_io_array.f90')
              endif
           end if
        else
           !           write(unit=*,fmt='(a)') ' + Your restart file exists:'//trim(pss_name)
           
        endif

        ! =================================
        ! Part II: Add the patches
        ! =================================


        ! Open file and read in patches
        open(12,file=trim(pss_name),form='formatted',status='old')
        read(12,*)  ! skip header
        if(ied_init_mode == 1) then
           read(12,*)cdum,nwater
           read(12,*)cdum,depth(1:nwater)
           read(12,*)
        else
           nwater = 1
           read(12,*)!water patch
        endif

        ! Note that if we are doing an ED1 restart we can 
        ! assume 1 site?
        
        ip = 1
        sitenum = 0
        count_patches: do
           
           if (ip>huge_patch) call fatal_error('IP too high','read_ed1_history_file_array','ed_history_io.f90')

           select case (ied_init_mode)
           case (1) !! read ED1 format files
              
              read(12,*,iostat=ierr)  time(ip),pname(ip),trk(ip),dage,darea,dfsc,dstsc,  &
                                      dstsl,dssc,dpsc,dmsn,dfsn,dwater(1:nwater)
              if(ierr /= 0)exit count_patches
              area(ip)   = sngl(max(snglmin,darea  ))
              age(ip)    = sngl(max(min_ok ,dage   ))
              fsc(ip)    = sngl(max(min_ok ,dfsc   ))
              stsc(ip)   = sngl(max(min_ok ,dstsc  ))
              stsl(ip)   = sngl(max(min_ok ,dstsl  ))
              ssc(ip)    = sngl(max(min_ok ,dssc   ))
              psc(ip)    = sngl(max(min_ok ,dpsc   ))
              msn(ip)    = sngl(max(min_ok ,dmsn   ))
              fsn(ip)    = sngl(max(min_ok ,dfsn   ))
              do nw=1,nwater
                 water(nw,ip)  = sngl(max(min_ok ,dwater(nw) ))
              end do
              
           case(2)  !! read ED2 format files
              
              read(12,*,iostat=ierr)time(ip),pname(ip),trk(ip),dwater(1),dage,darea,dfsc,dstsc  &
                                   ,dstsl,dssc,dpsc,dmsn,dfsn
              if(ierr /= 0)exit count_patches
              area(ip)   = sngl(max(snglmin,darea  ))
              age(ip)    = sngl(max(min_ok ,dage   ))
              fsc(ip)    = sngl(max(min_ok ,dfsc   ))
              stsc(ip)   = sngl(max(min_ok ,dstsc  ))
              stsl(ip)   = sngl(max(min_ok ,dstsl  ))
              ssc(ip)    = sngl(max(min_ok ,dssc   ))
              psc(ip)    = sngl(max(min_ok ,dpsc   ))
              msn(ip)    = sngl(max(min_ok ,dmsn   ))
              fsn(ip)    = sngl(max(min_ok ,dfsn   ))
              water(1,ip)  = sngl(max(min_ok ,dwater(1) ))
              
           case(3)
              
              read(12,*,iostat=ierr) sitenum(ip),time(ip),pname(ip),trk(ip),water(1,ip),age(ip), &
                   darea,fsc(ip),stsc(ip),stsl(ip),ssc(ip),psc(ip),msn(ip),fsn(ip)
              if(ierr /= 0)exit count_patches
              area(ip)=sngl(max(snglmin,darea))
              
              if(sitenum(ip)<= 0) continue !! check for valid site number
              
              !! check for valid year
              year = int(time(ip))
              if(use_target_year.eq.1 .and. year .ne. restart_target_year) continue
              
              site_match = .false.
              do isi = 1,cpoly%nsites
                 if(sitenum(ip).eq.cpoly%sitenum(isi)) then
                    site_match=.true.
                 endif
              enddo
              if(.not.site_match) then
                 print*,"error reading from patch file",pss_name
                 print*,"site number", sitenum,"not found"
                 stop
              endif
           case default !Nearly bare ground
              exit count_patches
           end select
           
           if (area(ip) > min_area) ip = ip + 1 ! This will remove patches with tiny area that often cause trouble
        
        enddo count_patches

        npatches = max(ip-1,0)
        ! Allocate memory
        
        if(ied_init_mode == 3) then
           
           ! If this is type 3, then we need to find the number of patches per site
           do isi = 1,cpoly%nsites
              
              npatch2 = 0
              do ip=1,npatches
                 if(sitenum(ip).eq.cpoly%sitenum(isi)) then
                    npatch2 = npatch2 + 1
                 endif
              enddo
              
              csite => cpoly%site(isi)
              
              ! Allocate the patches in this site
              call allocate_sitetype(csite,npatch2)
              
              ! Go through the patch list and fill the sites'
              ! patches with data
              
              ip2 = 0
              
              do ip=1,npatches
                 
                 if(sitenum(ip).eq.cpoly%sitenum(isi)) then
                    ip2 = ip2 + 1
                    csite%dist_type(ip2)          = trk(ip) + 1
                    csite%age(ip2)                = age(ip)
                    csite%area(ip2)               = area(ip)
                    csite%fast_soil_C(ip2)        = fsc(ip)
                    csite%slow_soil_C(ip2)        = ssc(ip)
                    csite%structural_soil_C(ip2)  = stsc(ip)
                    csite%structural_soil_L(ip2)  = stsl(ip)
                    csite%mineralized_soil_N(ip2) = msn(ip)
                    csite%fast_soil_N(ip2)        = fsn(ip)
                    csite%pname(ip2)              = trim(pname(ip))
                    csite%sum_dgd(ip2)            = 0.0
                    csite%sum_chd(ip2)            = 0.0
                    csite%plantation(ip2)         = 0
                    csite%cohort_count(ip2)       = 0
                 endif
              enddo
              ! Initialize the cohort counts per patch
              csite%cohort_count(:) = 0
           enddo
           
        else
           
           ! We can assume 1 site per polygon, so lets point to it and
           ! fill it up
           
           csite => cpoly%site(1)
           
           ! Allocate the patches in this site
           call allocate_sitetype(csite,npatches)
           
           do ip=1,npatches
              csite%dist_type(ip)          = trk(ip) + 1
              csite%age(ip)                = age(ip)
              csite%area(ip)               = area(ip)
              csite%fast_soil_C(ip)        = fsc(ip)
              csite%slow_soil_C(ip)        = ssc(ip)
              csite%structural_soil_C(ip)  = stsc(ip)
              csite%structural_soil_L(ip)  = stsl(ip)
              csite%mineralized_soil_N(ip) = msn(ip)
              csite%fast_soil_N(ip)        = fsn(ip)
              csite%pname(ip)              = trim(pname(ip))
              csite%sum_dgd(ip)            = 0.0+tiny(1.)
              csite%sum_chd(ip)            = 0.0+tiny(1.)
              csite%plantation(ip)         = 0.0+tiny(1.)
              csite%cohort_count(ip)       = 0.0+tiny(1.)
           enddo
           
           
           ! Initialize the cohort counts per patch
           csite%cohort_count(:) = 0
           
        endif
        
        close(12)

        ! =================================
        ! Part III: Add the cohorts
        ! =================================
        
        
        open(12,file=trim(css_name),form='formatted',status='old')
        read(12,*)  ! skip header
        read(12,*)  ! skip header
        
        ic = 0
        
        read_cohorts: do

           ic = ic + 1
           
           if (ic>huge_cohort) call fatal_error('IC too high','read_ed1_history_file_array','ed_history_io.f90')

           select case (ied_init_mode)
           case (1)
              read(12,*,iostat=ierr)ctime(ic),cpname(ic),cname(ic),dbh(ic),hite(ic),ipft(ic),nplant(ic),  &
                   bdead(ic),balive(ic),avgRg(ic),leaves_on(ic),cb(1:12,ic),cb_max(1:12,ic)
              if(ierr /= 0)exit read_cohorts  
           case (2,3)
              read(12,*,iostat=ierr)ctime(ic),cpname(ic),cname(ic),dbh(ic),hite(ic),ipft(ic),nplant(ic),  &
                   bdead(ic),balive(ic),avgRg(ic)
              if(ierr /= 0)exit read_cohorts
              cb(1:12,ic) = 1.0
              cb_max(1:12,ic) = 1.0           
           end select

           if(renumber_pfts) then
              if(ipft(ic) < 100)then
                 ipft(ic) = ipft(ic) + 1
                 if(ipft(ic) >= 5)ipft(ic) = ipft(ic) - 3
              else
                 ipft(ic) = ipft(ic) - 100
              endif
           endif
                     
           !! check if the year matches
           year = int(ctime(ic))
           if(use_target_year == 1 .and. year .ne. restart_target_year) continue
           
           ! Find site and patch and start counting how many to allocate
           
           put_cohort:do isi=1,cpoly%nsites
              csite => cpoly%site(isi)
              do ipa=1,csite%npatches
                 if (include_pft(ipft(ic)) == 0) then
                    select case (pft_1st_check)
                    case(0)
                       write (unit=*,fmt='(a,1x,i5,1x,a)') &
                            'I found a cohort with PFT=',ipft(ic),' and it is not in your include_these_pft...'
                       call fatal_error('Invalid PFT in history file','read_ed1_history_file_array','ed_history_io.f90')
                    case(1)
                       write (unit=*,fmt='(a,1x,i5,1x,a)') &
                            'I found a cohort with PFT=',ipft(ic),'... Including this PFT in your include_these_pft...'
                       include_pft(ipft(ic)) = 1
                       include_these_pft(sum(include_pft)) = ipft(ic)
                       call sort_up(include_these_pft,n_pft)
                       if (ipft(ic) == 1 .or. ipft(ic) == 5) include_pft_ag(ipft(ic)) = 1
                       add_this_cohort(ic) = .true.
                    case(2)
                       write (unit=*,fmt='(a,1x,i5,1x,a)') &
                            'I found a cohort with PFT=',ipft(ic),'... Ignoring it...'
                       add_this_cohort(ic) = .false.
                    end select
                 else
                    add_this_cohort(ic)=.true.
                 end if
                 
                 if(trim(csite%pname(ipa)) == trim(cpname(ic) ) .and. add_this_cohort(ic)) then
                    csite%cohort_count(ipa) = csite%cohort_count(ipa) + 1
                    exit put_cohort
                 endif
              enddo
           enddo put_cohort
           
        enddo read_cohorts
        
        ncohorts = max(ic-1,0)
        
        close(12)

        loop_sites: do isi=1,cpoly%nsites


           csite => cpoly%site(isi)
           loop_patches: do ipa=1,csite%npatches

              cpatch => csite%patch(ipa)
              
              if (csite%cohort_count(ipa) /= 0) then
                 call allocate_patchtype(cpatch,csite%cohort_count(ipa))
                 csite%plant_ag_biomass(ipa) = 0.
                 ic2 = 0
                 do ic = 1,ncohorts
                    
                    if(trim(csite%pname(ipa)) == trim(cpname(ic) ) .and. add_this_cohort(ic)) then
                       ic2 = ic2 + 1
                       
                       cpatch%pft(ic2) = ipft(ic)
                       if(ied_init_mode == 1) then
                          cpatch%nplant(ic2) = nplant(ic) / (csite%area(ipa) )
                       else
                          cpatch%nplant(ic2) = nplant(ic)
                       endif
                       
                       cpatch%dbh(ic2) = dbh(ic)
                       
                       
                       if(hite(ic) > 0.0)then
                          cpatch%hite(ic2) = hite(ic)
                       else
                          cpatch%hite(ic2) = dbh2h(ipft(ic), dbh(ic))
                       endif
                       
                       if(bdead(ic) > 0.0)then
                          cpatch%bdead(ic2) = bdead(ic)
                       else
                          cpatch%bdead(ic2) = dbh2bd(dbh(ic), cpatch%hite(ic2), ipft(ic))
                       endif

                       ! Setting balive to default instead of file value
                       
                       cpatch%bleaf(ic2) = dbh2bl(cpatch%dbh(ic2),ipft(ic))
                       
                       cpatch%balive(ic2) = cpatch%bleaf(ic2) * (1.0 + q(ipft(ic)) +  &
                            qsw(ipft(ic)) * cpatch%hite(ic2))
                       
                       cpatch%lai(ic2) = cpatch%bleaf(ic2) * cpatch%nplant(ic2) *   &
                            SLA(ipft(ic))

!                       print*,cpatch%lai(ic2),cpatch%bleaf(ic2),cpatch%nplant(ic2),SLA(ipft(ic)),ipft(ic)
                       
                       ! START COLD-DECIDUOUS TREES WITHOUT LEAVES.  ALL OTHER TREES
                       ! ARE FULLY FLUSHED.
                       if(phenology(ipft(ic)) /= 2)then
                          cpatch%phenology_status(ic2) = 0
                          cpatch%bstorage(ic2) = 0.0

                       else
                          cpatch%phenology_status(ic2) = 2
                          cpatch%bleaf(ic2) = 0.0
                          cpatch%lai(ic2) = 0.0
                          cpatch%bstorage(ic2) = 0.5 * cpatch%balive(ic2)
                       endif
                                  
                       cpatch%cb(1:12,ic2) = cb(1:12,ic)
                       cpatch%cb_max(1:12,ic2) = cb_max(1:12,ic)
                       cpatch%cb(13,ic2) = 0.0
                       cpatch%cb_max(13,ic2) = 0.0
                       
                       ! Initialize other cohort variables. Some of them won't be updated
                       ! unless the lai goes above lai_min
                       cpatch%fsw(ic2)   = 1.0
                       cpatch%gpp(ic2)   = 0.0
                       cpatch%par_v(ic2) = 0.0
                       
                       csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa) +         &
                           ed_biomass(cpatch%bdead(ic2),cpatch%balive(ic2)                 &
                                     ,cpatch%bleaf(ic2),cpatch%pft(ic2),cpatch%hite(ic2)   &
                                     ,cpatch%bstorage(ic2))* cpatch%nplant(ic2)
                    endif
                 enddo
              else ! if (csite%cohort_count(ipa) == 0) then

           ! MLO 5-27-08. Force the patch to have one cohort of each pft that should be included
           !              Considers whether this is agricultural or forest.
              
                 if (csite%dist_type(ipa) == 1) then
                   include_pft_ep = include_pft_ag
                 else 
                   include_pft_ep = include_pft
                end if
                ic = sum(include_pft_ep)
                 ! MLO - 5-27-08. "Phylosophical" question. If the patch has no cohort, shouldn't we 
                 !                reset its age to zero? It will behave as a near-bare ground patch...
                 !                I just set it up to zero here, if this is wrong please remove it...
                 csite%age(ipa) = 0.
                 ! MLO - 5-27-08. Another "phylosophical" question. If the patch has no cohort and 
                 !                it is not water, should it even exist?
                 ! Initialize aboveground biomass for this site.
                 csite%plant_ag_biomass(ipa) = 0.
                 call allocate_patchtype(cpatch,ic)
                 csite%cohort_count(ipa) = ic
                 ic = 0
                 do pft = 1,n_pft
                    if(include_pft_ep(pft) == 1)then
                       
                       ic = ic + 1
                       
                       ! Define the near-bare ground
                       cpatch%pft(ic)     = pft
                       cpatch%hite(ic)    = hgt_min(pft)
                       cpatch%dbh(ic)     = h2dbh(cpatch%hite(ic),pft)
                       cpatch%bdead(ic)   = dbh2bd(cpatch%dbh(ic),cpatch%hite(ic),pft)
                       cpatch%bleaf(ic)   = dbh2bl(cpatch%dbh(ic),pft)
                       cpatch%nplant(ic)  = 0.1
                       cpatch%phenology_status(ic) = 0
                       cpatch%balive(ic)  = cpatch%bleaf(ic) * (1.0 + q(pft) +  &
                            qsw(pft) * cpatch%hite(ic))
                       cpatch%lai(ic)      = cpatch%bleaf(ic) * cpatch%nplant(ic) * SLA(pft)
                       cpatch%bstorage(ic) = 0.0
                       csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa) +         &
                           ed_biomass(cpatch%bdead(ic),cpatch%balive(ic), cpatch%bleaf(ic) &
                                     ,cpatch%pft(ic), cpatch%hite(ic),cpatch%bstorage(ic)) &
                          * cpatch%nplant(ic)
                    endif
                 enddo
              
              end if
           enddo loop_patches
        enddo loop_sites

        close(12)

        !! Init sites, patches, and cohorts
        !! Check cohorts are not bare ground
        do isi = 1,cpoly%nsites
           
           area_sum = 0.0
           site_lai = 0.0
           ncohorts = 0


           ! Normalize the area
           csite => cpoly%site(isi)
           area_tot = sum(csite%area(1:csite%npatches))
           csite%area(:)=csite%area(:)/area_tot

           do ipa=1,csite%npatches
              area_sum = area_sum + csite%area(ipa)
              patch_lai = 0.0
              
              cpatch => csite%patch(ipa)
              do ico = 1,cpatch%ncohorts
                 patch_lai = patch_lai + cpatch%lai(ico)
                 ncohorts = ncohorts + 1
              enddo
              
              csite%lai(ipa) = patch_lai
              site_lai = site_lai + csite%area(ipa) * patch_lai
              
           enddo

           ! If there are no cohorts, set some up
           ! THERE ARE SOME DISTURBANCE TYPES,1,2 THAT HAVE NO COHORTS IN THEM
           ! SHOULD THESE BE POPULATED WITH STARTED COHORTS? THE LOGIC IN THE 
           ! NEXT FEW LINES CHECKS FOR COHORTS AT THE SITE LEVEL, NOT THE PATCH
           ! LEVEL.  IS THIS OK?  -RGK 4-3-08

           ! MLO 5-27-08. I don't think so, there are patches that are still with no cohorts...
           ! I am just switching to the patch level, so it forces these patches to have a minimum 
           ! number of cohorts. I moved this to the time it allocates the cohorts.
                     
           do ipa = 1,csite%npatches
              
              cpatch => csite%patch(ipa)
              do ico = 1,cpatch%ncohorts
                 call init_ed_cohort_vars_array(cpatch,ico,cpoly%lsl(isi))
              enddo
              
           enddo
           
           call init_ed_patch_vars_array(csite,1,csite%npatches)
           
        enddo
        
        call init_ed_site_vars_array(cpoly,cgrid%lat(ipy))

        !  Get a diagnostic on the polygon's vegetation
        
        poly_lai = 0.0
        ncohorts = 0
        
        do isi = 1,cpoly%nsites
           
           nsitepat = 0
           csite => cpoly%site(isi)
           
           do ipa = 1,csite%npatches
              
              csite%lai(ipa) = 0.0
              npatchco = 0
              cpatch => csite%patch(ipa)
              
              do ico = 1,cpatch%ncohorts
                 ncohorts=ncohorts+1
                 npatchco=npatchco+1
                 csite%lai(ipa) = csite%lai(ipa) + cpatch%lai(ico)

              enddo

              poly_lai = poly_lai + cpoly%area(isi)*csite%area(ipa)*csite%lai(ipa)
              csite%cohort_count(ipa) = npatchco
              nsitepat = nsitepat + 1
              
           enddo
           cpoly%patch_count(isi) = nsitepat
        enddo
        
     enddo


     
     !! need to check what's going on in here
     call init_ed_poly_vars_array(cgrid)

  enddo

  return
end subroutine read_ed1_history_file_array
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine init_full_history_restart()


  use max_dims, only: n_pft
  use pft_coms, only: SLA, q, qsw, hgt_min, include_pft, phenology
  use misc_coms, only: sfilin, ied_init_mode,current_time
  use mem_sites, only: grid_res,edres
  use consts_coms, only: pio180
  use ed_misc_coms, only: use_target_year, restart_target_year
  use ed_state_vars,only: polygontype,sitetype,patchtype,edtype, &
       edgrid_g,allocate_sitetype,allocate_patchtype,allocate_polygontype
  use soil_coms, only: alloc_soilgrid
  use grid_coms,only:ngrids
  use ed_node_coms, only: mynum,nmachs,nnodetot,mchnum,machs
  use hdf5
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       globdims,chnkdims,chnkoffs

  implicit none
  
  integer :: year
  real, external :: dbh2h
  character(len=1)  :: vnam
  character(len=2)  :: cgr
  character(len=128) :: hnamel
  type(edtype),pointer      :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer    :: csite
  type(patchtype),pointer   :: cpatch
  logical :: exists ! File existence
  real,allocatable :: file_lats(:),file_lons(:)
  integer,allocatable :: pysi_n(:),pysi_id(:)
  integer,allocatable :: sipa_n(:),sipa_id(:)
  integer,allocatable :: paco_n(:),paco_id(:)
  
  real :: ll_tolerance

  integer :: ngr,ifpy,ipft
  integer :: ipy,isi,ipa,ico
  integer :: py_index,si_index,pa_index

  ! HDF5 types are defined here
  integer :: hdferr
  include 'mpif.h'
  integer,       dimension(MPI_STATUS_SIZE) :: status
  integer                                   :: ierr
  integer :: igr
  real(kind=8) :: dbletime

 

  ! ------------------------------------------------------------------------------
  ! There are two types of history restarts that can be done.  An exact restart
  ! Assumes that you are continuing with the exact same configuration as a given
  ! model simulation that wrote the file in which you are using.  In this
  ! case, each node starts with a list of polygons, and searches the HDF history
  ! file for these polygons.  It expects to find at least one polygon in the file 
  ! within 250 meters of proximity.  If it does not find this it stops. It then fills
  ! each of these polygons with data, and traverses the data hierarchical tree
  ! that roots from each polygon, initializing the model to the exact same state
  ! as the end of the previous run.  A search based restart does not expect to 
  ! find exact matches, and may not use all of the polygons in the file. But none
  ! the less, it will traverse the tree from these polygons and populate the model 
  ! states with what is found in the files tree.
  ! -------------------------------------------------------------------------------
  ! Currently, we are not doing collective reads from the dataset, this may not
  ! be a feasible option, because at least during the polygon read-in, we are not
  ! sure which chunks to read in, so we take the whole vector of polygons for each
  ! node. Note that all nodes read the history restart separately.
  ! -------------------------------------------------------------------------------

  
  
  ! Set the tolerance on a matched latitude or longitude (100 meters)

  ll_tolerance = (1.0/115.0)*(1.0/10.0)   
  ! at equator: (1 degree / 115 kilometers)  (1 km / 10 100-meter invervals)

  ! Open the HDF environment

  call h5open_f(hdferr)


  ! Construct the file name for reinitiatlizing from
  ! The history file
  
  vnam = 'S'
  
  do ngr=1,ngrids
     
     cgrid => edgrid_g(1)

     print*,"================================================"
     print*,"      Entering Full History Initialization      "
     

     !=======================================
     ! 1) Open the HDF5 HISTORY FILE
     !=======================================

     write(cgr,'(a1,i1)') 'g',ngr

     dbletime=current_time%time
     
     call makefnam(hnamel,sfilin,dbletime,current_time%year, &
          current_time%month,current_time%date,0,vnam,cgr,'h5 ')

     inquire(file=trim(hnamel),exist=exists)

     if (.not.exists) then
        call fatal_error ('File '//trim(hnamel)//' not found.'         &
                         ,'init_full_history_restart','ed_history_io.f90')
     else
        call h5fopen_f(hnamel, H5F_ACC_RDONLY_F, file_id, hdferr)
        if (hdferr < 0) then
           print *, 'Error opening HDF5 file - error - ',hdferr
           print *, '   Filename: ',trim(hnamel)
           call fatal_error('Error opening HDF5 file - error - '//trim(hnamel) &
                           ,'init_full_history_restart','ed_history_io.f90')
        end if
     end if


     !=======================================
     ! 2) Retrieve global vector sizes
     !=======================================

     !
     ! TO DO!!!!  INCLUDE NZG,NPFT,NBBH AS PART OF
     !            THE HDF5 HEADER, AND COMPARE
     !            AS SANITY CHECK TO THOSE VALUES
     !            READ IN THE NAMELIST. IF THEY
     !            DONT MATCH THEN THE NAMELIST
     !            ENTRY SHOULD BE CHANGED....
     !


     
     globdims = 0
     chnkdims = 0
     chnkoffs = 0

     globdims(1) = 1

     call h5dopen_f(file_id,'NPOLYGONS_GLOBAL', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,cgrid%npolygons_global,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     call h5dopen_f(file_id,'NSITES_GLOBAL', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,cgrid%nsites_global,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     call h5dopen_f(file_id,'NPATCHES_GLOBAL', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,cgrid%npatches_global,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)

     call h5dopen_f(file_id,'NCOHORTS_GLOBAL', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,cgrid%ncohorts_global,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)

     !=======================================
     ! 3) Retrieve the mapping of the data tree
     !=======================================
     
     globdims = 0
     globdims(1) = cgrid%npolygons_global
     
     allocate(pysi_n(cgrid%npolygons_global))
     allocate(pysi_id(cgrid%npolygons_global))

     call h5dopen_f(file_id,'PYSI_N', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,pysi_n,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     call h5dopen_f(file_id,'PYSI_ID', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,pysi_id,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     globdims(1) = cgrid%nsites_global
     
     allocate(sipa_n(cgrid%nsites_global))
     allocate(sipa_id(cgrid%nsites_global))
     
     call h5dopen_f(file_id,'SIPA_N', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,sipa_n,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     call h5dopen_f(file_id,'SIPA_ID', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,sipa_id,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)

     globdims(1) = cgrid%npatches_global
     
     allocate(paco_n(cgrid%npatches_global))
     allocate(paco_id(cgrid%npatches_global))
     
     call h5dopen_f(file_id,'PACO_N', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,paco_n,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     call h5dopen_f(file_id,'PACO_ID', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,paco_id,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)


     ! ======================================
     ! 4) Retrieve the polygon coordinates data

     globdims(1) = cgrid%npolygons_global
     allocate(file_lats(cgrid%npolygons_global))
     allocate(file_lons(cgrid%npolygons_global))
     
     call h5dopen_f(file_id,'LATITUDE', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_REAL,file_lats,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)

     call h5dopen_f(file_id,'LONGITUDE', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_REAL,file_lons,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)

     ! ======================================
     ! 5) Loop the polygons in the model state
     !    and match them with those int he file
     !    After the match. Walk through the 
     !    data from that polygon and initialize.
     !    A polygon match must have both latitudes
     !    and longitudes within 100 meters
     
     do ipy = 1,cgrid%npolygons
        
        py_index = 0

        cpoly => cgrid%polygon(ipy)

        do ifpy = 1,cgrid%npolygons_global
           
           if ( abs(file_lats(ifpy)-cgrid%lat(ipy)) < ll_tolerance .and. &
                abs(file_lons(ifpy)-cgrid%lon(ipy)) < ll_tolerance ) py_index = ifpy

        enddo

        if (py_index==0) then
           print*,"COULD NOT MATCH A POLYGON WITH THE DATASET"
           print*,"STOPPING"
           print*,"GRID LATS: ",cgrid%lat
           print*,"GRID LONS: ",cgrid%lon
           print*,"FILE LATS: ",file_lats
           print*,"FILE LONS: ",file_lons
           call fatal_error('Mismatch between polygon and dataset'         &
                           ,'init_full_history_restart','ed_history_io.f90')
        endif

        
        ! ========================================
        ! Get all necessary polygon variables
        ! associated with this index for the
        ! current polygon, scalar reads

        call fill_history_grid(cgrid,ipy,py_index)

        if (pysi_n(py_index) > 0) then
           
  !         print*,"Allocating: ",pysi_n(py_index)," sites from polygon",py_index
           
           call allocate_polygontype(cpoly,pysi_n(py_index))
           
           ! ========================================
           ! Get all necessary site variables
           ! associated with this index for the
           ! current polygon, vector reads
           
           call fill_history_polygon(cpoly,pysi_id(py_index),cgrid%nsites_global)
           
           do isi = 1,cpoly%nsites
              csite => cpoly%site(isi)
              
              ! Calculate the index of this site's data in the HDF
              si_index = pysi_id(py_index) + isi - 1
              
              if (sipa_n(si_index) > 0) then
                 
 !                print*,"Allocating: ",sipa_n(si_index)," patches from site ",si_index
                 
                 call allocate_sitetype(csite,sipa_n(si_index))

                 ! ========================================
                 ! Get all necessary patch variables
                 ! associated with this index for the
                 ! current site
                 
                 call fill_history_site(csite,sipa_id(si_index),cgrid%npatches_global)

                 do ipa = 1,csite%npatches
                    cpatch => csite%patch(ipa)
                    
                    pa_index = sipa_id(si_index) + ipa - 1

                    if (paco_n(pa_index) > 0) then

!                       print*,"Allocating: ",paco_n(pa_index)," cohorts from patch ",pa_index
                       
                       call allocate_patchtype(cpatch,paco_n(pa_index))
                       
                       ! ========================================
                       ! Get all necessary patch variables
                       ! associated with this index for the
                       ! current site
                       
                       call fill_history_patch(cpatch,paco_id(pa_index),cgrid%ncohorts_global)

                       
                       do ipft = 1,n_pft
                          csite%old_stoma_data_max(ipft,ipa)%recalc = 1
                       enddo
                       

                    else

                       cpatch%ncohorts = 0
                       
                    endif
                 
                 enddo

              else

                 print*,"ATTEMPTING TO FILL SITE WITH PATCH VECTOR DATA"
                 print*,"NO PATCHES WERE FOUND in SIPA_N(SI_INDEX)"
                 print*,"THIS IS EXTREMELY UNLIKELY AND DOWNRIGHT WRONG,BOTH.."
                 stop
                 
              endif

           enddo
           
        else

           print*,"ATTEMPTING TO FILL A POLYGON WITH SITE VECTOR DATA"
           print*,"NO SITES WERE FOUND AT PYSI_N(PY_INDEX)"
           print*,"THIS IS EVEN MORE WRONG THAN HAVING NO PATCHES"
           print*,"AND THAT WAS REALLY REALLY WRONG"
           print*,"THIS IS WORSE THAN LETTING GIL BUY LUNCH FOR YOU"
           stop

        endif

        
     enddo


     call h5fclose_f(file_id, hdferr)
     if (hdferr.ne.0) then
         print*,"COULD NOT CLOSE THE HDF FILE"
         print*,hdferr
         stop	
     endif

     deallocate(file_lats,file_lons)
     deallocate(paco_n,paco_id)
     deallocate(sipa_n,sipa_id)
     deallocate(pysi_n,pysi_id )

  enddo

  
  ! Update some of the derived quantities (this may be redundant)
  ! Removing from here because it needs meteorological variables that aren't defined yet
  ! do ipy = 1,cgrid%npolygons
  !   
  !   cpoly => cgrid%polygon(ipy)
  !   
  !   do isi = 1,cpoly%nsites
  !      csite => cpoly%site(isi)
  !      do ipa = 1,csite%npatches
  !         call update_patch_derived_props_ar(csite, cpoly%lsl(isi), cpoly%met(isi)%rhos, ipa)
  !      enddo
  !      call update_site_derived_props_ar(cpoly, 0, isi)
  !   enddo
  !   call update_polygon_derived_props_ar(cgrid)
  ! enddo

  ! Close the HDF environment
  
  call h5close_f(hdferr)

  ! Initialize the disturbance transition rates
  
  write(*,'(a,i2.2)')'    Initializing anthropogenic disturbance forcing. Node: ',mynum
  call landuse_init_array


  return
end subroutine init_full_history_restart
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_history_grid(cgrid,ipy,py_index)

  use ed_state_vars,only: edtype,polygontype
  use grid_coms,only : nzg
  use max_dims,only : n_pft,n_dbh, n_dist_types
  use hdf5
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       globdims,chnkdims,chnkoffs,cnt,stride, &
       memdims,memoffs,memsize

  implicit none

  type(edtype),target ::       cgrid

  integer,intent(in) :: ipy,py_index
  integer :: iparallel,dsetrank

  iparallel = 0
  
  
  globdims = 0
  chnkdims = 0
  chnkoffs = 0
  memoffs  = 0
  memdims  = 0
  memsize  = 1

  dsetrank = 1

  ! These are the dimensions in the filespace
  ! itself. Global is the size of the dataset,
  ! chnkoffs is the offset of the chunk we
  ! are going to read.  Chnkdims is the size
  ! of the slab that is to be read.
  
  globdims(1) = cgrid%npolygons_global
  chnkdims(1) = 1
  chnkoffs(1) = py_index - 1

  ! These are the dimensions for the memory space
  ! this should essentially be the same dimensioning
  ! as the buffer that we are filling. This routine
  ! is just filling a scalar point in a vector
  ! of polygons.

  memdims(1)  = 1
  memoffs(1)  = 0
  memsize(1)  = 1


  call hdf_getslab_i(cgrid%lsl(ipy),'LSL ',dsetrank,iparallel)

  call hdf_getslab_r(cgrid%wbar(ipy),'WBAR ',dsetrank,iparallel)
  
  call hdf_getslab_r(cgrid%Te(ipy),'TE ',dsetrank,iparallel)

  call hdf_getslab_r(cgrid%zbar(ipy),'ZBAR ',dsetrank,iparallel)

  call hdf_getslab_r(cgrid%tau(ipy),'TAU ',dsetrank,iparallel)

  call hdf_getslab_r(cgrid%sheat(ipy),'SHEAT ',dsetrank,iparallel)

  call hdf_getslab_r(cgrid%baseflow(ipy),'BASEFLOW ',dsetrank,iparallel)

  call hdf_getslab_i(cgrid%load_adjacency(ipy),'LOAD_ADJACENCY ',dsetrank,iparallel)

  call hdf_getslab_r(cgrid%swliq(ipy),'SWLIQ ',dsetrank,iparallel)

  ! All daily and monthly variables need to be retrieved if you are loading there...
  
  if(associated(cgrid%dmean_gpp         )) &
       call hdf_getslab_r(cgrid%dmean_gpp(ipy)          ,'DMEAN_GPP '         ,dsetrank,iparallel)
  if(associated(cgrid%dmean_evap        )) &
       call hdf_getslab_r(cgrid%dmean_evap(ipy)         ,'DMEAN_EVAP '        ,dsetrank,iparallel)
  if(associated(cgrid%dmean_transp      )) &
       call hdf_getslab_r(cgrid%dmean_transp(ipy)       ,'DMEAN_TRANSP '      ,dsetrank,iparallel)
  if(associated(cgrid%dmean_sensible_vc )) &
       call hdf_getslab_r(cgrid%dmean_sensible_vc(ipy)  ,'DMEAN_SENSIBLE_VC ' ,dsetrank,iparallel)
  if(associated(cgrid%dmean_sensible_gc )) &
       call hdf_getslab_r(cgrid%dmean_sensible_gc(ipy)  ,'DMEAN_SENSIBLE_GC ' ,dsetrank,iparallel)
  if(associated(cgrid%dmean_sensible_ac )) &
       call hdf_getslab_r(cgrid%dmean_sensible_ac(ipy)  ,'DMEAN_SENSIBLE_AC ' ,dsetrank,iparallel)
  if(associated(cgrid%dmean_sensible    )) &
       call hdf_getslab_r(cgrid%dmean_sensible(ipy)     ,'DMEAN_SENSIBLE '    ,dsetrank,iparallel)
  if(associated(cgrid%dmean_plresp      )) &
       call hdf_getslab_r(cgrid%dmean_plresp(ipy)       ,'DMEAN_PLRESP '      ,dsetrank,iparallel)
  if(associated(cgrid%dmean_rh          )) &
       call hdf_getslab_r(cgrid%dmean_rh(ipy)           ,'DMEAN_RH '          ,dsetrank,iparallel)
  if(associated(cgrid%dmean_leaf_resp   )) &
       call hdf_getslab_r(cgrid%dmean_leaf_resp(ipy)    ,'DMEAN_LEAF_RESP '   ,dsetrank,iparallel)
  if(associated(cgrid%dmean_root_resp   )) &
       call hdf_getslab_r(cgrid%dmean_root_resp(ipy)    ,'DMEAN_ROOT_RESP '   ,dsetrank,iparallel)
  if(associated(cgrid%dmean_growth_resp )) &
       call hdf_getslab_r(cgrid%dmean_growth_resp(ipy)  ,'DMEAN_GROWTH_RESP ' ,dsetrank,iparallel)
  if(associated(cgrid%dmean_storage_resp)) &
       call hdf_getslab_r(cgrid%dmean_storage_resp(ipy) ,'DMEAN_STORAGE_RESP ',dsetrank,iparallel)
  if(associated(cgrid%dmean_vleaf_resp  )) &
       call hdf_getslab_r(cgrid%dmean_vleaf_resp(ipy)   ,'DMEAN_VLEAF_RESP '  ,dsetrank,iparallel)
  if(associated(cgrid%dmean_nep         )) &
       call hdf_getslab_r(cgrid%dmean_nep(ipy)          ,'DMEAN_NEP '         ,dsetrank,iparallel)
  if(associated(cgrid%dmean_fsw         )) &
       call hdf_getslab_r(cgrid%dmean_fsw(ipy)          ,'DMEAN_FSW '         ,dsetrank,iparallel)
  if(associated(cgrid%dmean_fsn         )) &
       call hdf_getslab_r(cgrid%dmean_fsn(ipy)          ,'DMEAN_FSN '         ,dsetrank,iparallel)
  if(associated(cgrid%mmean_gpp         )) &
       call hdf_getslab_r(cgrid%mmean_gpp(ipy)          ,'MMEAN_GPP '         ,dsetrank,iparallel)
  if(associated(cgrid%mmean_evap        )) &
       call hdf_getslab_r(cgrid%mmean_evap(ipy)         ,'MMEAN_EVAP '        ,dsetrank,iparallel)
  if(associated(cgrid%mmean_transp      )) &
       call hdf_getslab_r(cgrid%mmean_transp(ipy)       ,'MMEAN_TRANSP '      ,dsetrank,iparallel)
  if(associated(cgrid%mmean_sensible    )) &
       call hdf_getslab_r(cgrid%mmean_sensible(ipy)     ,'MMEAN_SENSIBLE '    ,dsetrank,iparallel)
  if(associated(cgrid%mmean_sensible_ac )) &
       call hdf_getslab_r(cgrid%mmean_sensible_ac(ipy)  ,'MMEAN_SENSIBLE_AC ' ,dsetrank,iparallel)
  if(associated(cgrid%mmean_sensible_gc )) &
       call hdf_getslab_r(cgrid%mmean_sensible_gc(ipy)  ,'MMEAN_SENSIBLE_GC ' ,dsetrank,iparallel)
  if(associated(cgrid%mmean_sensible_vc )) &
       call hdf_getslab_r(cgrid%mmean_sensible_vc(ipy)  ,'MMEAN_SENSIBLE_VC ' ,dsetrank,iparallel)
  if(associated(cgrid%mmean_nep         )) &
       call hdf_getslab_r(cgrid%mmean_nep(ipy)          ,'MMEAN_NEP '         ,dsetrank,iparallel)
  if(associated(cgrid%mmean_plresp      )) &
       call hdf_getslab_r(cgrid%mmean_plresp(ipy)       ,'MMEAN_PLRESP '      ,dsetrank,iparallel)
  if(associated(cgrid%mmean_rh          )) &
       call hdf_getslab_r(cgrid%mmean_rh(ipy)           ,'MMEAN_RH '          ,dsetrank,iparallel)
  if(associated(cgrid%mmean_leaf_resp   )) &
       call hdf_getslab_r(cgrid%mmean_leaf_resp(ipy)    ,'MMEAN_LEAF_RESP '   ,dsetrank,iparallel)
  if(associated(cgrid%mmean_root_resp   )) &
       call hdf_getslab_r(cgrid%mmean_root_resp(ipy)    ,'MMEAN_ROOT_RESP '   ,dsetrank,iparallel)
  if(associated(cgrid%mmean_growth_resp )) &
       call hdf_getslab_r(cgrid%mmean_growth_resp(ipy)  ,'MMEAN_GROWTH_RESP ' ,dsetrank,iparallel)
  if(associated(cgrid%mmean_storage_resp)) &
       call hdf_getslab_r(cgrid%mmean_storage_resp(ipy) ,'MMEAN_STORAGE_RESP ',dsetrank,iparallel)
  if(associated(cgrid%mmean_vleaf_resp  )) &
       call hdf_getslab_r(cgrid%mmean_vleaf_resp(ipy)   ,'MMEAN_VLEAF_RESP '  ,dsetrank,iparallel)
  if(associated(cgrid%stdev_gpp         )) &
       call hdf_getslab_r(cgrid%stdev_gpp(ipy)          ,'STDEV_GPP '         ,dsetrank,iparallel)
  if(associated(cgrid%stdev_evap        )) &
       call hdf_getslab_r(cgrid%stdev_evap(ipy)         ,'STDEV_EVAP '        ,dsetrank,iparallel)
  if(associated(cgrid%stdev_transp      )) &
       call hdf_getslab_r(cgrid%stdev_transp(ipy)       ,'STDEV_TRANSP '      ,dsetrank,iparallel)
  if(associated(cgrid%stdev_sensible    )) &
       call hdf_getslab_r(cgrid%stdev_sensible(ipy)     ,'STDEV_SENSIBLE '    ,dsetrank,iparallel)
  if(associated(cgrid%stdev_nep         )) &
       call hdf_getslab_r(cgrid%stdev_nep(ipy)          ,'STDEV_NEP '         ,dsetrank,iparallel)
  if(associated(cgrid%stdev_rh          )) &
       call hdf_getslab_r(cgrid%stdev_rh(ipy)           ,'STDEV_RH '          ,dsetrank,iparallel)

  ! Variables with 2 dimensions (nzg,npolygons)
  dsetrank    = 2
  globdims(1) = nzg
  chnkdims(1) = nzg
  memdims(1)  = nzg
  memsize(1)  = nzg
  chnkoffs(1) = 0
  memoffs(1)  = 0

  globdims(2)  = cgrid%npolygons_global
  chnkdims(2)  = 1
  chnkoffs(2)  = py_index - 1
  memdims(2)   = 1
  memsize(2)   = 1
  memoffs(2)   = 0

  ! Ryan - This one was (1,1) before, I changed to (1,ipy), is that correct?
  call hdf_getslab_i(cgrid%ntext_soil(1,ipy)       ,'NTEXT_SOIL '       ,dsetrank,iparallel)
  if(associated(cgrid%dmean_soil_temp)) &
     call hdf_getslab_r(cgrid%dmean_soil_temp(1,ipy)  ,'DMEAN_SOIL_TEMP '  ,dsetrank,iparallel)
  if(associated(cgrid%dmean_soil_water)) &
     call hdf_getslab_r(cgrid%dmean_soil_water(1,ipy) ,'DMEAN_SOIL_WATER ' ,dsetrank,iparallel)
  if(associated(cgrid%mmean_soil_temp)) &
     call hdf_getslab_r(cgrid%mmean_soil_temp(1,ipy)  ,'MMEAN_SOIL_TEMP '  ,dsetrank,iparallel)
  if(associated(cgrid%dmean_soil_water)) &
     call hdf_getslab_r(cgrid%mmean_soil_water(1,ipy) ,'MMEAN_SOIL_WATER ' ,dsetrank,iparallel)


  ! Variables with 2 dimensions (n_pft,npolygons)
  dsetrank    = 2
  globdims(1) = n_pft
  chnkdims(1) = n_pft
  memdims(1)  = n_pft
  memsize(1)  = n_pft
  chnkoffs(1) = 0
  memoffs(1)  = 0

  globdims(2)  = cgrid%npolygons_global
  chnkdims(2)  = 1
  chnkoffs(2)  = py_index - 1
  memdims(2)   = 1
  memsize(2)   = 1
  memoffs(2)   = 0

  if(associated(cgrid%lai_pft)) call hdf_getslab_r(cgrid%lai_pft(1,ipy) ,'LAI_PFT '       ,dsetrank,iparallel)
  if(associated(cgrid%lai_pft)) call hdf_getslab_r(cgrid%lai_pft(1,ipy) ,'MMEAN_LAI_PFT ' ,dsetrank,iparallel)
  if(associated(cgrid%agb_pft)) call hdf_getslab_r(cgrid%agb_pft(1,ipy) ,'AGB_PFT '       ,dsetrank,iparallel)
  if(associated(cgrid%ba_pft)) call hdf_getslab_r(cgrid%ba_pft(1,ipy) ,'BA_PFT '        ,dsetrank,iparallel)

  ! Variables with 2 dimensions (n_pft,npolygons)
  dsetrank    = 2
  globdims(1) = n_dist_types
  chnkdims(1) = n_dist_types
  memdims(1)  = n_dist_types
  memsize(1)  = n_dist_types
  chnkoffs(1) = 0
  memoffs(1)  = 0

  globdims(2)  = cgrid%npolygons_global
  chnkdims(2)  = 1
  chnkoffs(2)  = py_index - 1
  memdims(2)   = 1
  memsize(2)   = 1
  memoffs(2)   = 0

  if(associated(cgrid%dmean_gpp_lu))  call hdf_getslab_r(cgrid%dmean_gpp_lu(1,ipy) ,'DMEAN_GPP_LU ' ,dsetrank,iparallel)
  if(associated(cgrid%dmean_rh_lu ))  call hdf_getslab_r(cgrid%dmean_rh_lu(1,ipy)  ,'DMEAN_RH_LU '  ,dsetrank,iparallel)
  if(associated(cgrid%dmean_nep_lu))  call hdf_getslab_r(cgrid%dmean_nep_lu(1,ipy) ,'DMEAN_NEP_LU ' ,dsetrank,iparallel)
  if(associated(cgrid%mmean_gpp_lu))  call hdf_getslab_r(cgrid%mmean_gpp_lu(1,ipy) ,'MMEAN_GPP_LU ' ,dsetrank,iparallel)
  if(associated(cgrid%mmean_rh_lu ))  call hdf_getslab_r(cgrid%mmean_rh_lu(1,ipy)  ,'MMEAN_RH_LU '  ,dsetrank,iparallel)
  if(associated(cgrid%mmean_nep_lu))  call hdf_getslab_r(cgrid%mmean_nep_lu(1,ipy) ,'MMEAN_NEP_LU ' ,dsetrank,iparallel)


  ! Variables with 2 dimensions (n_dbh,npolygons)
  dsetrank    = 2
  globdims(1) = n_dbh
  chnkdims(1) = n_dbh
  memdims(1)  = n_dbh
  memsize(1)  = n_dbh
  chnkoffs(1) = 0
  memoffs(1)  = 0

  globdims(2)  = cgrid%npolygons_global
  chnkdims(2)  = 1
  chnkoffs(2)  = py_index - 1
  memdims(2)   = 1
  memsize(2)   = 1
  memoffs(2)   = 0

  if(associated(cgrid%dmean_gpp_dbh)) call hdf_getslab_r(cgrid%dmean_gpp_dbh(1,ipy) ,'DMEAN_GPP_DBH ' ,dsetrank,iparallel)
  if(associated(cgrid%mmean_gpp_dbh)) call hdf_getslab_r(cgrid%mmean_gpp_dbh(1,ipy) ,'MMEAN_GPP_DBH ' ,dsetrank,iparallel)

  ! SITE_ADJACENCY


  return
end subroutine fill_history_grid
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_history_polygon(cpoly,pysi_index,nsites_global)
  
  use ed_state_vars,only: polygontype
  use hdf5
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       globdims,chnkdims,chnkoffs,cnt,stride, &
       memdims,memoffs,memsize

  use grid_coms,only : nzg
  use max_dims,only : n_pft,n_dbh,n_dist_types
  
  implicit none
  
  type(polygontype),target :: cpoly
  integer,intent(in) :: pysi_index
  integer,intent(in) :: nsites_global
  integer :: iparallel
  integer :: dsetrank

  iparallel = 0
  
  dsetrank = 1
  globdims = 0
  chnkdims = 0
  chnkoffs = 0
  memoffs  = 0
  memdims  = 0
  memsize  = 1
  
  globdims(1) = nsites_global
  chnkdims(1) = cpoly%nsites
  chnkoffs(1) = pysi_index - 1
  memdims(1)  = cpoly%nsites
  memsize(1)  = cpoly%nsites
  memoffs(1)  = 0
  
  call hdf_getslab_i(cpoly%patch_count(1),'PATCH_COUNT ',dsetrank,iparallel)  
  call hdf_getslab_i(cpoly%sitenum(1),'SITENUM ',dsetrank,iparallel)
  call hdf_getslab_i(cpoly%fia_forestry(1),'FIA_FORESTRY ',dsetrank,iparallel)
  call hdf_getslab_i(cpoly%agri_species(1),'AGRI_SPECIES ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%agri_stocking(1),'AGRI_STOCKING ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%lambda_primary(1),'LAMBDA_PRIMARY ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%lambda_secondary(1),'LAMBDA_SECONDARY ',dsetrank,iparallel)
  call hdf_getslab_i(cpoly%plantation_species(1),'PLANTATION_SPECIES ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%plantation_stocking(1),'PLANTATION_STOCKING ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%reference_agb(1),'REFERENCE_AGB ',dsetrank,iparallel)
 ! call hdf_getslab_r(cpoly%first_lutime(1),'FIRST_LUTIME ',dsetrank,iparallel)
 ! call hdf_getslab_r(cpoly%last_lutime(1),'LAST_LUTIME ',dsetrank,iparallel)
 ! call hdf_getslab_r(cpoly%clutime(1),' ',dsetrank,iparallel)
  call hdf_getslab_i(cpoly%lsl(1),'LSL_SI ',dsetrank,iparallel)   
  call hdf_getslab_r(cpoly%area(1),'AREA_SI ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%patch_area(1),'PATCH_AREA ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%elevation(1),'ELEVATION ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%slope(1),'SLOPE ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%aspect(1),'ASPECT ',dsetrank,iparallel)
    
  call hdf_getslab_i(cpoly%num_landuse_years(1),'NUM_LANDUSE_YEARS ',dsetrank,iparallel)
  call hdf_getslab_i(cpoly%soi(1),'SOI ',dsetrank,iparallel)
!  call hdf_getslab_r(cpoly%soi_name(1),' ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%TCI(1),'TCI ',dsetrank,iparallel)      
  call hdf_getslab_i(cpoly%hydro_next(1),'HYDRO_NEXT ',dsetrank,iparallel)
  call hdf_getslab_i(cpoly%hydro_prev(1),'HYDRO_PREV ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%moist_W(1),'MOIST_W ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%moist_f(1),'MOIST_F ',dsetrank,iparallel)  
  call hdf_getslab_r(cpoly%moist_tau(1),'MOIST_TAU ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%moist_zi(1),'MOIST_ZI ',dsetrank,iparallel) 
  call hdf_getslab_r(cpoly%baseflow(1),'BASEFLOW_SI ',dsetrank,iparallel) 
  call hdf_getslab_i(cpoly%metplex_beg_month(1),'METPLEX_BEG_MONTH ',dsetrank,iparallel)
  call hdf_getslab_i(cpoly%metplex_beg_year(1),'METPLEX_BEG_YEAR ',dsetrank,iparallel)
  call hdf_getslab_i(cpoly%metplex_end_year(1),'METPLEX_END_YEAR ',dsetrank,iparallel)

  call hdf_getslab_r(cpoly%min_monthly_temp(1),'MIN_MONTHLY_TEMP ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%removed_biomass(1),'REMOVED_BIOMASS ',dsetrank,iparallel) 
  call hdf_getslab_r(cpoly%harvested_biomass(1),'HARVESTED_BIOMASS ',dsetrank,iparallel) 
  call hdf_getslab_i(cpoly%plantation(1),'PLANTATION_SI ',dsetrank,iparallel) 
  call hdf_getslab_i(cpoly%agri_stocking_pft(1),'AGRI_STOCKING_PFT ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%agri_stocking_density(1),'AGRI_STOCKING_DENSITY ',dsetrank,iparallel)
  call hdf_getslab_i(cpoly%plantation_stocking_pft(1),'PLANTATION_STOCKING_PFT ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%plantation_stocking_density(1),'PLANTATION_STOCKING_DENSITY ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%primary_harvest_memory(1),'PRIMARY_HARVEST_MEMORY ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%secondary_harvest_memory(1),'SECONDARY_HARVEST_MEMORY ',dsetrank,iparallel)
  call hdf_getslab_i(cpoly%fire_flag(1),'FIRE_FLAG ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%fire_disturbance_rate(1),'FIRE_DISTURBANCE_RATE ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%fuel(1),'FUEL ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%ignition_rate(1),'IGNITION_RATE ',dsetrank,iparallel)
 ! call hdf_getslab_r(cpoly%phen_pars(1),'PHEN_PARS ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%treefall_disturbance_rate(1),'TREEFALL_DISTURBANCE_RATE ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%nat_disturbance_rate(1),'NAT_DISTURBANCE_RATE ',dsetrank,iparallel)
  call hdf_getslab_i(cpoly%nat_dist_type(1),'NAT_DIST_TYPE ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%disturbance_rate(1),'DISTURBANCE_RATE ',dsetrank,iparallel)

  dsetrank    = 2
  globdims(1) = n_pft
  chnkdims(1) = n_pft
  memdims(1)  = n_pft
  memsize(1)  = n_pft
  chnkoffs(1) = 0
  memoffs(1)  = 0
  globdims(2)  = nsites_global
  chnkdims(2)  = cpoly%nsites
  chnkoffs(2)  = pysi_index - 1
  memdims(2)   = cpoly%nsites
  memsize(2)   = cpoly%nsites
  memoffs(2)   = 0

  call hdf_getslab_r(cpoly%elongation_factor(1,1),'ELONGATION_FACTOR ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%delta_elongf(1,1),'DELTA_ELONGF ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%gee_phen_delay(1,1),'GEE_PHEN_DELAY ',dsetrank,iparallel)
  if (associated(cpoly%lai_pft)) call hdf_getslab_r(cpoly%lai_pft(1,1),'LAI_PFT_SI ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%green_leaf_factor(1,1),'GREEN_LEAF_FACTOR ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%leaf_aging_factor(1,1),'LEAF_AGING_FACTOR ',dsetrank,iparallel)

  dsetrank    = 2
  globdims(1) = nzg
  chnkdims(1) = nzg
  memdims(1)  = nzg
  memsize(1)  = nzg
  chnkoffs(1) = 0
  memoffs(1)  = 0
  globdims(2)  = nsites_global
  chnkdims(2)  = cpoly%nsites
  chnkoffs(2)  = pysi_index - 1
  memdims(2)   = cpoly%nsites
  memsize(2)   = cpoly%nsites
  memoffs(2)   = 0

  call hdf_getslab_i(cpoly%ntext_soil(1,1),'NTEXT_SOIL_SI ',dsetrank,iparallel)

  dsetrank    = 2
  globdims(1) = 12
  chnkdims(1) = 12
  memdims(1)  = 12
  memsize(1)  = 12
  chnkoffs(1) = 0
  memoffs(1)  = 0
  globdims(2)  = nsites_global
  chnkdims(2)  = cpoly%nsites
  chnkoffs(2)  = pysi_index - 1
  memdims(2)   = cpoly%nsites
  memsize(2)   = cpoly%nsites
  memoffs(2)   = 0

  call hdf_getslab_r(cpoly%lambda1(1,1),'LAMBDA1 ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%lambda_fire(1,1),'LAMBDA_FIRE ',dsetrank,iparallel)

  dsetrank    = 2
  globdims(1) = n_dist_types
  chnkdims(1) = n_dist_types
  memdims(1)  = n_dist_types
  memsize(1)  = n_dist_types
  chnkoffs(1) = 0
  memoffs(1)  = 0
  globdims(2)  = nsites_global
  chnkdims(2)  = cpoly%nsites
  chnkoffs(2)  = pysi_index - 1
  memdims(2)   = cpoly%nsites
  memsize(2)   = cpoly%nsites
  memoffs(2)   = 0

  call hdf_getslab_r(cpoly%lu_dist_area(1,1),'LU_DIST_AREA ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%loss_fraction(1,1),'LOSS_FRACTION ',dsetrank,iparallel)

  dsetrank    = 3
  globdims(1:2) = n_dist_types
  chnkdims(1:2) = n_dist_types
  memdims(1:2)  = n_dist_types
  memsize(1:2)  = n_dist_types
  chnkoffs(1:2) = 0
  memoffs(1:2)  = 0
  globdims(3)  = nsites_global
  chnkdims(3)  = cpoly%nsites
  chnkoffs(3)  = pysi_index - 1
  memdims(3)   = cpoly%nsites
  memsize(3)   = cpoly%nsites
  memoffs(3)   = 0
  
  call hdf_getslab_r(cpoly%disturbance_memory(1,1,1),'DISTURBANCE_MEMORY ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%disturbance_rates(1,1,1),'DISTURBANCE_RATES ',dsetrank,iparallel)
  
  dsetrank    = 3
  globdims(3) = nsites_global
  chnkdims(3) = cpoly%nsites
  chnkoffs(3) = pysi_index - 1
  memdims(3)  = cpoly%nsites
  memsize(3)  = cpoly%nsites
  memoffs(3)  = 0
  globdims(2) = n_dbh
  chnkdims(2) = n_dbh
  memdims(2)  = n_dbh
  memsize(2)  = n_dbh
  chnkoffs(2) = 0
  memoffs(2)  = 0
  globdims(1) = n_pft
  chnkdims(1) = n_pft
  memdims(1)  = n_pft
  memsize(1)  = n_pft
  chnkoffs(1) = 0
  memoffs(1)  = 0

  call hdf_getslab_r(cpoly%basal_area(1,1,1),'BASAL_AREA_SI ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%agb(1,1,1),'AGB_SI ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%basal_area_growth(1,1,1),'BASAL_AREA_GROWTH ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%agb_growth(1,1,1),'AGB_GROWTH ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%basal_area_mort(1,1,1),'BASAL_AREA_MORT ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%basal_area_cut(1,1,1),'BASAL_AREA_CUT ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%agb_mort(1,1,1),'AGB_MORT ',dsetrank,iparallel)
  call hdf_getslab_r(cpoly%agb_cut(1,1,1),'AGB_CUT ',dsetrank,iparallel)


  return
end subroutine fill_history_polygon
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_history_site(csite,sipa_index,npatches_global)

  use ed_state_vars,only: sitetype
  use grid_coms,only : nzg,nzs
  use max_dims,only : n_pft,n_dbh
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       globdims,chnkdims,chnkoffs,cnt,stride, &
       memdims,memoffs,memsize
  use fusion_fission_coms, only: ff_ndbh

  implicit none

  type(sitetype),target :: csite
  integer,intent(in) :: sipa_index
  integer,intent(in) :: npatches_global
  integer :: iparallel
  integer :: dsetrank
  
  iparallel = 0
  
  dsetrank = 1
  globdims = 0
  chnkdims = 0
  chnkoffs = 0
  memoffs  = 0
  memdims  = 0
  memsize  = 1

  ! These are the dimensions in the filespace
  ! itself. Global is the size of the dataset,
  ! chnkoffs is the offset of the chunk we
  ! are going to read.  Chnkdims is the size
  ! of the slab that is to be read.
  
  globdims(1) = npatches_global
  chnkdims(1) = csite%npatches
  chnkoffs(1) = sipa_index - 1

  memdims(1)  = csite%npatches
  memsize(1)  = csite%npatches
  memoffs(1)  = 0

  call hdf_getslab_i(csite%dist_type(1),'DIST_TYPE ',dsetrank,iparallel)
  call hdf_getslab_r(csite%age(1),'AGE ',dsetrank,iparallel)
  call hdf_getslab_r(csite%area(1),'AREA ',dsetrank,iparallel)
  call hdf_getslab_r(csite%fast_soil_C(1),'FAST_SOIL_C ',dsetrank,iparallel)
  call hdf_getslab_r(csite%slow_soil_C(1),'SLOW_SOIL_C ',dsetrank,iparallel)
  call hdf_getslab_r(csite%structural_soil_C(1),'STRUCTURAL_SOIL_C ',dsetrank,iparallel)
  call hdf_getslab_r(csite%structural_soil_L(1),'STRUCTURAL_SOIL_L ',dsetrank,iparallel)
  call hdf_getslab_r(csite%mineralized_soil_N(1),'MINERALIZED_SOIL_N ',dsetrank,iparallel)
  call hdf_getslab_r(csite%fast_soil_N(1),'FAST_SOIL_N ',dsetrank,iparallel)
  call hdf_getslab_r(csite%sum_dgd(1),'SUM_DGD ',dsetrank,iparallel)
  call hdf_getslab_r(csite%sum_chd(1),'SUM_CHD ',dsetrank,iparallel)
  call hdf_getslab_i(csite%plantation(1),'PLANTATION ',dsetrank,iparallel)
!  call hdf_getslab_i(csite%cohort_count(1),'COHORT_COUNT ',dsetrank,iparallel)
  call hdf_getslab_r(csite%can_temp(1),'CAN_TEMP ',dsetrank,iparallel)
  call hdf_getslab_r(csite%can_shv(1),'CAN_SHV ',dsetrank,iparallel)
  call hdf_getslab_r(csite%can_co2(1),'CAN_CO2 ',dsetrank,iparallel)
  call hdf_getslab_r(csite%can_depth(1),'CAN_DEPTH ',dsetrank,iparallel)
!  call hdf_getslab_i(csite%pname(1),'PNAME ',dsetrank,iparallel)
  call hdf_getslab_r(csite%lai(1),'LAI_PA ',dsetrank,iparallel)
  call hdf_getslab_i(csite%nlev_sfcwater(1),'NLEV_SFCWATER ',dsetrank,iparallel)
  call hdf_getslab_r(csite%ground_shv(1),'GROUND_SHV ',dsetrank,iparallel)
  call hdf_getslab_r(csite%surface_ssh(1),'SURFACE_SSH ',dsetrank,iparallel)
  call hdf_getslab_r(csite%rough(1),'ROUGH ',dsetrank,iparallel)
  call hdf_getslab_r(csite%avg_daily_temp(1),'AVG_DAILY_TEMP ',dsetrank,iparallel)  
  call hdf_getslab_r(csite%mean_rh(1),'MEAN_RH ',dsetrank,iparallel)
  call hdf_getslab_r(csite%mean_nep(1),'MEAN_NEP ',dsetrank,iparallel)
  call hdf_getslab_r(csite%wbudget_loss2atm(1),'WBUDGET_LOSS2ATM ',dsetrank,iparallel)
  call hdf_getslab_r(csite%wbudget_precipgain(1),'WBUDGET_PRECIPGAIN ',dsetrank,iparallel)
  call hdf_getslab_r(csite%wbudget_loss2runoff(1),'WBUDGET_LOSS2RUNOFF ',dsetrank,iparallel)
  call hdf_getslab_r(csite%wbudget_initialstorage(1),'WBUDGET_INITIALSTORAGE ',dsetrank,iparallel)
  call hdf_getslab_r(csite%ebudget_latent(1),'EBUDGET_LATENT ',dsetrank,iparallel)
  call hdf_getslab_r(csite%ebudget_loss2atm(1),'EBUDGET_LOSS2ATM ',dsetrank,iparallel)
  call hdf_getslab_r(csite%ebudget_loss2runoff(1),'EBUDGET_LOSS2RUNOFF ',dsetrank,iparallel)
  call hdf_getslab_r(csite%ebudget_netrad(1),'EBUDGET_NETRAD ',dsetrank,iparallel)
  call hdf_getslab_r(csite%ebudget_precipgain(1),'EBUDGET_PRECIPGAIN ',dsetrank,iparallel)
  call hdf_getslab_r(csite%ebudget_initialstorage(1),'EBUDGET_INITIALSTORAGE ',dsetrank,iparallel)
  call hdf_getslab_r(csite%co2budget_initialstorage(1),'CO2BUDGET_INITIALSTORAGE ',dsetrank,iparallel)
  call hdf_getslab_r(csite%co2budget_loss2atm(1),'CO2BUDGET_LOSS2ATM ',dsetrank,iparallel)
  call hdf_getslab_r(csite%co2budget_gpp(1),'CO2BUDGET_GPP ',dsetrank,iparallel)
  call hdf_getslab_r(csite%co2budget_plresp(1),'CO2BUDGET_PLRESP ',dsetrank,iparallel)
  call hdf_getslab_r(csite%co2budget_rh(1),'CO2BUDGET_RH ',dsetrank,iparallel)
  call hdf_getslab_r(csite%dmean_A_decomp(1),'DMEAN_A_DECOMP ',dsetrank,iparallel)
  call hdf_getslab_r(csite%dmean_Af_decomp(1),'DMEAN_AF_DECOMP ',dsetrank,iparallel)
  call hdf_getslab_r(csite%veg_rough(1),'VEG_ROUGH ',dsetrank,iparallel)
  call hdf_getslab_r(csite%veg_height (1),'VEG_HEIGHT ',dsetrank,iparallel)
  call hdf_getslab_r(csite%fsc_in(1),'FSC_IN ',dsetrank,iparallel)
  call hdf_getslab_r(csite%ssc_in(1),'SSC_IN ',dsetrank,iparallel)
  call hdf_getslab_r(csite%ssl_in(1),'SSL_IN ',dsetrank,iparallel)
  call hdf_getslab_r(csite%fsn_in(1),'FSN_IN ',dsetrank,iparallel)
  call hdf_getslab_r(csite%total_plant_nitrogen_uptake(1),'TOTAL_PLANT_NITROGEN_UPTAKE ',dsetrank,iparallel)

!  call hdf_getslab_r(csite%rshort_g(1),'RSHORT_G ',dsetrank,iparallel)
!  call hdf_getslab_r(csite%rshort_g_beam(1),'RSHORT_G_BEAM ',dsetrank,iparallel)
!  call hdf_getslab_r(csite%rshort_g_diffuse(1),'RSHORT_G_DIFFUSE ',dsetrank,iparallel)
!  call hdf_getslab_r(csite%rlong_g(1),'RLONG_G ',dsetrank,iparallel)
!  call hdf_getslab_r(csite%rlong_g_surf(1),'RLONG_G_SURF ',dsetrank,iparallel)
!  call hdf_getslab_r(csite%rlong_g_incid(1),'RLONG_G_INCID ',dsetrank,iparallel)
!  call hdf_getslab_r(csite%rlong_s(1),'RLONG_S ',dsetrank,iparallel)
!  call hdf_getslab_r(csite%rlong_s_surf(1),'RLONG_S_SURF ',dsetrank,iparallel)
!  call hdf_getslab_r(csite%rlong_s_incid(1),'RLONG_S_INCID ',dsetrank,iparallel)

!  call hdf_getslab_r(csite%albedt(1),'ALBEDT ',dsetrank,iparallel)
!  call hdf_getslab_r(csite%albedo_beam(1),'ALBEDO_BEAM ',dsetrank,iparallel)
!  call hdf_getslab_r(csite%albedo_diffuse(1),'ALBEDO_DIFFUSE ',dsetrank,iparallel)
!  call hdf_getslab_r(csite%rlongup(1),'RLONGUP ',dsetrank,iparallel)
!  call hdf_getslab_r(csite%rlong_albedo(1),'RLONGUP_ALBEDO ',dsetrank,iparallel)
  call hdf_getslab_r(csite%total_snow_depth(1),'TOTAL_SNOW_DEPTH ',dsetrank,iparallel)
  call hdf_getslab_r(csite%snowfac(1),'SNOWFAC ',dsetrank,iparallel)
  call hdf_getslab_r(csite%A_decomp(1),'A_DECOMP ',dsetrank,iparallel)
  call hdf_getslab_r(csite%f_decomp(1),'F_DECOMP ',dsetrank,iparallel)
  call hdf_getslab_r(csite%rh(1),'RH ',dsetrank,iparallel)
  call hdf_getslab_r(csite%cwd_rh(1),'CWD_RH ',dsetrank,iparallel)
  call hdf_getslab_i(csite%fuse_flag(1),'FUSE_FLAG ',dsetrank,iparallel)
  call hdf_getslab_r(csite%plant_ag_biomass(1),'PLANT_AG_BIOMASS ',dsetrank,iparallel)
  call hdf_getslab_r(csite%mean_wflux(1),'MEAN_WFLUX ',dsetrank,iparallel)
  call hdf_getslab_r(csite%mean_latflux(1),'MEAN_LATFLUX ',dsetrank,iparallel)
  call hdf_getslab_r(csite%mean_hflux(1),'MEAN_HFLUX ',dsetrank,iparallel)
  call hdf_getslab_r(csite%mean_runoff(1),'MEAN_RUNOFF ',dsetrank,iparallel)
  call hdf_getslab_r(csite%mean_qrunoff(1),'MEAN_QRUNOFF ',dsetrank,iparallel)
  call hdf_getslab_r(csite%htry(1),'HTRY ',dsetrank,iparallel)

  dsetrank    = 2
  globdims(1) = nzs
  chnkdims(1) = nzs
  memdims(1)  = nzs
  memsize(1)  = nzs  
  chnkoffs(1) = 0
  memoffs(1)  = 0
  globdims(2) = npatches_global
  chnkdims(2) = csite%npatches
  chnkoffs(2) = sipa_index - 1
  memdims(2)  = csite%npatches
  memsize(2)  = csite%npatches
  memoffs(2)  = 0

  call hdf_getslab_r(csite%sfcwater_mass(1,1),'SFCWATER_MASS ',dsetrank,iparallel)
  call hdf_getslab_r(csite%sfcwater_energy(1,1),'SFCWATER_ENERGY ',dsetrank,iparallel)
  call hdf_getslab_r(csite%sfcwater_depth(1,1),'SFCWATER_DEPTH ',dsetrank,iparallel)
  call hdf_getslab_r(csite%rshort_s(1,1),'RSHORT_S ',dsetrank,iparallel)
  call hdf_getslab_r(csite%rshort_s_beam(1,1),'RSHORT_S_BEAM ',dsetrank,iparallel)
  call hdf_getslab_r(csite%rshort_s_diffuse(1,1),'RSHORT_S_DIFFUSE ',dsetrank,iparallel)
  call hdf_getslab_r(csite%sfcwater_tempk(1,1),'SFCWATER_TEMPK ',dsetrank,iparallel)
  call hdf_getslab_r(csite%sfcwater_fracliq(1,1),'SFCWATER_FRACLIQ ',dsetrank,iparallel)

  dsetrank    = 2
  globdims(1) = nzg
  chnkdims(1) = nzg
  memdims(1)  = nzg
  memsize(1)  = nzg
  chnkoffs(1) = 0
  memoffs(1)  = 0
  globdims(2) = npatches_global
  chnkdims(2) = csite%npatches
  chnkoffs(2) = sipa_index - 1
  memdims(2)  = csite%npatches
  memsize(2)  = csite%npatches
  memoffs(2)  = 0

  call hdf_getslab_i(csite%ntext_soil(1,1),'NTEXT_SOIL_PA ',dsetrank,iparallel)
  call hdf_getslab_r(csite%soil_energy(1,1),'SOIL_ENERGY_PA ',dsetrank,iparallel)
  call hdf_getslab_r(csite%soil_water(1,1),'SOIL_WATER_PA ',dsetrank,iparallel)
  call hdf_getslab_r(csite%soil_tempk(1,1),'SOIL_TEMPK_PA ',dsetrank,iparallel)
  call hdf_getslab_r(csite%soil_fracliq(1,1),'SOIL_FRACLIQ_PA ',dsetrank,iparallel)

  dsetrank    = 2
  globdims(1) = n_pft
  chnkdims(1) = n_pft
  memdims(1)  = n_pft
  memsize(1)  = n_pft
  chnkoffs(1) = 0
  memoffs(1)  = 0
  globdims(2) = npatches_global
  chnkdims(2) = csite%npatches
  chnkoffs(2) = sipa_index - 1
  memdims(2)  = csite%npatches
  memsize(2)  = csite%npatches
  memoffs(2)  = 0
  
  call hdf_getslab_r(csite%A_o_max(1,1),'A_O_MAX ',dsetrank,iparallel) 
  call hdf_getslab_r(csite%A_c_max(1,1),'A_C_MAX ',dsetrank,iparallel) 
  call hdf_getslab_r(csite%repro(1,1),'REPRO_PA ',dsetrank,iparallel)


  dsetrank    = 2
  globdims(1) = n_dbh
  chnkdims(1) = n_dbh
  memdims(1)  = n_dbh
  memsize(1)  = n_dbh
  chnkoffs(1) = 0
  memoffs(1)  = 0
  globdims(2) = npatches_global
  chnkdims(2) = csite%npatches
  chnkoffs(2) = sipa_index - 1
  memdims(2)  = csite%npatches
  memsize(2)  = csite%npatches
  memoffs(2)  = 0
  call hdf_getslab_r(csite%co2budget_gpp_dbh(1,1),'CO2BUDGET_GPP_DBH ',dsetrank,iparallel)

!!!! MAY NEED TO ADD THIS ONE
!  call hdf_getslab_r(csite%old_stoma_data_max(1,1),'OLD ',dsetrank,iparallel)

  dsetrank    = 3
  globdims(3) = npatches_global
  chnkdims(3) = csite%npatches
  chnkoffs(3) = sipa_index - 1

  memdims(3)  = csite%npatches
  memsize(3)  = csite%npatches
  memoffs(3)  = 0
  
  globdims(2) = ff_ndbh
  chnkdims(2) = ff_ndbh
  memdims(2)  = ff_ndbh
  memsize(2)  = ff_ndbh
  chnkoffs(2) = 0
  memoffs(2)  = 0

  globdims(1) = n_pft
  chnkdims(1) = n_pft
  memdims(1)  = n_pft
  memsize(1)  = n_pft
  chnkoffs(1) = 0
  memoffs(1)  = 0

  call hdf_getslab_r(csite%pft_density_profile,'PFT_DENSITY_PROFILE ',dsetrank,iparallel)

  return
end subroutine fill_history_site
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_history_patch(cpatch,paco_index,ncohorts_global)
  
  use ed_state_vars,only: patchtype
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       filespace,memspace, &
       globdims,chnkdims,chnkoffs,cnt,stride, &
       memdims,memoffs,memsize

  implicit none
  
  type(patchtype),target :: cpatch
  integer,intent(in) :: paco_index
  integer,intent(in) :: ncohorts_global
  integer :: iparallel,dsetrank

  iparallel = 0
  
  dsetrank = 1
  globdims = 0
  chnkdims = 0
  chnkoffs = 0
  memoffs  = 0
  memdims  = 0
  memsize  = 1
  
  globdims(1) = ncohorts_global
  chnkdims(1) = cpatch%ncohorts
  chnkoffs(1) = paco_index - 1

  memdims(1)  = cpatch%ncohorts
  memsize(1)  = cpatch%ncohorts
  memoffs(1)  = 0

  call hdf_getslab_i(cpatch%pft(1),'PFT ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%nplant(1),'NPLANT ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%hite(1),'HITE ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%dbh(1),'DBH ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%bdead(1),'BDEAD ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%bleaf(1),'BLEAF ',dsetrank,iparallel)
  call hdf_getslab_i(cpatch%phenology_status(1),'PHENOLOGY_STATUS ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%balive(1),'BALIVE ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%lai(1),'LAI_CO ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%bstorage(1),'BSTORAGE ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%cbr_bar(1),'CBR_BAR ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%veg_temp(1),'VEG_TEMP ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%veg_water(1),'VEG_WATER ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%mean_gpp(1),'MEAN_GPP ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%mean_leaf_resp(1),'MEAN_LEAF_RESP ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%mean_root_resp(1),'MEAN_ROOT_RESP ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%dmean_leaf_resp(1),'DMEAN_LEAF_RESP_CO ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%dmean_root_resp(1),'DMEAN_ROOT_RESP_CO ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%dmean_gpp(1),'DMEAN_GPP_CO ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%dmean_gpp_pot(1),'DMEAN_GPP_POT ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%dmean_gpp_max(1),'DMEAN_GPP_MAX ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%growth_respiration(1),'GROWTH_RESPIRATION ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%storage_respiration(1),'STORAGE_RESPIRATION ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%vleaf_respiration(1),'VLEAF_RESPIRATION ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%fsn(1),'FSN ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%monthly_dndt(1),'MONTHLY_DNDT ',dsetrank,iparallel)
  
  call hdf_getslab_r(cpatch%Psi_open(1),'PSI_OPEN ',dsetrank,iparallel)
  call hdf_getslab_i(cpatch%krdepth(1),'KRDEPTH ',dsetrank,iparallel)
  call hdf_getslab_i(cpatch%first_census(1),'FIRST_CENSUS ',dsetrank,iparallel)
  call hdf_getslab_i(cpatch%new_recruit_flag(1),'NEW_RECRUIT_FLAG ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%par_v(1),'PAR_V ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%par_v_beam(1),'PAR_V_BEAM ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%par_v_diffuse(1),'PAR_V_DIFFUSE ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%rshort_v(1),'RSHORT_V ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%rshort_v_beam(1),'RSHORT_V_BEAM ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%rshort_v_diffuse(1),'RSHORT_V_DIFFUSE ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%rlong_v(1),'RLONG_V ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%rlong_v_surf(1),'RLONG_V_SURF ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%rlong_v_incid(1),'RLONG_V_INCID ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%rb(1),'RB ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%A_open(1),'A_OPEN ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%A_closed(1),'A_CLOSED ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%Psi_closed(1),'PSI_CLOSED ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%rsw_open(1),'RSW_OPEN ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%rsw_closed(1),'RSW_CLOSED ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%fsw(1),'FSW ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%fs_open(1),'FS_OPEN ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%stomatal_resistance(1),'STOMATAL_RESISTANCE ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%maintenance_costs(1),'MAINTENANCE_COSTS ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%bseeds(1),'BSEEDS ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%leaf_respiration(1),'LEAF_RESPIRATION ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%root_respiration(1),'ROOT_RESPIRATION ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%hcapveg(1),'HCAPVEG ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%gpp(1),'GPP ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%paw_avg10d(1),'PAW_AVG10D ',dsetrank,iparallel)

  dsetrank    = 2
  globdims(1) = 13
  chnkdims(1) = 13
  chnkoffs(1) = 0
  memdims(1)  = 13
  memsize(1)  = 13
  memoffs(2)  = 0

  globdims(2) = ncohorts_global
  chnkdims(2) = cpatch%ncohorts
  chnkoffs(2) = paco_index - 1

  memdims(2)  = cpatch%ncohorts
  memsize(2)  = cpatch%ncohorts
  memoffs(2)  = 0


  call hdf_getslab_r(cpatch%cb(1,1),'CB ',dsetrank,iparallel)
  call hdf_getslab_r(cpatch%cb_max(1,1),'CB_MAX ',dsetrank,iparallel)

!  call hdf_getslab_i(cpatch%old_stoma_data(1),' ',dsetrank,iparallel)

  return
end subroutine fill_history_patch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel)
  
  use hdf5
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       filespace,memspace, &
       globdims,chnkdims,chnkoffs,cnt,stride, &
       memdims,memoffs,memsize
  
  
  implicit none
  
  real,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff

  integer :: hdferr,dsetrank
  integer :: iparallel
  character(len=*),intent(in) :: varn
  

  call h5dopen_f(file_id,trim(varn), dset_id, hdferr)
  if (hdferr /= 0) then
     write(unit=*,fmt=*) 'File_ID = ',file_id
     write(unit=*,fmt=*) 'Dset_ID = ',dset_id
     call fatal_error('Could not get the dataset for '//trim(varn)//'!!!' &
                     ,'hdf_getslab_r','ed_history_io.f90')
  end if
  
  call h5dget_space_f(dset_id,filespace,hdferr)
  if (hdferr /= 0) then
     call fatal_error('Could not get the hyperslabs filespace for '//trim(varn)//'!!!' &
                     ,'hdf_getslab_r','ed_history_io.f90')
  end if

  call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,chnkoffs, &
       chnkdims,hdferr)
  if (hdferr /= 0) then
     call fatal_error('Could not assign the hyperslabs filespace for '//trim(varn)//'!!!' &
                     ,'hdf_getslab_r','ed_history_io.f90')
  end if

  call h5screate_simple_f(dsetrank,memsize,memspace,hdferr)
  if (hdferr /= 0) then
     write(unit=*,fmt=*) 'Chnkdims = ',chnkdims
     write(unit=*,fmt=*) 'Dsetrank = ',dsetrank
     call fatal_error('Could not create the hyperslabs memspace for '//trim(varn)//'!!!' &
                     ,'hdf_getslab_r','ed_history_io.f90')
  end if

  call h5sselect_hyperslab_f(memspace,H5S_SELECT_SET_F,memoffs, &
       memdims,hdferr)
  if (hdferr /= 0) then
     call fatal_error('Could not assign the hyperslabs filespace for '//trim(varn)//'!!!' &
                     ,'hdf_getslab_r','ed_history_io.f90')
  end if

  if (iparallel == 1) then
     
     call h5dread_f(dset_id, H5T_NATIVE_REAL,buff,globdims, hdferr, &
          mem_space_id = memspace, file_space_id = filespace, &
          xfer_prp = plist_id)
     if (hdferr /= 0) then
        select case (trim(varn))
        case('DMEAN_GPP','DMEAN_EVAP','DMEAN_TRANSP','DMEAN_SENSIBLE_VC','DMEAN_SENSIBLE_GC'    &
            ,'DMEAN_SENSIBLE_AC','DMEAN_SENSIBLE','DMEAN_PLRESP','DMEAN_RH','DMEAN_LEAF_RESP'   &
            ,'DMEAN_ROOT_RESP','DMEAN_GROWTH_RESP','DMEAN_STORAGE_RESP','DMEAN_VLEAF_RESP'      &
            ,'DMEAN_NEP','DMEAN_FSW','DMEAN_FSN','DMEAN_SOIL_TEMP','DMEAN_SOIL_WATER'           &
            ,'DMEAN_GPP_LU','DMEAN_RH_LU','DMEAN_NEP_LU','DMEAN_GPP_DBH')
           write (unit=*,fmt='(a)') '-----------------------------------------------------------'
           write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
           write (unit=*,fmt='(a)') '                                                           '
           write (unit=*,fmt='(a)') ' + Variable '//trim(varn)//' not found in your history.'
           write (unit=*,fmt='(a)') ' + Initializing this variable with zero. '
           write (unit=*,fmt='(a)') ' + This may cause your first daily and first monthly       '
           write (unit=*,fmt='(a)') '   output to be incorrect for this variable.'
           write (unit=*,fmt='(a)') '-----------------------------------------------------------'
           write (unit=*,fmt='(a)') '                                                           '
           buff=0.
        case('MMEAN_GPP','MMEAN_EVAP','MMEAN_TRANSP','MMEAN_SENSIBLE','MMEAN_NEP','MMEAN_PLRESP'&
            ,'MMEAN_RH','MMEAN_LEAF_RESP','MMEAN_ROOT_RESP','MMEAN_GROWTH_RESP'                 &
            ,'MMEAN_STORAGE_RESP','MMEAN_VLEAF_RESP','STDEV_GPP','STDEV_EVAP','STDEV_TRANSP'    &
            ,'STDEV_SENSIBLE','STDEV_NEP','STDEV_RH','MMEAN_LAI_PFT','MMEAN_GPP_LU'             &
            ,'MMEAN_SOIL_TEMP','MMEAN_SOIL_WATER'                                               &
            ,'MMEAN_RH_LU','MMEAN_NEP_LU','MMEAN_GPP_DBH')
           write (unit=*,fmt='(a)') '-----------------------------------------------------------'
           write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
           write (unit=*,fmt='(a)') '                                                           '
           write (unit=*,fmt='(a)') ' + Variable '//trim(varn)//' not found in your history.'
           write (unit=*,fmt='(a)') ' + Initializing this variable with zero. '
           write (unit=*,fmt='(a)') ' + This may cause your first monthly output to be incorrect'
           write (unit=*,fmt='(a)') '   for this variable.'
           write (unit=*,fmt='(a)') '-----------------------------------------------------------'
           write (unit=*,fmt='(a)') '                                                           '
           buff=0.
        case default
           call fatal_error('Could not read in the hyperslab dataset for '//trim(varn)//'!!!' &
                           ,'hdf_getslab_r','ed_history_io.f90')
        end select
     end if
     
  else
     call h5dread_f(dset_id, H5T_NATIVE_REAL,buff,globdims, hdferr, &
          mem_space_id = memspace, file_space_id = filespace )
     if (hdferr /= 0) then
        select case (trim(varn))
        case('DMEAN_GPP','DMEAN_EVAP','DMEAN_TRANSP','DMEAN_SENSIBLE_VC','DMEAN_SENSIBLE_GC'    &
            ,'DMEAN_SENSIBLE_AC','DMEAN_SENSIBLE','DMEAN_PLRESP','DMEAN_RH','DMEAN_LEAF_RESP'   &
            ,'DMEAN_ROOT_RESP','DMEAN_GROWTH_RESP','DMEAN_STORAGE_RESP','DMEAN_VLEAF_RESP'      &
            ,'DMEAN_NEP','DMEAN_FSW','DMEAN_FSN','DMEAN_SOIL_TEMP','DMEAN_SOIL_WATER'           &
            ,'DMEAN_GPP_LU','DMEAN_RH_LU','DMEAN_NEP_LU','DMEAN_GPP_DBH')
           write (unit=*,fmt='(a)') '-----------------------------------------------------------'
           write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
           write (unit=*,fmt='(a)') '                                                           '
           write (unit=*,fmt='(a)') ' + Variable '//trim(varn)//' not found in your history.'
           write (unit=*,fmt='(a)') ' + Initializing this variable with zero. '
           write (unit=*,fmt='(a)') ' + This may cause your first daily and first monthly       '
           write (unit=*,fmt='(a)') '   output to be incorrect for this variable.'
           write (unit=*,fmt='(a)') '-----------------------------------------------------------'
           write (unit=*,fmt='(a)') '                                                           '
           buff=0.
        case('MMEAN_GPP','MMEAN_EVAP','MMEAN_TRANSP','MMEAN_SENSIBLE','MMEAN_NEP','MMEAN_PLRESP'&
            ,'MMEAN_RH','MMEAN_LEAF_RESP','MMEAN_ROOT_RESP','MMEAN_GROWTH_RESP'                 &
            ,'MMEAN_STORAGE_RESP','MMEAN_VLEAF_RESP','STDEV_GPP','STDEV_EVAP','STDEV_TRANSP'    &
            ,'STDEV_SENSIBLE','STDEV_NEP','STDEV_RH','MMEAN_LAI_PFT','MMEAN_GPP_LU'             &
            ,'MMEAN_SOIL_TEMP','MMEAN_SOIL_WATER'                                               &
            ,'MMEAN_RH_LU','MMEAN_NEP_LU','MMEAN_GPP_DBH')
           write (unit=*,fmt='(a)') '-----------------------------------------------------------'
           write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
           write (unit=*,fmt='(a)') '                                                           '
           write (unit=*,fmt='(a)') ' + Variable '//trim(varn)//' not found in your history.'
           write (unit=*,fmt='(a)') ' + Initializing this variable with zero. '
           write (unit=*,fmt='(a)') ' + This may cause your first monthly output to be incorrect'
           write (unit=*,fmt='(a)') '   for this variable.'
           write (unit=*,fmt='(a)') '-----------------------------------------------------------'
           write (unit=*,fmt='(a)') '                                                           '
           buff=0.
        case default
           call fatal_error('Could not read in the hyperslab dataset for '//trim(varn)//'!!!' &
                           ,'hdf_getslab_r','ed_history_io.f90')
        end select
     end if
  endif

!  write(unit=*,fmt='(a)') 'History start: Loading '//trim(varn)//'...'

  call h5sclose_f(filespace, hdferr)
  call h5sclose_f(memspace , hdferr)
  call h5dclose_f(dset_id  , hdferr)
  

  return
end subroutine hdf_getslab_r
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel)
  
  use hdf5
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       filespace,memspace, &
       globdims,chnkdims,chnkoffs,cnt,stride, &
       memdims,memoffs,memsize
  
  implicit none
  
  integer,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff

  integer :: hdferr,dsetrank
  integer :: iparallel
  character(len=*),intent(in) :: varn

  call h5dopen_f(file_id,trim(varn), dset_id, hdferr)
  if (hdferr /= 0) then
     write(unit=*,fmt=*) 'File_ID = ',file_id
     write(unit=*,fmt=*) 'Dset_ID = ',dset_id
     call fatal_error('Could not get the dataset for '//trim(varn)//'!!!' &
                     ,'hdf_getslab_i','ed_history_io.f90')
  endif
  
  call h5dget_space_f(dset_id,filespace,hdferr)
  if (hdferr /= 0) then
     call fatal_error('Could not get the hyperslabs filespace for '//trim(varn)//'!!!' &
                     ,'hdf_getslab_i','ed_history_io.f90')
  endif

  call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,chnkoffs, &
       chnkdims,hdferr)
  if (hdferr /= 0) then
     call fatal_error('Could not assign the hyperslabs filespace for '//trim(varn)//'!!!' &
                     ,'hdf_getslab_i','ed_history_io.f90')
  endif

  call h5screate_simple_f(dsetrank,memsize,memspace,hdferr)
  if (hdferr /= 0) then
     write(unit=*,fmt=*) 'Chnkdims = ',chnkdims
     write(unit=*,fmt=*) 'Dsetrank = ',dsetrank
     call fatal_error('Could not create the hyperslabs memspace for '//trim(varn)//'!!!' &
                     ,'hdf_getslab_i','ed_history_io.f90')
  endif

  call h5sselect_hyperslab_f(memspace,H5S_SELECT_SET_F,memoffs, &
       memdims,hdferr)
  if (hdferr /= 0) then
     call fatal_error('Could not assign the hyperslabs filespace for '//trim(varn)//'!!!' &
                     ,'hdf_getslab_i','ed_history_io.f90')
  end if

  if (iparallel == 1) then
     
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,buff,globdims, hdferr, &
          mem_space_id = memspace, file_space_id = filespace, &
          xfer_prp = plist_id)
     if (hdferr /= 0) then
        call fatal_error('Could not read in the hyperslab dataset for '//trim(varn)//'!!!' &
                        ,'hdf_getslab_i','ed_history_io.f90')
     end if
     
  else
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,buff,globdims, hdferr, &
          mem_space_id = memspace, file_space_id = filespace )
     if (hdferr /= 0) then
        call fatal_error('Could not read in the hyperslab dataset for '//trim(varn)//'!!!' &
                        ,'hdf_getslab_r','ed_history_io.f90')
     end if
  end if

!  write(unit=*,fmt='(a)') 'History start: Loading '//trim(varn)//'...'

  call h5sclose_f(filespace, hdferr)
  call h5sclose_f(memspace, hdferr)
  call h5dclose_f(dset_id, hdferr)
  
  return
end subroutine hdf_getslab_i
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine create_ed1_fname(lat, ed_res, lon, sfilin, pss_name, css_name,   &
     site_name)
  implicit none

  real, intent(in) :: lat
  real, intent(in) :: lon
  real, intent(in) :: ed_res
  real :: flon
  real :: flat
  character(len=256) :: ed_fname
  character(len=*), intent(in) :: sfilin
  character(len=256), intent(out) :: pss_name
  character(len=256), intent(out) :: css_name
  character(len=256), intent(out) :: site_name
 
  ! Make file name
  if(lat >= 0.0)then
     flat = ed_res * int(lat / ed_res) + 0.5 * ed_res 
  else
     flat = - ed_res * int(-lat / ed_res) - 0.5 * ed_res
  endif
  
  if(lon >= 0.0)then
     flon = ed_res * int(lon / ed_res) + 0.5 * ed_res 
  else
     flon = - ed_res * int(-lon / ed_res) - 0.5 * ed_res
  endif

  if(ed_res > 0.999 .and. ed_res < 1.001)then

     if(flat <= -10.0)then
        write(ed_fname,'(a,f5.1)')trim(sfilin)//'lat',flat
     elseif(flat < 0.0 .or. flat >= 10.0)then
        write(ed_fname,'(a,f4.1)')trim(sfilin)//'lat',flat
     else
        write(ed_fname,'(a,f3.1)')trim(sfilin)//'lat',flat
     endif
     if(flon <= -100.0)then
        write(ed_fname,'(a,f6.1)')trim(ed_fname)//'lon',flon
     elseif(flon <= -10.0 .or. flon >= 100.0)then
        write(ed_fname,'(a,f5.1)')trim(ed_fname)//'lon',flon
     elseif(flon < 0.0)then
        write(ed_fname,'(a,f4.1)')trim(ed_fname)//'lon',flon
     elseif(flon < 10.0)then
        write(ed_fname,'(a,f3.1)')trim(ed_fname)//'lon',flon
     else
        write(ed_fname,'(a,f4.1)')trim(ed_fname)//'lon',flon
     endif
  else
     print*,'bad ed_res in create_ed1_fname'
     stop
  endif

  pss_name = trim(ed_fname)//'.pss'
  css_name = trim(ed_fname)//'.css'
  site_name = trim(ed_fname)//'.site'

  return
end subroutine create_ed1_fname
