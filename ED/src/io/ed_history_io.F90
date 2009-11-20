subroutine read_ed1_history_file


  use ed_max_dims, only: n_pft,huge_patch,huge_cohort,max_water,str_len,maxfiles,maxlist
  use pft_coms, only: SLA, q, qsw, hgt_min, include_pft, include_pft_ag, &
       phenology,pft_1st_check,include_these_pft
  use ed_misc_coms, only: sfilin, ied_init_mode
  use mem_sites, only: grid_res,edres
  use consts_coms, only: pio180,pio4
  use ed_misc_coms, only: use_target_year, restart_target_year
  use ed_state_vars,only:polygontype,sitetype,patchtype,edtype, &
       edgrid_g,allocate_sitetype,allocate_patchtype
  use grid_coms,only:ngrids
  use allometry, only: dbh2h,h2dbh,dbh2bd,dbh2bl, ed_biomass,area_indices
  use fuse_fiss_utils, only: sort_cohorts
  use disturb_coms , only : min_new_patch_area ! ! intent(in)
  implicit none

  integer :: year

  type(edtype),pointer      :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer    :: csite
  type(patchtype),pointer   :: cpatch
  
  integer :: pft
  logical :: renumber_pfts
  character(len=str_len) :: pss_name
  character(len=str_len) :: css_name
  integer :: igr,ipy,isi,ipa,ico
  integer :: ip,ip2,ic,ic2

  character(len=str_len) :: cdum
  real :: dummy
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
  real(kind=8), parameter :: min_area = 1.d-7     ! Doesn't need to be the machine epsilon,
                                                  !     just chose a small number.
  real(kind=8), parameter :: min_ok   = 1.d-20    ! Chose a number small enough, but with some
                                                  !     room for multiplying by a small area.
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
  real :: patch_lai,patch_wpa,patch_wai,poly_lai
  real :: site_lai
  integer :: ncohorts,npatchco
  integer :: npatches,nsitepat,npatch2
  integer, parameter :: harvard_override = 0
  logical :: site_match
  real   , external  :: dist_gc
  integer :: nw
  integer :: ied_init_mode_local

  !----- Variables for new method to find the closest file --------------------------------!
  integer                                               :: nf,nflist,nflsite,nflpss,nflcss
  integer                                               :: nclosest
  character(len=str_len), dimension(maxlist)            :: full_list
  character(len=str_len), dimension(maxfiles)           :: site_list,pss_list,css_list
  real                  , dimension(maxfiles)           :: slon_list,slat_list
  real                  , dimension(maxfiles)           :: plon_list,plat_list
  real                  , dimension(maxfiles)           :: clon_list,clat_list
  real                  , dimension(maxfiles)           :: file_pdist,file_cdist
  !----- External function. ---------------------------------------------------------------!
  real                  , external                      :: sngloff
  !----------------------------------------------------------------------------------------!
  
  
  ied_init_mode_local = ied_init_mode 

  !----- Retrieve all files with the specified prefix. ------------------------------------!
  call ed_filelist(full_list,sfilin,nflist)

  
  !----- Retrieve LON/LAT information for sites -------------------------------------------!
  if (ied_init_mode == 3) then
     renumber_pfts = .false.
     call ed1_fileinfo('.site',nflist,full_list,nflsite,site_list,slon_list,slat_list)
  else
     renumber_pfts = .true.
  end if
  
  !----- Retrieve LON/LAT information for patches and cohorts -----------------------------!
  call ed1_fileinfo('.pss',nflist,full_list,nflpss,pss_list,plon_list,plat_list)
  call ed1_fileinfo('.css',nflist,full_list,nflcss,css_list,clon_list,clat_list)
  
  ! Loop over all grids, polygons, and sites

  gridloop: do igr = 1,ngrids

     cgrid => edgrid_g(igr)

     polyloop: do ipy = 1,cgrid%npolygons
        
        cpoly => cgrid%polygon(ipy)

        ! Override of LSL, ntext_soil
        if(harvard_override == 1)then

           cpoly%lsl(1)   = 1
           cgrid%lsl(ipy) = 1
           cgrid%ntext_soil(1,ipy) = 2
           cgrid%ntext_soil(2,ipy) = 2
           cgrid%ntext_soil(3,ipy) = 3
           cgrid%ntext_soil(4,ipy) = 3
        end if

        ! Intialize the distances as very large
        file_pdist = 1e20
        file_cdist = 1e20
        

        ! =================================
        ! Part I: Find the restart files
        ! =================================
        do nf=1,nflpss
           file_pdist(nf)=dist_gc(cgrid%lon(ipy),plon_list(nf),cgrid%lat(ipy),plat_list(nf))
        end do

        do nf=1,nflcss
           file_cdist(nf)=dist_gc(cgrid%lon(ipy),clon_list(nf),cgrid%lat(ipy),clat_list(nf))
        end do

        ! ===================================
        ! Loop through the closest files
        ! until one is determined suitable. 
        ! Non suitable files are likely those 
        ! that are water patches.
        ! ===================================

        find_nonwater: do nf=1,nflpss

           ! nclosest is the file with the closest information. 
           nclosest=minloc(file_pdist,dim=1)
           pss_name = trim(pss_list(nclosest))
           
           write (unit=*,fmt='(2a)') 'Using patch file: ',trim(pss_name)
           
           ! =================================
           ! Part II: Add the patches
           ! =================================
           
           
           ! Open file and read in patches
           open(12,file=trim(pss_name),form='formatted',status='old')
           read(12,'(a4)')  cdum ! skip header
           
           !----- If running mixed ED-1/ED-2, check whether the file is ED1 or ED2. ----------!
           if (ied_init_mode == -1) then
              if (trim(cdum) == 'time') then
                 ied_init_mode_local = 2
              else
                 ied_init_mode_local = 1
              end if
           end if
           
           nwater = 1
           if(ied_init_mode_local == 1) then
              read(12,*)cdum,nwater
              read(12,*)cdum,depth(1:nwater)
              read(12,*)
           elseif(ied_init_mode_local == 2) then
           !   read(12,*)!water patch
           endif
           
           ! Note that if we are doing an ED1 restart we can 
           ! assume 1 site?
           
           ip = 1
           sitenum = 0
           count_patches: do
              
              if (ip>huge_patch) call fatal_error('IP too high,increase array size huge_patch', &
                   'read_ed1_history_file','ed_history_io.f90')
              
              select case (ied_init_mode_local)
              case (1) !! read ED1 format files
                 
                 read(12,*,iostat=ierr)  time(ip),pname(ip),trk(ip),dage,darea,dfsc,dstsc,  &
                      dstsl,dssc,dpsc,dmsn,dfsn,dwater(1:nwater)
                 if(ierr /= 0)exit count_patches
             
                 area(ip)   = sngloff(darea      ,min_area)
                 age(ip)    = sngloff(dage       ,min_ok  )
                 fsc(ip)    = sngloff(dfsc       ,min_ok  )
                 stsc(ip)   = sngloff(dstsc      ,min_ok  )
                 stsl(ip)   = sngloff(dstsl      ,min_ok  )
                 ssc(ip)    = sngloff(dssc       ,min_ok  )
                 psc(ip)    = sngloff(dpsc       ,min_ok  )
                 msn(ip)    = sngloff(dmsn       ,min_ok  )
                 fsn(ip)    = sngloff(dfsn       ,min_ok  )
                 do nw=1,nwater
                    water(nw,ip) = sngloff(dwater(nw) ,min_ok  )
                 end do
                 
              case(2)  !! read ED2 format files
                 read(12,*,iostat=ierr)time(ip),pname(ip),trk(ip),dage,darea,dwater(1),dfsc,dstsc  &
                      ,dstsl,dssc,dummy,dmsn,dfsn
                 if(ierr /= 0)exit count_patches
              
                 area(ip)    = sngloff(darea    ,min_area)
                 age(ip)     = sngloff(dage     ,min_ok  )
                 fsc(ip)     = sngloff(dfsc     ,min_ok  )
                 stsc(ip)    = sngloff(dstsc    ,min_ok  )
                 stsl(ip)    = sngloff(dstsl    ,min_ok  )
                 ssc(ip)     = sngloff(dssc     ,min_ok  )
                 msn(ip)     = sngloff(dmsn     ,min_ok  )
                 fsn(ip)     = sngloff(dfsn     ,min_ok  )
                 water(1,ip) = sngloff(dwater(1),min_ok  )
                 
              case(3)
                 
                 read(12,*,iostat=ierr) sitenum(ip),time(ip),pname(ip),trk(ip),age(ip), &
                      darea,water(1,ip),fsc(ip),stsc(ip),stsl(ip),ssc(ip),psc(ip),msn(ip),fsn(ip)
                 if(ierr /= 0)exit count_patches
              
                 area(ip)=sngloff(darea, min_area)
                 
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
                print*,isi,cpoly%nsites,ip
                    print*,"error reading from patch file",trim(pss_name)
                    print*,cpoly%sitenum
                    print*,"site number", sitenum,"not found"
                    stop
                 endif
              case default !Nearly bare ground
                 exit count_patches
              end select
              
              if (area(ip) > min_area) ip = ip + 1 ! This will remove patches with tiny area that often cause trouble
              
           enddo count_patches
           
           npatches = max(ip-1,0)
           
           ! If there are no patches, that is really strange, or, it is possible that
           ! the file it grabbed from was all water.. In that case, then we should use
           ! The next closest file.

           close(12)
           if (npatches>0) then

              ! We have found a suitable file, we will break from the 
              ! loop, and find the name of the corresponding cohort file

              nclosest=minloc( abs(file_pdist(nclosest)-file_cdist),dim=1)
              css_name = trim(css_list(nclosest))
              write (unit=*,fmt='(2a)') 'Using cohort file: ',trim(css_name)
              exit find_nonwater
           else
              file_pdist(nclosest) = 1e20  !The closest file was no good, so we make it far away for now
           endif
           
           

        enddo find_nonwater


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
              csite%plantation(ip)         = 0
              csite%cohort_count(ip)       = 0
           enddo
           
           
           ! Initialize the cohort counts per patch
           csite%cohort_count(:) = 0
           
        endif
        
        close(12)

        ! =================================
        ! Part III: Add the cohorts
        ! =================================
        
        open(12,file=trim(css_name),form='formatted',status='old')
        read(12,'(a4)')  cdum ! skip header

        !----- If running mixed ED-1/ED-2, check whether the file is ED1 or ED2. ----------!
        if (ied_init_mode == -1) then
           if (trim(cdum) == 'time') then
              ied_init_mode_local = 2
           else
              ied_init_mode_local = 1
           end if
        else
           ied_init_mode_local = ied_init_mode
        end if

        if (ied_init_mode_local == 1) then
           read(12,*) ! Skip second line.
        end if
        ic = 0

!        open(12,file=trim(css_name),form='formatted',status='old')
!        read(12,*)  ! skip header
!        if(ied_init_mode_local /= 3) then
!           read(12,*)  ! skip header
!        end if
!        
!        ic = 0
        
        read_cohorts: do

           ic = ic + 1
           add_this_cohort(ic) = .true. ! .false.
           
           if (ic>huge_cohort) call fatal_error('IC too high','read_ed1_history_file','ed_history_io.f90')

           select case (ied_init_mode_local)
           case (1)
              read(12,*,iostat=ierr)ctime(ic),cpname(ic),cname(ic),&
                   dbh(ic),hite(ic),ipft(ic),nplant(ic),  &
                   bdead(ic),balive(ic),avgRg(ic),leaves_on(ic),&
                   cb(1:12,ic),cb_max(1:12,ic)
              if(ierr /= 0) exit read_cohorts  
           case (2,3)
              read(12,*,iostat=ierr)ctime(ic),cpname(ic),cname(ic),&
                   dbh(ic),hite(ic),ipft(ic),nplant(ic),  &
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
           if(use_target_year == 1 .and. year .ne. restart_target_year)then
              add_this_cohort(ic) = .false. 
              continue
           endif

         !! check that nplant > 0 
          if(nplant(ic) < tiny(1.0)) then
            add_this_cohort(ic) = .false.
            continue
          endif

           ! Find site and patch and start counting how many to allocate
           
           put_cohort:do isi=1,cpoly%nsites
              csite => cpoly%site(isi)
              do ipa=1,csite%npatches
                 if (include_pft(ipft(ic)) == 0) then
                    select case (pft_1st_check)
                    case(0)
                       write (unit=*,fmt='(a,1x,i5,1x,a)') &
                            'I found a cohort with PFT=',ipft(ic),' and it is not in your include_these_pft...'
                       call fatal_error('Invalid PFT in history file','read_ed1_history_file','ed_history_io.f90')
                    case(1)
                       write (unit=*,fmt='(a,1x,i5,1x,a)') &
                            'I found a cohort with PFT=',ipft(ic),'... Including this PFT in your include_these_pft...'
                       include_pft(ipft(ic)) = 1
                       include_these_pft(sum(include_pft)) = ipft(ic)
                       call sort_up(include_these_pft,n_pft)
                       if (ipft(ic) == 1 .or. ipft(ic) == 5) include_pft_ag(ipft(ic)) = 1
                    case(2)
                       write (unit=*,fmt='(a,1x,i5,1x,a)') &
                            'I found a cohort with PFT=',ipft(ic),'... Ignoring it...'
                       add_this_cohort(ic) = .false.
                    end select
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
                       if(ied_init_mode_local == 1) then
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
                       
                       cpatch%broot(ic2)  = cpatch%balive(ic2) * q(ipft(ic))               &
                                          /(1.0 + q(ipft(ic)) + qsw(ipft(ic2)) * cpatch%hite(ic2))

                       cpatch%bsapwood(ic2) = cpatch%balive(ic2) * qsw(ipft(ic)) * cpatch%hite(ic2) &
                                            /(1.0 + q(ipft(ic)) + qsw(ipft(ic2)) * cpatch%hite(ic2))

                       !print*,cpatch%lai(ic2),cpatch%bleaf(ic2),cpatch%nplant(ic2),SLA(ipft(ic)),ipft(ic)
                       
                       cpatch%lai(ic2) = cpatch%bleaf(ic2) * cpatch%nplant(ic2) *   &
                            SLA(ipft(ic))
if(cpatch%lai(ic2) /= cpatch%lai(ic2)) then
print*,"invalid initial LAI" 
                       print*,cpatch%lai(ic2),cpatch%bleaf(ic2),cpatch%nplant(ic2),SLA(ipft(ic)),ipft(ic)
                       stop
endif
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
                                  
                       !----- Assign LAI, WPA, and WAI ------------------------------------!
                       call area_indices(cpatch%nplant(ic2),cpatch%bleaf(ic2)              &
                                        ,cpatch%bdead(ic2),cpatch%balive(ic2)              &
                                        ,cpatch%dbh(ic2), cpatch%hite(ic2)                 &
                                        ,cpatch%pft(ic2), SLA(cpatch%pft(ic2))             &
                                        ,cpatch%lai(ic2),cpatch%wpa(ic2), cpatch%wai(ic2))
                       
                       cpatch%cb(1:12,ic2) = cb(1:12,ic)
                       cpatch%cb_max(1:12,ic2) = cb_max(1:12,ic)
                       cpatch%cb(13,ic2) = 0.0
                       cpatch%cb_max(13,ic2) = 0.0
                       
                       cpatch%agb(ic2) = ed_biomass(cpatch%bdead(ic2),cpatch%balive(ic2)   &
                                                   ,cpatch%bleaf(ic2),cpatch%pft(ic2)      &
                                                   ,cpatch%hite(ic2),cpatch%bstorage(ic2)) 
                       cpatch%basarea(ic2)  = cpatch%nplant(ic2) * pio4                    &
                                            * cpatch%dbh(ic2) * cpatch%dbh(ic2)
                       cpatch%dagb_dt(ic2)  = 0.
                       cpatch%dba_dt(ic2)   = 0.
                       cpatch%ddbh_dt(ic2)  = 0.
                       
                       ! Initialize other cohort variables. Some of them won't be updated
                       ! unless the lai goes above lai_min
                       cpatch%fsw(ic2)   = 1.0
                       cpatch%gpp(ic2)   = 0.0
                       cpatch%par_v(ic2) = 0.0
                       
                       csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)           &
                                                   + cpatch%agb(ic2) * cpatch%nplant(ic2)
                    endif
                 enddo
              end if
           enddo loop_patches
        enddo loop_sites

        close(12)

        !! Init sites, patches, and cohorts
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
              csite%lai(ipa)  = 0.0
              csite%wpa(ipa)  = 0.0
              csite%wai(ipa)  = 0.0
              
              cpatch => csite%patch(ipa)
              do ico = 1,cpatch%ncohorts
                 csite%lai(ipa)  = csite%lai(ipa) + cpatch%lai(ico)
                 csite%wpa(ipa)  = csite%wpa(ipa) + cpatch%wpa(ico)
                 csite%wai(ipa)  = csite%wai(ipa) + cpatch%wai(ico)
                 ncohorts  = ncohorts + 1
              end do
              site_lai = site_lai + csite%area(ipa) * csite%lai(ipa)
              
           enddo

           do ipa = 1,csite%npatches
              
              cpatch => csite%patch(ipa)
              do ico = 1,cpatch%ncohorts                 
                 call init_ed_cohort_vars(cpatch,ico,cpoly%lsl(isi))
              enddo
              !----- Need to sort cohorts by size. ----------------------------------------!
              call sort_cohorts(cpatch)
              
           enddo
           
           call init_ed_patch_vars(csite,1,csite%npatches,cpoly%lsl(isi))
           
        enddo
        
        call init_ed_site_vars(cpoly,cgrid%lat(ipy))
        
        !  Get a diagnostic on the polygon's vegetation
        
        poly_lai = 0.0
        ncohorts = 0

        do isi = 1,cpoly%nsites

           nsitepat = 0
           csite => cpoly%site(isi)

           do ipa = 1,csite%npatches
              
              csite%lai(ipa)  = 0.0
              csite%wpa(ipa)  = 0.0
              csite%wai(ipa)  = 0.0
              npatchco = 0
              cpatch => csite%patch(ipa)
              
              do ico = 1,cpatch%ncohorts
                 ncohorts=ncohorts+1
                 npatchco=npatchco+1
                 csite%lai(ipa)  = csite%lai(ipa)  + cpatch%lai(ico)
                 csite%wpa(ipa)  = csite%wpa(ipa)  + cpatch%wpa(ico)
                 csite%wai(ipa)  = csite%wai(ipa)  + cpatch%wai(ico)

              end do

              poly_lai = poly_lai + cpoly%area(isi)*csite%area(ipa)*csite%lai(ipa)
              csite%cohort_count(ipa) = npatchco
              nsitepat = nsitepat + 1
              
           enddo
           cpoly%patch_count(isi) = nsitepat
        enddo
        
        
     end do polyloop

     !! need to check what's going on in here
     call init_ed_poly_vars(cgrid)

  end do gridloop

  return
end subroutine read_ed1_history_file

!==========================================================================================!
!==========================================================================================!


subroutine read_ed21_history_file

  use ed_max_dims, only: n_pft,str_len
  use pft_coms, only: SLA, q, qsw, hgt_min, include_pft, include_pft_ag,&
       phenology,pft_1st_check,include_these_pft
  use ed_misc_coms, only: sfilin,current_time, imonthh,iyearh,idateh,itimeh
  use ed_state_vars,only:polygontype,sitetype,patchtype,edtype, &
       edgrid_g,allocate_polygontype,allocate_sitetype,allocate_patchtype
  use grid_coms,only:ngrids,nzg
  use hdf5
  use consts_coms,only:pio4
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       globdims,chnkdims,chnkoffs,memdims,memoffs,memsize
  use allometry, only : area_indices, ed_biomass
  use fuse_fiss_utils, only : terminate_cohorts

  implicit none

  integer :: year

  type(edtype),pointer      :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer    :: csite
  type(patchtype),pointer   :: cpatch
  
  integer :: igr,ipy,isi,ipa,ico
  integer :: k
  integer :: dset_npolygons_global
  integer :: dset_nsites_global
  integer :: dset_npatches_global
  integer :: dset_ncohorts_global
  integer :: dset_nzg

  character(len=1)  :: vnam
  character(len=3)  :: cgr
  character(len=str_len) :: hnamel

  integer,allocatable :: pysi_n(:),pysi_id(:)
  integer,allocatable :: sipa_n(:),sipa_id(:)
  integer,allocatable :: paco_n(:),paco_id(:)
  real,allocatable :: file_lats(:),file_lons(:)

  real,parameter :: ll_tolerance = 25.0
  real    :: minrad, currad
  real    :: elim_nplant, elim_lai
  integer :: ngr,ifpy,ipft
  integer :: py_index,si_index,pa_index
  integer :: dsetrank,iparallel

  ! HDF5 types are defined here
  logical :: exists ! File existence
  integer :: hdferr
  include 'mpif.h'
  real(kind=8) :: dbletime

  ! Function that computes the distance between two points.
  real, external :: dist_gc

  ! Open the HDF environment

  call h5open_f(hdferr)

  ! Initialize the dimensional control variables for the H5 slabs
  globdims = 0_8
  chnkdims = 0_8
  chnkoffs = 0_8
  memoffs  = 0_8
  memdims  = 0_8
  memsize  = 1_8  

  ! The file name should be the exact file that
  
  vnam = 'S'

  ! ======================================
  ! Walk the tree and pull data from the dataset
  
  gridloop: do igr = 1,ngrids
     
     cgrid => edgrid_g(igr)

     !=======================================
     ! 1) Open the HDF5 HISTORY FILE
     !=======================================

     write(cgr,'(a1,i2.2)') 'g',igr

     dbletime=0.d0
     
     !!call makefnam(hnamel,sfilin,dbletime,iyearh,imonthh,idateh,itimeh*100,vnam,cgr,'h5 ')

     hnamel = trim(sfilin)//"-"//trim(cgr)//".h5"


     inquire(file=trim(hnamel),exist=exists)

     if (.not.exists) then
        write (unit=*,fmt='(a,1x,a)')    'SFILIN  = ',trim(sfilin)
        write (unit=*,fmt='(a,1x,i4.4)') 'IYEARH  = ',iyearh
        write (unit=*,fmt='(a,1x,i2.2)') 'IMONTHH = ',imonthh
        write (unit=*,fmt='(a,1x,i2.2)') 'IDATEH  = ',idateh
        write (unit=*,fmt='(a,1x,i4.4)') 'ITIMEH  = ',itimeh
        call fatal_error ('File '//trim(hnamel)//' not found.'         &
                         ,'read_ed21_history_fill','ed_history_io.f90')
     else
        call h5fopen_f(hnamel, H5F_ACC_RDONLY_F, file_id, hdferr)
        if (hdferr < 0) then
           print *, 'Error opening HDF5 file - error - ',hdferr
           print *, '   Filename: ',trim(hnamel)
           call fatal_error('Error opening HDF5 file - error - '//trim(hnamel) &
                           ,'read_ed21_history_fill','ed_history_io.f90')
        end if
     end if


     !=================================================
     ! 2) Retrieve global vector sizes and mapping tree
     !=================================================
     
     globdims = 0_8
     chnkdims = 0_8
     chnkoffs = 0_8

     globdims(1) = 1_8

     call h5dopen_f(file_id,'NZG', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_nzg,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     call h5dopen_f(file_id,'NPOLYGONS_GLOBAL', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_npolygons_global,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     call h5dopen_f(file_id,'NSITES_GLOBAL', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_nsites_global,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     call h5dopen_f(file_id,'NPATCHES_GLOBAL', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_npatches_global,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     call h5dopen_f(file_id,'NCOHORTS_GLOBAL', dset_id, hdferr)
     call h5dget_space_f(dset_id, dspace_id, hdferr)
     call h5dread_f(dset_id, H5T_NATIVE_INTEGER,dset_ncohorts_global,globdims, hdferr)
     call h5sclose_f(dspace_id, hdferr)
     call h5dclose_f(dset_id, hdferr)
     
     globdims(1) = int(dset_npolygons_global,8)
     
     allocate(pysi_n(dset_npolygons_global))
     allocate(pysi_id(dset_npolygons_global))
     
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
     
     globdims(1) = int(dset_nsites_global,8)
     
     allocate(sipa_n(dset_nsites_global))
     allocate(sipa_id(dset_nsites_global))
     
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
     
     globdims(1) = int(dset_npatches_global,8)
     allocate(paco_n(dset_npatches_global))
     allocate(paco_id(dset_npatches_global))
     
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
     
     globdims(1) = int(dset_npolygons_global,8)
     allocate(file_lats(dset_npolygons_global))
     allocate(file_lons(dset_npolygons_global))
     
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

     polyloop: do ipy = 1,cgrid%npolygons
        
        py_index = 0
        cpoly => cgrid%polygon(ipy)
        minrad = 1.e20
        
        do ifpy = 1,dset_npolygons_global
           
           currad = dist_gc(file_lons(ifpy),cgrid%lon(ipy),file_lats(ifpy),cgrid%lat(ipy))
           
           if ( abs(file_lats(ifpy)-cgrid%lat(ipy)) < ll_tolerance .and. &
                abs(file_lons(ifpy)-cgrid%lon(ipy)) < ll_tolerance .and. &
                (currad <  minrad) ) then
              py_index = ifpy
              minrad   = currad
           end if
           
        end do
        

        if ( py_index.eq.0 .or. pysi_n(py_index)<1) then
           print*,"COULD NOT MATCH A POLYGON WITH THE DATASET"
           print*,"STOPPING"
           print*,"THIS IS THE ",ipy,"th POLYGON"
           print*,"GRID LATS: ",cgrid%lat(ipy)
           print*,"GRID LONS: ",cgrid%lon(ipy)

           call fatal_error('Mismatch between polygon and dataset'         &
                   ,'read_ed21_history_file','ed_history_io.f90')
        else
           
           ! A suitably close polygon has been located in the datasets
           ! Use these values, and its children values in sites, patchs and cohorts
           ! =======================================================================
           iparallel = 0
           
           ! Load 1D dataset
           dsetrank = 1
           globdims(1) = int(dset_npolygons_global,8)
           chnkdims(1) = 1_8
           chnkoffs(1) = int(py_index - 1,8)
           memdims(1)  = 1_8
           memoffs(1)  = 0_8
           memsize(1)  = 1_8
           
           call hdf_getslab_i(cgrid%load_adjacency(ipy),'LOAD_ADJACENCY ',&
                dsetrank,iparallel,.true.)
           call hdf_getslab_r(cgrid%wbar(ipy),'WBAR ',dsetrank,iparallel,.true.)
           
           ! Load the workload (2D)
           dsetrank    = 2
           globdims(1) = int(13,8)
           chnkdims(1) = int(13,8)
           memdims(1)  = int(13,8)
           memsize(1)  = int(13,8)
           chnkoffs(1) = 0_8
           memoffs(1)  = 0_8
           globdims(2) = int(pysi_n(py_index),8)
           chnkdims(2) = int(pysi_n(py_index),8)
           memdims(2)  = int(pysi_n(py_index),8)
           memsize(2)  = int(pysi_n(py_index),8)
           chnkoffs(2) = 0_8
           memoffs(2)  = 0_8
           
           call hdf_getslab_r(cgrid%workload(:,ipy),'WORKLOAD ',dsetrank,iparallel,.false.)



           ! Load the site adjacency dataset (3D)
           dsetrank    = 3
           globdims(1) = int(pysi_n(py_index),8)
           chnkdims(1) = int(pysi_n(py_index),8)
           memdims(1)  = int(pysi_n(py_index),8)
           memsize(1)  = int(pysi_n(py_index),8)
           chnkoffs(1) = 0_8
           memoffs(1)  = 0_8
           globdims(2) = int(pysi_n(py_index),8)
           chnkdims(2) = int(pysi_n(py_index),8)
           memdims(2)  = int(pysi_n(py_index),8)
           memsize(2)  = int(pysi_n(py_index),8)
           chnkoffs(2) = 0_8
           memoffs(2)  = 0_8
           globdims(3)  = int(dset_npolygons_global,8)
           chnkdims(3)  = 1_8
           chnkoffs(3)  = int(py_index - 1,8)
           memdims(3)   = 1_8
           memsize(3)   = 1_8
           memoffs(3)   = 0_8
           
!!!           call hdf_getslab_i(cgrid%site_adjacency(:,:,ipy),'SITE_ADJACENCY ',&
!!!                dsetrank,iparallel,.true.)
           
           ! Allocate the vector of sites in the polygon
           
           call allocate_polygontype(cpoly,pysi_n(py_index))     
           ! But, these sites may have variable soil types, and we 
           ! want to preserve that
           
           dsetrank    = 2_8
           globdims(1) = int(dset_nzg,8)  ! How many layers in the dataset?
           chnkdims(1) = int(1,8)         ! We are only extracting one layer
           memdims(1)  = int(1,8)       ! We only need memory for one layer
           memsize(1)  = int(1,8)       ! On both sides
           chnkoffs(1) = int(dset_nzg - 1,8) ! Take the top layer, not the bottom
           memoffs(1)  = 0_8
           globdims(2)  = int(dset_nsites_global,8)
           chnkdims(2)  = int(cpoly%nsites,8)
           chnkoffs(2)  = int(pysi_id(py_index) - 1,8)
           memdims(2)   = int(cpoly%nsites,8)
           memsize(2)   = int(cpoly%nsites,8)
           memoffs(2)   = 0_8
           
           call hdf_getslab_i(cpoly%ntext_soil(nzg,:),'NTEXT_SOIL_SI '&
                ,dsetrank,iparallel,.true.)
           
           globdims = 0_8
           chnkdims = 0_8
           chnkoffs = 0_8
           memoffs  = 0_8
           memdims  = 0_8
           memsize  = 1_8
           
           dsetrank = 1_8
           globdims(1) = int(dset_nsites_global,8)
           chnkdims(1) = int(cpoly%nsites,8)
           chnkoffs(1) = int(pysi_id(py_index) - 1,8)
           memdims(1)  = int(cpoly%nsites,8)
           memsize(1)  = int(cpoly%nsites,8)
           memoffs(1)  = 0_8
           
           call hdf_getslab_r(cpoly%area,'AREA_SI ',dsetrank,iparallel,.true.)
           call hdf_getslab_r(cpoly%moist_f,'MOIST_F ',dsetrank,iparallel,.true.)  
           call hdf_getslab_r(cpoly%moist_W,'MOIST_W ',dsetrank,iparallel,.true.)
           call hdf_getslab_r(cpoly%elevation,'ELEVATION ',dsetrank,iparallel,.true.)
           call hdf_getslab_r(cpoly%slope,'SLOPE ',dsetrank,iparallel,.true.)
           call hdf_getslab_r(cpoly%aspect,'ASPECT ',dsetrank,iparallel,.true.)
           call hdf_getslab_r(cpoly%TCI,'TCI ',dsetrank,iparallel,.true.)  
           call hdf_getslab_i(cpoly%patch_count,'PATCH_COUNT ',dsetrank,iparallel,.true.)  
           
           siteloop: do isi = 1,cpoly%nsites
              csite => cpoly%site(isi)
              
              ! Calculate the index of this site's data in the HDF
              si_index = pysi_id(py_index) + isi - 1
              
              if (sipa_n(si_index) > 0) then
                 
                 ! The soil layer in this case is use defined,
                 ! so take this from the grid level
                 ! variable, and not from the dataset.
                 cpoly%lsl(isi)  = cgrid%lsl(ipy)  ! Initialize lowest soil layer
                 
                 
                 ! Now fill the soil column based on the top layer data
                 do k=1,nzg
                    cpoly%ntext_soil(k,isi) = cpoly%ntext_soil(nzg,isi)
                 enddo
                 
                 ! Fill 1D polygon (site unique) level variables
                 
                 call allocate_sitetype(csite,sipa_n(si_index))   
                 
                 iparallel = 0
                 
                 dsetrank = 1
                 globdims(1) = int(dset_npatches_global,8)
                 chnkdims(1) = int(csite%npatches,8)
                 chnkoffs(1) = int(sipa_id(si_index) - 1,8)
                 memdims(1)  = int(csite%npatches,8)
                 memsize(1)  = int(csite%npatches,8)
                 memoffs(1)  = 0

                 ! Assign patch soils based off of the site level soils data
                 do k=1,nzg
                    csite%ntext_soil(k,:) = cpoly%ntext_soil(k,isi)
                 enddo
                 
                 call hdf_getslab_i(csite%dist_type,'DIST_TYPE ',dsetrank,iparallel,.true.)
                 call hdf_getslab_r(csite%age,'AGE ',dsetrank,iparallel,.true.)
                 call hdf_getslab_r(csite%area,'AREA ',dsetrank,iparallel,.true.)
                 call hdf_getslab_r(csite%fast_soil_C,'FAST_SOIL_C ',dsetrank,iparallel,.true.)
                 call hdf_getslab_r(csite%slow_soil_C,'SLOW_SOIL_C ',dsetrank,iparallel,.true.)
                 call hdf_getslab_r(csite%structural_soil_C,'STRUCTURAL_SOIL_C ',&
                      dsetrank,iparallel,.true.)
                 call hdf_getslab_r(csite%structural_soil_L,'STRUCTURAL_SOIL_L ',&
                      dsetrank,iparallel,.true.)
                 call hdf_getslab_r(csite%mineralized_soil_N,'MINERALIZED_SOIL_N ', &
                      dsetrank,iparallel,.true.)
                 call hdf_getslab_r(csite%fast_soil_N,'FAST_SOIL_N ',dsetrank,iparallel,.true.)
                 call hdf_getslab_r(csite%sum_dgd,'SUM_DGD ',dsetrank,iparallel,.true.)
                 call hdf_getslab_r(csite%sum_chd,'SUM_CHD ',dsetrank,iparallel,.true.)
                 call hdf_getslab_i(csite%plantation,'PLANTATION ',dsetrank,iparallel,.true.)
                 
                 patchloop: do ipa = 1,csite%npatches
                    cpatch => csite%patch(ipa)
                    
                    csite%lai(ipa)  = 0.0
                    csite%wpa(ipa)  = 0.0
                    csite%wai(ipa)  = 0.0
                    csite%plant_ag_biomass(ipa)  = 0.0

                    pa_index = sipa_id(si_index) + ipa - 1
                    
                    call allocate_patchtype(cpatch,paco_n(pa_index))

                    if (cpatch%ncohorts > 0) then
                       
                       dsetrank = 1
                       globdims(1) = int(dset_ncohorts_global,8)
                       chnkdims(1) = int(cpatch%ncohorts,8)
                       chnkoffs(1) = int(paco_id(pa_index) - 1,8)
                       memdims(1)  = int(cpatch%ncohorts,8)
                       memsize(1)  = int(cpatch%ncohorts,8)
                       memoffs(1)  = 0_8
                       
                       call hdf_getslab_r(cpatch%dbh,'DBH ',dsetrank,iparallel,.true.)
                       call hdf_getslab_r(cpatch%hite,'HITE ',dsetrank,iparallel,.true.)
                       call hdf_getslab_i(cpatch%pft,'PFT ',dsetrank,iparallel,.true.)
                       call hdf_getslab_r(cpatch%nplant,'NPLANT ',dsetrank,iparallel,.true.)
                       call hdf_getslab_r(cpatch%bdead,'BDEAD ',dsetrank,iparallel,.true.)
                       call hdf_getslab_r(cpatch%balive,'BALIVE ',dsetrank,iparallel,.true.)
                       call hdf_getslab_i(cpatch%phenology_status,'PHENOLOGY_STATUS ',dsetrank,iparallel,.true.)
                       call hdf_getslab_r(cpatch%bleaf,'BLEAF ',dsetrank,iparallel,.true.)
                       call hdf_getslab_r(cpatch%bstorage,'BSTORAGE ',dsetrank,iparallel,.true.)
                                            
                       dsetrank    = 2
                       globdims(1) = 13_8
                       chnkdims(1) = 13_8
                       chnkoffs(1) = 0_8
                       memdims(1)  = 13_8
                       memsize(1)  = 13_8
                       memoffs(2)  = 0_8
                       globdims(2) = int(dset_ncohorts_global,8)
                       chnkdims(2) = int(cpatch%ncohorts,8)
                       chnkoffs(2) = int(paco_id(pa_index) - 1,8)
                       memdims(2)  = int(cpatch%ncohorts,8)
                       memsize(2)  = int(cpatch%ncohorts,8)
                       memoffs(2)  = 0_8
                       
                       call hdf_getslab_r(cpatch%cb,'CB ',dsetrank,iparallel,.true.)
                       call hdf_getslab_r(cpatch%cb_max,'CB_MAX ',dsetrank,iparallel,.true.)
                       
                       
                       cpatch%dagb_dt              = 0.
                       cpatch%dba_dt               = 0.
                       cpatch%ddbh_dt              = 0.
                       cpatch%fsw                  = 1.0
                       cpatch%gpp                  = 0.0
                       cpatch%par_v                = 0.0
                       
                       cohortloop: do ico=1,cpatch%ncohorts
                          !----------------------------------------------------------------!
                          !    We will now check the PFT of each cohort, so we determine   !
                          ! if this is a valid PFT.  If not, then we must decide what we   !
                          ! should do...                                                   !
                          !----------------------------------------------------------------!
                          if (include_pft(cpatch%pft(ico)) == 0) then
                             select case(pft_1st_check)
                             case (0)
                                write (unit=*,fmt='(a,1x,i5,1x,a)')                        &
                                      'I found a cohort with PFT=',cpatch%pft(ico)         &
                                     ,' and it is not in your include_these_pft...'
                                call fatal_error('Invalid PFT in history file'             &
                                                ,'read_ed21_history_file'                  &
                                                ,'ed_history_io.f90')
                             case (1)
                                write (unit=*,fmt='(a,1x,i5,1x,a)')                        &
                                     'I found a cohort with PFT=',cpatch%pft(ico)          &
                                    ,'... Including this PFT in your include_these_pft...'
                                include_pft(cpatch%pft(ico)) = 1
                                include_these_pft(sum(include_pft)) = cpatch%pft(ico)
                                call sort_up(include_these_pft,n_pft)
                                if (cpatch%pft(ico) == 1 .or. cpatch%pft(ico) == 5) then
                                   include_pft_ag(cpatch%pft(ico)) = 1
                                end if
                             case (2)
                                write (unit=*,fmt='(a,1x,i5,1x,a)')                        &
                                      'I found a cohort with PFT=',cpatch%pft(ico)         &
                                     ,'... Ignoring it...'
                                !----------------------------------------------------------!
                                !    The way we will ignore this cohort is by setting its  !
                                ! nplant to zero, and calling the "terminate_cohorts"      !
                                ! subroutine right after this.                             !
                                !----------------------------------------------------------!
                                cpatch%nplant(ico) = 0.
                             end select
                          end if
                          cpatch%agb(ico) = ed_biomass(cpatch%bdead(ico),cpatch%balive(ico)   &
                               ,cpatch%bleaf(ico),cpatch%pft(ico)      &
                               ,cpatch%hite(ico),cpatch%bstorage(ico))

                          cpatch%basarea(ico)  = cpatch%nplant(ico) * pio4                    &
                               * cpatch%dbh(ico) * cpatch%dbh(ico)
                          
                          cpatch%broot(ico) = q(cpatch%pft(ico)) * cpatch%balive(ico)         &
                                            / ( 1.0 + q(cpatch%pft(ico))                      &
                                              + qsw(cpatch%pft(ico)) * cpatch%hite(ico))
                          
                          cpatch%bsapwood(ico) = qsw(cpatch%pft(ico)) * cpatch%balive(ico)    &
                                               * cpatch%hite(ico)                             &
                                               / ( 1.0 + q(cpatch%pft(ico))                   &
                                                 + qsw(cpatch%pft(ico)) * cpatch%hite(ico))
                          
                          !----- Assign LAI, WPA, and WAI ------------------------------------!
                          call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico)              &
                               ,cpatch%bdead(ico),cpatch%balive(ico)              &
                               ,cpatch%dbh(ico), cpatch%hite(ico)                 &
                               ,cpatch%pft(ico), SLA(cpatch%pft(ico)), cpatch%lai(ico) &
                               ,cpatch%wpa(ico), cpatch%wai(ico))

                          csite%lai(ipa)  = csite%lai(ipa) + cpatch%lai(ico)
                          csite%wpa(ipa)  = csite%wpa(ipa) + cpatch%wpa(ico)
                          csite%wai(ipa)  = csite%wai(ipa) + cpatch%wai(ico)
                          csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)        &
                                                      + cpatch%agb(ico)*cpatch%nplant(ico)

                          call init_ed_cohort_vars(cpatch,ico,cpoly%lsl(isi))
                          
                       end do cohortloop
                       call terminate_cohorts(csite,ipa,elim_nplant,elim_lai)

                    end if
                    
                 enddo patchloop
              else
                 call fatal_error('No patches?','read_ed21_history_file','ed_history_io.f90')
              endif
              
              call init_ed_patch_vars(csite,1,csite%npatches,cpoly%lsl(isi))
              
           enddo siteloop
        end if
        
        call init_ed_site_vars(cpoly,cgrid%lat(ipy))
        
     end do polyloop


     !! need to check what's going on in here
     call init_ed_poly_vars(cgrid)
     
     call h5fclose_f(file_id, hdferr)
     if (hdferr.ne.0) then
        print*,"COULD NOT CLOSE THE HDF FILE"
        print*,hdferr
        call fatal_error('Could not close the HDF file','read_ed21_history_file'&
                        ,'ed_history_io.F90')
     end if
     
     deallocate(file_lats,file_lons)
     deallocate(paco_n,paco_id)
     deallocate(sipa_n,sipa_id)
     deallocate(pysi_n,pysi_id )
     
  end do gridloop
  
  ! Close the HDF environment
  
  call h5close_f(hdferr)
  
  return
end subroutine read_ed21_history_file
!==========================================================================================!
!==========================================================================================!




!==========================================================================================!
!==========================================================================================!
subroutine init_full_history_restart()


  use ed_max_dims, only: n_pft
  use pft_coms, only: SLA, q, qsw, hgt_min, include_pft, phenology
  use ed_misc_coms, only: sfilin, ied_init_mode,current_time
  use mem_sites, only: grid_res,edres
  use consts_coms, only: pio180
  use ed_misc_coms, only: use_target_year, restart_target_year,ied_init_mode,runtype
  use ed_state_vars,only: polygontype,sitetype,patchtype,edtype, &
       edgrid_g,allocate_sitetype,allocate_patchtype,allocate_polygontype
  use soil_coms, only: alloc_soilgrid
  use grid_coms,only:ngrids
  use ed_node_coms, only: mynum,nmachs,nnodetot,mchnum,machs
  use hdf5
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       globdims,chnkdims,chnkoffs
  use allometry, only: dbh2h

  implicit none
  
  character(len=1)  :: vnam
  character(len=3)  :: cgr
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
  real :: minrad, currad

  integer :: ngr,ifpy,ipft
  integer :: ipy,isi,ipa,ico
  integer :: py_index,si_index,pa_index

  ! HDF5 types are defined here
  integer :: hdferr
  include 'mpif.h'
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
  ! If this is a true history restart, the tolerance should be very small
  ! If this is an initialization, perhaps a new grid, using the closest,
  ! then there is no tolerance, ie 90.0 degrees, or some value which indicates
  ! the destination location is no-where close to the donor location

!!!  if (trim(runtype) == 'HISTORY' ) then

     ll_tolerance = (1.0/115.0)*(1.0/10.0)     ! 1/10th km
!!     ll_tolerance = (1.0/115.0)*(2.0/1.0)      ! 2km

     print*,"====================================================="
     print*,"         Entering Full State Initialization      "

!!!  else if (trim(runtype)=='INITIAL' .and. ied_init_mode==4 ) then

!!!     ll_tolerance = 20.0                       ! 20 km

!!     print*,"====================================================="
!!     print*," Entering Nearest Neighbor State File Initialization "


!!!  else
!!     call fatal_error ('Innapropriate run type encountered here'         &
 !!!         ,'init_full_history_restart','ed_history_io.f90')
!!!  end if


   
  ! at equator: (1 degree / 115 kilometers)  (1 km / 10 100-meter invervals)

  ! Open the HDF environment

  call h5open_f(hdferr)


  ! Turn off automatic error printing. This is done because there
  ! may be datasets that are not in the file, but it is OK. If
  ! data is missing that should be there, ED2 error reporting
  ! will detect it. If something is truly missing, the following
  ! call can be bypassed. Note, that automatic error reporting
  ! is turned back on at the end.
  
  !  call h5eset_auto_f(0,hdferr)


  ! Construct the file name for reinitiatlizing from
  ! The history file
  
  vnam = 'S'
  
  do ngr=1,ngrids
     
     cgrid => edgrid_g(ngr)
     

     !=======================================
     ! 1) Open the HDF5 HISTORY FILE
     !=======================================

     write(cgr,'(a1,i2.2)') 'g',ngr

     dbletime=dble(current_time%time)
     
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


     ! Only read in global grid data if this is a history restart
     ! otherwise, this data is not correct
     ! ==========================================================
     
     globdims = 0_8
     chnkdims = 0_8
     chnkoffs = 0_8
     
     globdims(1) = 1_8
     
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
     
     globdims = 0_8
     globdims(1) = int(cgrid%npolygons_global,8)
     
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
     
     globdims(1) = int(cgrid%nsites_global,8)
     
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
     
     globdims(1) = int(cgrid%npatches_global,8)
     
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
     
     globdims(1) = int(cgrid%npolygons_global,8)
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

!!!        minrad = sqrt(2*(ll_tolerance**2))

        do ifpy = 1,cgrid%npolygons_global
           
!!!           currad = sqrt( (file_lats(ifpy)-cgrid%lat(ipy))**2 + (file_lons(ifpy)-cgrid%lon(ipy))**2 )

           if ( abs(file_lats(ifpy)-cgrid%lat(ipy)) < ll_tolerance .and. &
                abs(file_lons(ifpy)-cgrid%lon(ipy)) < ll_tolerance ) then !!! .and. &
!!!                (currad <  minrad) ) then
              py_index = ifpy
!!!              minrad   = currad
           end if
           
        enddo

        if (py_index==0) then
           print*,"COULD NOT MATCH A POLYGON WITH THE DATASET"
           print*,"STOPPING"
           print*,"THIS IS THE ",ipy,"th POLYGON"
           print*,"GRID LATS: ",cgrid%lat(ipy)
           print*,"GRID LONS: ",cgrid%lon(ipy)
!           print*,"FILE LATS: ",file_lats
!           print*,"FILE LONS: ",file_lons
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

                 csite%hcapveg = 0.

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
                       
                       call fill_history_patch(cpatch,paco_id(pa_index),cgrid%ncohorts_global &
                                              ,cpoly%green_leaf_factor(:,isi))
                       
!                       do ipft = 1,n_pft
!                          csite%old_stoma_data_max(ipft,ipa)%recalc = 1
!                       enddo
                       do ico = 1,cpatch%ncohorts
                          csite%hcapveg(ipa) = csite%hcapveg(ipa) + cpatch%hcapveg(ico)
                       end do
                       

                    else

                       cpatch%ncohorts = 0
                       
                    endif
                 
                 enddo

              else

                 print*,"ATTEMPTING TO FILL SITE WITH PATCH VECTOR DATA"
                 print*,"NO PATCHES WERE FOUND in SIPA_N(SI_INDEX)"
                 print*,"THIS IS UNLIKELY AND MORALLY QUESTIONABLE."
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

  ! Turn automatic error reporting back on.
  ! This is probably unecessary, because the environment
  ! is about to be flushed.

  !  call h5eset_auto_f(1,hdferr)


  ! Close the HDF environment
  
  call h5close_f(hdferr)

  ! Initialize the disturbance transition rates
  
  write(*,'(a,i2.2)')'    Initializing anthropogenic disturbance forcing. Node: ',mynum
  call landuse_init


  return
end subroutine init_full_history_restart
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine fill_history_grid(cgrid,ipy,py_index)

  use ed_state_vars,only: edtype,polygontype
  use grid_coms,only : nzg
  use ed_max_dims,only : n_pft,n_dbh,n_age,n_dist_types
  use hdf5
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       globdims,chnkdims,chnkoffs,cnt,stride, &
       memdims,memoffs,memsize,datatype_id

  implicit none

  
#if USE_INTERF
  interface
     subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required)
       use hdf5_coms,only:memsize
       real,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
       integer :: dsetrank,iparallel
       character(len=*),intent(in) :: varn
       logical,intent(in) :: required
     end subroutine hdf_getslab_r
     subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required)
       use hdf5_coms,only:memsize
       real(kind=8),dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
       integer :: dsetrank,iparallel
       character(len=*),intent(in) :: varn
       logical,intent(in) :: required
     end subroutine hdf_getslab_d
     subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required)
       use hdf5_coms,only:memsize
       integer,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
       integer :: dsetrank,iparallel
       character(len=*),intent(in) :: varn
       logical,intent(in) :: required
     end subroutine hdf_getslab_i
  end interface
#endif

  type(edtype),target ::       cgrid

  integer,intent(in) :: ipy,py_index
  integer :: iparallel,dsetrank
  integer(SIZE_T) :: sz
  integer :: hdferr

  iparallel = 0
 

  globdims = 0_8
  chnkdims = 0_8
  chnkoffs = 0_8
  memoffs  = 0_8
  memdims  = 0_8
  memsize  = 1_8

  dsetrank = 1

  ! These are the dimensions in the filespace
  ! itself. Global is the size of the dataset,
  ! chnkoffs is the offset of the chunk we
  ! are going to read.  Chnkdims is the size
  ! of the slab that is to be read.
  
  globdims(1) = int(cgrid%npolygons_global,8)
  chnkdims(1) = 1_8
  chnkoffs(1) = int(py_index - 1,8)

  ! These are the dimensions for the memory space
  ! this should essentially be the same dimensioning
  ! as the buffer that we are filling. This routine
  ! is just filling a scalar point in a vector
  ! of polygons.

  memdims(1)  = 1_8
  memoffs(1)  = 0_8
  memsize(1)  = 1_8


  call hdf_getslab_d(cgrid%walltime_py(ipy:ipy),'WALLTIME_PY ',dsetrank,iparallel,.false.)

  call hdf_getslab_i(cgrid%lsl(ipy:ipy),'LSL ',dsetrank,iparallel,.true.)

  call hdf_getslab_r(cgrid%wbar(ipy:ipy),'WBAR ',dsetrank,iparallel,.true.)
  
  call hdf_getslab_r(cgrid%Te(ipy:ipy),'TE ',dsetrank,iparallel,.true.)

  call hdf_getslab_r(cgrid%zbar(ipy:ipy),'ZBAR ',dsetrank,iparallel,.true.)

!!  call hdf_getslab_r(cgrid%tau(ipy:ipy),'TAU ',dsetrank,iparallel,.true.)

  call hdf_getslab_r(cgrid%sheat(ipy:ipy),'SHEAT ',dsetrank,iparallel,.true.)

  call hdf_getslab_r(cgrid%baseflow(ipy:ipy),'BASEFLOW ',dsetrank,iparallel,.true.)

  call hdf_getslab_i(cgrid%load_adjacency(ipy:ipy),'LOAD_ADJACENCY ',dsetrank,iparallel,.true.)

  call hdf_getslab_r(cgrid%swliq(ipy:ipy),'SWLIQ ',dsetrank,iparallel,.true.)
  
  ! All daily and monthly variables need to be retrieved if you are loading there...
  
  if (associated(cgrid%dmean_pcpg           ))                                             &
     call hdf_getslab_r(cgrid%dmean_pcpg           (ipy:ipy) ,'DMEAN_PCPG            '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_runoff         ))                                             &
     call hdf_getslab_r(cgrid%dmean_runoff         (ipy:ipy) ,'DMEAN_RUNOFF          '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_drainage       ))                                             &
     call hdf_getslab_r(cgrid%dmean_drainage       (ipy:ipy) ,'DMEAN_DRAINAGE        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_gpp            ))                                             &
     call hdf_getslab_r(cgrid%dmean_gpp            (ipy:ipy) ,'DMEAN_GPP             '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_evap           ))                                             &
     call hdf_getslab_r(cgrid%dmean_evap           (ipy:ipy) ,'DMEAN_EVAP            '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_transp         ))                                             &
     call hdf_getslab_r(cgrid%dmean_transp         (ipy:ipy) ,'DMEAN_TRANSP          '     &
                       ,dsetrank,iparallel,.false.)                                        

  if (associated(cgrid%dmean_sensible_vc    ))                                             &
     call hdf_getslab_r(cgrid%dmean_sensible_vc    (ipy:ipy) ,'DMEAN_SENSIBLE_VC     '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_sensible_gc    ))                                             &
     call hdf_getslab_r(cgrid%dmean_sensible_gc    (ipy:ipy) ,'DMEAN_SENSIBLE_GC     '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_sensible_ac    ))                                             &
     call hdf_getslab_r(cgrid%dmean_sensible_ac    (ipy:ipy) ,'DMEAN_SENSIBLE_AC     '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_vapor_vc       ))                                             &
     call hdf_getslab_r(cgrid%dmean_vapor_vc       (ipy:ipy) ,'DMEAN_VAPOR_VC        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_vapor_gc       ))                                             &
     call hdf_getslab_r(cgrid%dmean_vapor_gc       (ipy:ipy) ,'DMEAN_VAPOR_GC        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_vapor_ac       ))                                             &
     call hdf_getslab_r(cgrid%dmean_vapor_ac       (ipy:ipy) ,'DMEAN_VAPOR_AC        '     &
                       ,dsetrank,iparallel,.false.)                                        

  if (associated(cgrid%dmean_nep            ))                                             &
     call hdf_getslab_r(cgrid%dmean_nep            (ipy:ipy) ,'DMEAN_NEP             '     &
                       ,dsetrank,iparallel,.false.)                                        

  if (associated(cgrid%dmean_plresp         ))                                             &
     call hdf_getslab_r(cgrid%dmean_plresp         (ipy:ipy) ,'DMEAN_PLRESP          '     &
                       ,dsetrank,iparallel,.false.)                                        

  if (associated(cgrid%dmean_rh             ))                                             &
     call hdf_getslab_r(cgrid%dmean_rh             (ipy:ipy) ,'DMEAN_RH              '     &
                       ,dsetrank,iparallel,.false.)                                        

  if (associated(cgrid%dmean_leaf_resp      ))                                             &
     call hdf_getslab_r(cgrid%dmean_leaf_resp      (ipy:ipy) ,'DMEAN_LEAF_RESP       '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_root_resp      ))                                             &
     call hdf_getslab_r(cgrid%dmean_root_resp      (ipy:ipy) ,'DMEAN_ROOT_RESP       '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_growth_resp    ))                                             &
     call hdf_getslab_r(cgrid%dmean_growth_resp    (ipy:ipy) ,'DMEAN_GROWTH_RESP     '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_storage_resp   ))                                             &
     call hdf_getslab_r(cgrid%dmean_storage_resp   (ipy:ipy) ,'DMEAN_STORAGE_RESP    '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_vleaf_resp     ))                                             &
     call hdf_getslab_r(cgrid%dmean_vleaf_resp     (ipy:ipy) ,'DMEAN_VLEAF_RESP      '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_fs_open        ))                                             &
     call hdf_getslab_r(cgrid%dmean_fs_open        (ipy:ipy) ,'DMEAN_FS_OPEN         '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_fsw            ))                                             &
     call hdf_getslab_r(cgrid%dmean_fsw            (ipy:ipy) ,'DMEAN_FSW             '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_fsn            ))                                             &
     call hdf_getslab_r(cgrid%dmean_fsn            (ipy:ipy) ,'DMEAN_FSN             '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_can_temp       ))                                             &
     call hdf_getslab_r(cgrid%dmean_can_temp       (ipy:ipy) ,'DMEAN_CAN_TEMP        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_can_shv        ))                                             &
     call hdf_getslab_r(cgrid%dmean_can_shv        (ipy:ipy) ,'DMEAN_CAN_SHV         '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_can_prss       ))                                             &
     call hdf_getslab_r(cgrid%dmean_can_prss       (ipy:ipy) ,'DMEAN_CAN_PRSS        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_can_theta      ))                                             &
     call hdf_getslab_r(cgrid%dmean_can_theta      (ipy:ipy) ,'DMEAN_CAN_THETA       '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_can_co2        ))                                             &
     call hdf_getslab_r(cgrid%dmean_can_co2        (ipy:ipy) ,'DMEAN_CAN_CO2         '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_can_rhos       ))                                             &
     call hdf_getslab_r(cgrid%dmean_can_rhos       (ipy:ipy) ,'DMEAN_CAN_RHOS        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_veg_energy     ))                                             &
     call hdf_getslab_r(cgrid%dmean_veg_energy     (ipy:ipy) ,'DMEAN_VEG_ENERGY      '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_veg_water      ))                                             &
     call hdf_getslab_r(cgrid%dmean_veg_water      (ipy:ipy) ,'DMEAN_VEG_WATER       '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_veg_hcap       ))                                             &
     call hdf_getslab_r(cgrid%dmean_veg_hcap       (ipy:ipy) ,'DMEAN_VEG_HCAP        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_veg_temp       ))                                             &
     call hdf_getslab_r(cgrid%dmean_veg_temp       (ipy:ipy) ,'DMEAN_VEG_TEMP        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_atm_temp       ))                                             &
     call hdf_getslab_r(cgrid%dmean_atm_temp       (ipy:ipy) ,'DMEAN_ATM_TEMP        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_atm_shv        ))                                             &
     call hdf_getslab_r(cgrid%dmean_atm_shv        (ipy:ipy) ,'DMEAN_ATM_SHV         '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_atm_prss       ))                                             &
     call hdf_getslab_r(cgrid%dmean_atm_prss       (ipy:ipy) ,'DMEAN_ATM_PRSS        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_atm_vels       ))                                             &
     call hdf_getslab_r(cgrid%dmean_atm_vels       (ipy:ipy) ,'DMEAN_ATM_VELS        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_co2_residual   ))                                             &
     call hdf_getslab_r(cgrid%dmean_co2_residual   (ipy:ipy) ,'DMEAN_CO2_RESIDUAL    '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_energy_residual))                                             &
     call hdf_getslab_r(cgrid%dmean_energy_residual(ipy:ipy) ,'DMEAN_ENERGY_RESIDUAL '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%dmean_water_residual ))                                             &
     call hdf_getslab_r(cgrid%dmean_water_residual (ipy:ipy) ,'DMEAN_WATER_RESIDUAL  '     &
                       ,dsetrank,iparallel,.false.)


  if (associated(cgrid%mmean_co2_residual   ))                                             &
     call hdf_getslab_r(cgrid%mmean_co2_residual   (ipy:ipy) ,'MMEAN_CO2_RESIDUAL    '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_energy_residual))                                             &
     call hdf_getslab_r(cgrid%mmean_energy_residual(ipy:ipy) ,'MMEAN_ENERGY_RESIDUAL '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_water_residual ))                                             &
     call hdf_getslab_r(cgrid%mmean_water_residual (ipy:ipy) ,'MMEAN_WATER_RESIDUAL  '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_gpp            ))                                             &
     call hdf_getslab_r(cgrid%mmean_gpp            (ipy:ipy) ,'MMEAN_GPP             '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_evap           ))                                             &
     call hdf_getslab_r(cgrid%mmean_evap           (ipy:ipy) ,'MMEAN_EVAP            '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_transp         ))                                             &
     call hdf_getslab_r(cgrid%mmean_transp         (ipy:ipy) ,'MMEAN_TRANSP          '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_sensible_vc    ))                                             &
     call hdf_getslab_r(cgrid%mmean_sensible_vc    (ipy:ipy) ,'MMEAN_SENSIBLE_VC     '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_sensible_gc    ))                                             &
     call hdf_getslab_r(cgrid%mmean_sensible_gc    (ipy:ipy) ,'MMEAN_SENSIBLE_GC     '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_sensible_ac    ))                                             &
     call hdf_getslab_r(cgrid%mmean_sensible_ac    (ipy:ipy) ,'MMEAN_SENSIBLE_AC     '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_vapor_vc       ))                                             &
     call hdf_getslab_r(cgrid%mmean_vapor_vc       (ipy:ipy) ,'MMEAN_VAPOR_VC        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_vapor_gc       ))                                             &
     call hdf_getslab_r(cgrid%mmean_vapor_gc       (ipy:ipy) ,'MMEAN_VAPOR_GC        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_vapor_ac       ))                                             &
     call hdf_getslab_r(cgrid%mmean_vapor_ac       (ipy:ipy) ,'MMEAN_VAPOR_AC     '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_nep            ))                                             &
     call hdf_getslab_r(cgrid%mmean_nep            (ipy:ipy) ,'MMEAN_NEP             '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_plresp         ))                                             &
     call hdf_getslab_r(cgrid%mmean_plresp         (ipy:ipy) ,'MMEAN_PLRESP          '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_rh             ))                                             &
     call hdf_getslab_r(cgrid%mmean_rh             (ipy:ipy) ,'MMEAN_RH              '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_leaf_resp      ))                                             &
     call hdf_getslab_r(cgrid%mmean_leaf_resp      (ipy:ipy) ,'MMEAN_LEAF_RESP       '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_root_resp      ))                                             &
     call hdf_getslab_r(cgrid%mmean_root_resp      (ipy:ipy) ,'MMEAN_ROOT_RESP       '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_growth_resp    ))                                             &
     call hdf_getslab_r(cgrid%mmean_growth_resp    (ipy:ipy) ,'MMEAN_GROWTH_RESP     '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_storage_resp   ))                                             &
     call hdf_getslab_r(cgrid%mmean_storage_resp   (ipy:ipy) ,'MMEAN_STORAGE_RESP    '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_vleaf_resp     ))                                             &
     call hdf_getslab_r(cgrid%mmean_vleaf_resp     (ipy:ipy) ,'MMEAN_VLEAF_RESP      '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_can_temp       ))                                             &
     call hdf_getslab_r(cgrid%mmean_can_temp       (ipy:ipy) ,'MMEAN_CAN_TEMP        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_can_shv        ))                                             &
     call hdf_getslab_r(cgrid%mmean_can_shv        (ipy:ipy) ,'MMEAN_CAN_SHV         '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_can_co2        ))                                             &
     call hdf_getslab_r(cgrid%mmean_can_co2        (ipy:ipy) ,'MMEAN_CAN_CO2         '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_can_rhos       ))                                             &
     call hdf_getslab_r(cgrid%mmean_can_rhos       (ipy:ipy) ,'MMEAN_CAN_RHOS        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%mmean_can_prss       ))                                             &
     call hdf_getslab_r(cgrid%mmean_can_prss       (ipy:ipy) ,'MMEAN_CAN_PRSS        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%mmean_can_theta      ))                                             &
     call hdf_getslab_r(cgrid%mmean_can_theta      (ipy:ipy) ,'MMEAN_CAN_THETA       '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%mmean_veg_energy     ))                                             &
     call hdf_getslab_r(cgrid%mmean_veg_energy     (ipy:ipy) ,'MMEAN_VEG_ENERGY      '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_veg_water      ))                                             &
     call hdf_getslab_r(cgrid%mmean_veg_water      (ipy:ipy) ,'MMEAN_VEG_WATER       '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_veg_temp       ))                                             &
     call hdf_getslab_r(cgrid%mmean_veg_temp       (ipy:ipy) ,'MMEAN_VEG_TEMP        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_veg_hcap       ))                                             &
     call hdf_getslab_r(cgrid%mmean_veg_hcap       (ipy:ipy) ,'MMEAN_VEG_HCAP        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_atm_temp       ))                                             &
     call hdf_getslab_r(cgrid%mmean_atm_temp       (ipy:ipy) ,'MMEAN_ATM_TEMP        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_atm_shv        ))                                             &
     call hdf_getslab_r(cgrid%mmean_atm_shv        (ipy:ipy) ,'MMEAN_ATM_SHV         '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_atm_prss       ))                                             &
     call hdf_getslab_r(cgrid%mmean_atm_prss       (ipy:ipy) ,'MMEAN_ATM_PRSS        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_atm_vels       ))                                             &
     call hdf_getslab_r(cgrid%mmean_atm_vels       (ipy:ipy) ,'MMEAN_ATM_VELS        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_pcpg           ))                                             &
     call hdf_getslab_r(cgrid%mmean_pcpg           (ipy:ipy) ,'MMEAN_PCPG            '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_runoff         ))                                             &
     call hdf_getslab_r(cgrid%mmean_runoff         (ipy:ipy) ,'MMEAN_RUNOFF          '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%mmean_drainage       ))                                             &
     call hdf_getslab_r(cgrid%mmean_drainage       (ipy:ipy) ,'MMEAN_DRAINAGE        '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%mmean_fs_open        ))                                             &
     call hdf_getslab_r(cgrid%mmean_fs_open        (ipy:ipy) ,'MMEAN_FS_OPEN         '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%mmean_fsw            ))                                             &
     call hdf_getslab_r(cgrid%mmean_fsw            (ipy:ipy) ,'MMEAN_FSW             '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%mmean_fsn            ))                                             &
     call hdf_getslab_r(cgrid%mmean_fsn            (ipy:ipy) ,'MMEAN_FSN             '     &
                       ,dsetrank,iparallel,.false.)

  if (associated(cgrid%stdev_gpp            ))                                             &
     call hdf_getslab_r(cgrid%stdev_gpp            (ipy:ipy) ,'STDEV_GPP             '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%stdev_evap           ))                                             &
     call hdf_getslab_r(cgrid%stdev_evap           (ipy:ipy) ,'STDEV_EVAP            '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%stdev_transp         ))                                             &
     call hdf_getslab_r(cgrid%stdev_transp         (ipy:ipy) ,'STDEV_TRANSP          '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%stdev_sensible       ))                                             &
     call hdf_getslab_r(cgrid%stdev_sensible       (ipy:ipy) ,'STDEV_SENSIBLE        '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%stdev_nep            ))                                             &
     call hdf_getslab_r(cgrid%stdev_nep            (ipy:ipy) ,'STDEV_NEP             '     &
                       ,dsetrank,iparallel,.false.)
 
  if (associated(cgrid%stdev_rh             ))                                             &
     call hdf_getslab_r(cgrid%stdev_rh             (ipy:ipy) ,'STDEV_RH              '     &
                       ,dsetrank,iparallel,.false.)

   ! Variables with 2 dimensions (nzg,npolygons)
   dsetrank    = 2
   globdims(1) = int(nzg,8)
   chnkdims(1) = int(nzg,8)
   memdims(1)  = int(nzg,8)
   memsize(1)  = int(nzg,8)
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8

   globdims(2)  = int(cgrid%npolygons_global,8)
   chnkdims(2)  = 1_8
   chnkoffs(2)  = int(py_index - 1,8)
   memdims(2)   = 1_8
   memsize(2)   = 1_8
   memoffs(2)   = 0_8

   call hdf_getslab_i(cgrid%ntext_soil(:,ipy)          ,'NTEXT_SOIL '       ,&
        dsetrank,iparallel,.false.)
   if(associated(cgrid%dmean_soil_temp)) &
      call hdf_getslab_r(cgrid%dmean_soil_temp(:,ipy)  ,'DMEAN_SOIL_TEMP '  ,&
      dsetrank,iparallel,.false.)
   if(associated(cgrid%dmean_soil_water)) &
      call hdf_getslab_r(cgrid%dmean_soil_water(:,ipy) ,'DMEAN_SOIL_WATER ' ,&
      dsetrank,iparallel,.false.)
   if(associated(cgrid%mmean_soil_temp)) &
      call hdf_getslab_r(cgrid%mmean_soil_temp(:,ipy)  ,'MMEAN_SOIL_TEMP '  ,&
      dsetrank,iparallel,.false.)
   if(associated(cgrid%dmean_soil_water)) &
      call hdf_getslab_r(cgrid%mmean_soil_water(:,ipy) ,'MMEAN_SOIL_WATER ' ,&
      dsetrank,iparallel,.false.)


   ! Variables with 2 dimensions (n_pft,npolygons)
   dsetrank    = 2
   globdims(1) = int(n_pft,8)
   chnkdims(1) = int(n_pft,8)
   memdims(1)  = int(n_pft,8)
   memsize(1)  = int(n_pft,8)
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8

   globdims(2)  = int(cgrid%npolygons_global,8)
   chnkdims(2)  = 1_8
   chnkoffs(2)  = int(py_index - 1,8)
   memdims(2)   = 1_8
   memsize(2)   = 1_8
   memoffs(2)   = 0_8

   if(associated(cgrid%bseeds_pft)) call hdf_getslab_r(cgrid%bseeds_pft(:,ipy) ,'BSEEDS_PFT '    , &
        dsetrank,iparallel,.false.)
   if(associated(cgrid%lai_pft)) call hdf_getslab_r(cgrid%lai_pft(:,ipy) ,'LAI_PFT '       , &
        dsetrank,iparallel,.false.)
   if(associated(cgrid%wpa_pft)) call hdf_getslab_r(cgrid%wpa_pft(:,ipy) ,'WPA_PFT '       , &
        dsetrank,iparallel,.false.)
   if(associated(cgrid%wai_pft)) call hdf_getslab_r(cgrid%wai_pft(:,ipy) ,'WAI_PFT '       , &
        dsetrank,iparallel,.false.)
   if(associated(cgrid%mmean_lai_pft)) call hdf_getslab_r(cgrid%mmean_lai_pft(:,ipy) ,'MMEAN_LAI_PFT ' , &
        dsetrank,iparallel,.false.)
   if(associated(cgrid%mmean_wpa_pft)) call hdf_getslab_r(cgrid%mmean_wpa_pft(:,ipy) ,'MMEAN_WPA_PFT ' , &
        dsetrank,iparallel,.false.)
   if(associated(cgrid%mmean_wai_pft)) call hdf_getslab_r(cgrid%mmean_wai_pft(:,ipy) ,'MMEAN_WAI_PFT ' , &
        dsetrank,iparallel,.false.)
   if(associated(cgrid%agb_pft)) call hdf_getslab_r(cgrid%agb_pft(:,ipy) ,'AGB_PFT '       , &
        dsetrank,iparallel,.true.)
   if(associated(cgrid%ba_pft)) call hdf_getslab_r(cgrid%ba_pft(:,ipy) ,'BA_PFT '        ,   &
        dsetrank,iparallel,.true.)


   ! Variables with 2 dimensions (n_dbh,npolygons)
   dsetrank    = 2
   globdims(1) = int(n_dbh,8)
   chnkdims(1) = int(n_dbh,8)
   memdims(1)  = int(n_dbh,8)
   memsize(1)  = int(n_dbh,8)
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8

   globdims(2)  = int(cgrid%npolygons_global,8)
   chnkdims(2)  = 1_8
   chnkoffs(2)  = int(py_index - 1,8)
   memdims(2)   = 1_8
   memsize(2)   = 1_8
   memoffs(2)   = 0_8

   if(associated(cgrid%dmean_gpp_dbh)) call hdf_getslab_r(cgrid%dmean_gpp_dbh(:,ipy) , &
        'DMEAN_GPP_DBH ' ,dsetrank,iparallel,.false.)
   if(associated(cgrid%mmean_gpp_dbh)) call hdf_getslab_r(cgrid%mmean_gpp_dbh(:,ipy) , &
        'MMEAN_GPP_DBH ' ,dsetrank,iparallel,.false.)

   ! Variables with 2 dimensions (13,npolygons)
   dsetrank    = 2
   globdims(1) = int(13,8)
   chnkdims(1) = int(13,8)
   memdims(1)  = int(13,8)
   memsize(1)  = int(13,8)
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8

   globdims(2)  = int(cgrid%npolygons_global,8)
   chnkdims(2)  = 1_8
   chnkoffs(2)  = int(py_index - 1,8)
   memdims(2)   = 1_8
   memsize(2)   = 1_8
   memoffs(2)   = 0_8

   if(associated(cgrid%workload)) call hdf_getslab_r(cgrid%workload(:,ipy) , &
        'WORKLOAD ' ,dsetrank,iparallel,.true.)


   ! Variables with three dimensions(n_dist_types,n_dist_types,npolygons)
   dsetrank    = 3
   globdims(1) = int(n_dist_types,8)
   chnkdims(1) = int(n_dist_types,8)
   memdims(1)  = int(n_dist_types,8)
   memsize(1)  = int(n_dist_types,8)
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8

   globdims(2) = int(n_dist_types,8)
   chnkdims(2) = int(n_dist_types,8)
   memdims(2)  = int(n_dist_types,8)
   memsize(2)  = int(n_dist_types,8)
   chnkoffs(2) = 0_8
   memoffs(2)  = 0_8

   globdims(3)  = int(cgrid%npolygons_global,8)
   chnkdims(3)  = 1_8
   chnkoffs(3)  = int(py_index - 1,8)
   memdims(3)   = 1_8
   memsize(3)   = 1_8
   memoffs(3)   = 0_8
   if (associated(cgrid%disturbance_rates))                                                &
       call hdf_getslab_r(cgrid%disturbance_rates(:,:,ipy),'DISTURBANCE_RATES '            &
                         ,dsetrank,iparallel,.false.)

   return
 end subroutine fill_history_grid
 !=========================================================================================!
 !=========================================================================================!






 !==========================================================================================!
 !==========================================================================================!
 subroutine fill_history_polygon(cpoly,pysi_index,nsites_global)

   use ed_state_vars,only: polygontype
   use hdf5
   use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
        globdims,chnkdims,chnkoffs,cnt,stride, &
        memdims,memoffs,memsize

   use grid_coms,only : nzg
   use ed_max_dims,only : n_pft,n_dbh,n_dist_types

   implicit none

#if USE_INTERF
     interface
      subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required)
        use hdf5_coms,only:memsize
        real,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
        integer :: dsetrank,iparallel
        character(len=*),intent(in) :: varn
        logical,intent(in):: required
      end subroutine hdf_getslab_r
      subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required)
        use hdf5_coms,only:memsize
        real(kind=8),dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
        integer :: dsetrank,iparallel
        character(len=*),intent(in) :: varn
        logical,intent(in) :: required
      end subroutine hdf_getslab_d
      subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required)
        use hdf5_coms,only:memsize
        integer,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
        integer :: dsetrank,iparallel
        character(len=*),intent(in) :: varn
        logical,intent(in) :: required
      end subroutine hdf_getslab_i
   end interface
#endif

   type(polygontype),target :: cpoly
   integer,intent(in) :: pysi_index
   integer,intent(in) :: nsites_global
   integer :: iparallel
   integer :: dsetrank

   iparallel = 0

   dsetrank = 1_8
   globdims = 0_8
   chnkdims = 0_8
   chnkoffs = 0_8
   memoffs  = 0_8
   memdims  = 0_8
   memsize  = 1_8

   globdims(1) = int(nsites_global,8)
   chnkdims(1) = int(cpoly%nsites,8)
   chnkoffs(1) = int(pysi_index - 1,8)
   memdims(1)  = int(cpoly%nsites,8)
   memsize(1)  = int(cpoly%nsites,8)
   memoffs(1)  = 0_8

   call hdf_getslab_i(cpoly%patch_count,'PATCH_COUNT ',dsetrank,iparallel,.true.)  
   call hdf_getslab_i(cpoly%sitenum,'SITENUM ',dsetrank,iparallel,.true.)

   call hdf_getslab_i(cpoly%lsl,'LSL_SI ',dsetrank,iparallel,.true.)   
   call hdf_getslab_r(cpoly%area,'AREA_SI ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%patch_area,'PATCH_AREA ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%elevation,'ELEVATION ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%slope,'SLOPE ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%aspect,'ASPECT ',dsetrank,iparallel,.true.)

   call hdf_getslab_i(cpoly%num_landuse_years,'NUM_LANDUSE_YEARS ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%TCI,'TCI ',dsetrank,iparallel,.true.)      
   call hdf_getslab_r(cpoly%pptweight,'pptweight ',dsetrank,iparallel,.true.)      
   call hdf_getslab_i(cpoly%hydro_next,'HYDRO_NEXT ',dsetrank,iparallel,.true.)
   call hdf_getslab_i(cpoly%hydro_prev,'HYDRO_PREV ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%moist_W,'MOIST_W ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%moist_f,'MOIST_F ',dsetrank,iparallel,.true.)  
   call hdf_getslab_r(cpoly%moist_tau,'MOIST_TAU ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%moist_zi,'MOIST_ZI ',dsetrank,iparallel,.true.) 
   call hdf_getslab_r(cpoly%baseflow,'BASEFLOW_SI ',dsetrank,iparallel,.true.) 
!   call hdf_getslab_i(cpoly%metplex_beg_month,'METPLEX_BEG_MONTH ',dsetrank,iparallel,.true.)
!   call hdf_getslab_i(cpoly%metplex_beg_year,'METPLEX_BEG_YEAR ',dsetrank,iparallel,.true.)
!   call hdf_getslab_i(cpoly%metplex_end_year,'METPLEX_END_YEAR ',dsetrank,iparallel,.true.)

   call hdf_getslab_r(cpoly%min_monthly_temp,'MIN_MONTHLY_TEMP ',dsetrank,iparallel,.true.)
!   call hdf_getslab_r(cpoly%removed_biomass,'REMOVED_BIOMASS ',dsetrank,iparallel,.true.) 
!   call hdf_getslab_r(cpoly%harvested_biomass,'HARVESTED_BIOMASS ', &
!        dsetrank,iparallel,.true.) 
   call hdf_getslab_i(cpoly%plantation,'PLANTATION_SI ',dsetrank,iparallel,.true.) 
   call hdf_getslab_i(cpoly%agri_stocking_pft,'AGRI_STOCKING_PFT ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%agri_stocking_density,'AGRI_STOCKING_DENSITY ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_i(cpoly%plantation_stocking_pft,'PLANTATION_STOCKING_PFT ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%plantation_stocking_density,'PLANTATION_STOCKING_DENSITY ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%primary_harvest_memory,'PRIMARY_HARVEST_MEMORY ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%secondary_harvest_memory,'SECONDARY_HARVEST_MEMORY ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%fire_disturbance_rate,'FIRE_DISTURBANCE_RATE ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%ignition_rate,'IGNITION_RATE ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%treefall_disturbance_rate,'TREEFALL_DISTURBANCE_RATE ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%nat_disturbance_rate,'NAT_DISTURBANCE_RATE ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_i(cpoly%nat_dist_type,'NAT_DIST_TYPE ',dsetrank,iparallel,.true.)

   dsetrank    = 2_8
   globdims(1) = int(n_pft,8)
   chnkdims(1) = int(n_pft,8)
   memdims(1)  = int(n_pft,8)
   memsize(1)  = int(n_pft,8)
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8
   globdims(2)  = int(nsites_global,8)
   chnkdims(2)  = int(cpoly%nsites,8)
   chnkoffs(2)  = int(pysi_index - 1,8)
   memdims(2)   = int(cpoly%nsites,8)
   memsize(2)   = int(cpoly%nsites,8)
   memoffs(2)   = 0_8

   if (associated(cpoly%lai_pft)) call hdf_getslab_r(cpoly%lai_pft,'LAI_PFT_SI ', &
        dsetrank,iparallel,.false.)
   if (associated(cpoly%wpa_pft)) call hdf_getslab_r(cpoly%wpa_pft,'WPA_PFT_SI ', &
        dsetrank,iparallel,.false.)
   if (associated(cpoly%wai_pft)) call hdf_getslab_r(cpoly%wai_pft,'WAI_PFT_SI ', &
        dsetrank,iparallel,.false.)
   call hdf_getslab_r(cpoly%green_leaf_factor,'GREEN_LEAF_FACTOR ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%leaf_aging_factor,'LEAF_AGING_FACTOR ', &
        dsetrank,iparallel,.true.)

   dsetrank    = 2_8
   globdims(1) = int(nzg,8)
   chnkdims(1) = int(nzg,8)
   memdims(1)  = int(nzg,8)
   memsize(1)  = int(nzg,8)
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8
   globdims(2)  = int(nsites_global,8)
   chnkdims(2)  = int(cpoly%nsites,8)
   chnkoffs(2)  = int(pysi_index - 1,8)
   memdims(2)   = int(cpoly%nsites,8)
   memsize(2)   = int(cpoly%nsites,8)
   memoffs(2)   = 0_8

   call hdf_getslab_i(cpoly%ntext_soil,'NTEXT_SOIL_SI ',dsetrank,iparallel,.true.)

   dsetrank    = 2_8
   globdims(1) = 12_8
   chnkdims(1) = 12_8
   memdims(1)  = 12_8
   memsize(1)  = 12_8
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8
   globdims(2)  = int(nsites_global,8)
   chnkdims(2)  = int(cpoly%nsites,8)
   chnkoffs(2)  = int(pysi_index - 1,8)
   memdims(2)   = int(cpoly%nsites,8)
   memsize(2)   = int(cpoly%nsites,8)
   memoffs(2)   = 0_8

   call hdf_getslab_r(cpoly%lambda_fire,'LAMBDA_FIRE ',dsetrank,iparallel,.true.)

   dsetrank    = 2_8
   globdims(1) = int(n_dist_types,8)
   chnkdims(1) = int(n_dist_types,8)
   memdims(1)  = int(n_dist_types,8)
   memsize(1)  = int(n_dist_types,8)
   chnkoffs(1) = 0_8
   memoffs(1)  = 0_8
   globdims(2)  = int(nsites_global,8)
   chnkdims(2)  = int(cpoly%nsites,8)
   chnkoffs(2)  = int(pysi_index - 1,8)
   memdims(2)   = int(cpoly%nsites,8)
   memsize(2)   = int(cpoly%nsites,8)
   memoffs(2)   = 0_8

   call hdf_getslab_r(cpoly%loss_fraction,'LOSS_FRACTION ',dsetrank,iparallel,.true.)

   dsetrank    = 3_8
   globdims(1:2) = int(n_dist_types,8)
   chnkdims(1:2) = int(n_dist_types,8)
   memdims(1:2)  = int(n_dist_types,8)
   memsize(1:2)  = int(n_dist_types,8)
   chnkoffs(1:2) = 0
   memoffs(1:2)  = 0
   globdims(3)  = int(nsites_global,8)
   chnkdims(3)  = int(cpoly%nsites,8)
   chnkoffs(3)  = int(pysi_index - 1,8)
   memdims(3)   = int(cpoly%nsites,8)
   memsize(3)   = int(cpoly%nsites,8)
   memoffs(3)   = 0

   call hdf_getslab_r(cpoly%disturbance_memory,'DISTURBANCE_MEMORY ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%disturbance_rates,'DISTURBANCE_RATES_SI ', &
        dsetrank,iparallel,.true.)

   dsetrank    = 3
   globdims(3) = int(nsites_global,8)
   chnkdims(3) = int(cpoly%nsites,8)
   chnkoffs(3) = int(pysi_index - 1,8)
   memdims(3)  = int(cpoly%nsites,8)
   memsize(3)  = int(cpoly%nsites,8)
   memoffs(3)  = 0
   globdims(2) = int(n_dbh,8)
   chnkdims(2) = int(n_dbh,8)
   memdims(2)  = int(n_dbh,8)
   memsize(2)  = int(n_dbh,8)
   chnkoffs(2) = 0
   memoffs(2)  = 0
   globdims(1) = int(n_pft,8)
   chnkdims(1) = int(n_pft,8)
   memdims(1)  = int(n_pft,8)
   memsize(1)  = int(n_pft,8)
   chnkoffs(1) = 0
   memoffs(1)  = 0

   call hdf_getslab_r(cpoly%basal_area,'BASAL_AREA_SI ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%agb,'AGB_SI ',dsetrank,iparallel,.true.)
!   call hdf_getslab_r(cpoly%pldens,'PLDENS_SI ',dsetrank,iparallel,.false.)
   call hdf_getslab_r(cpoly%basal_area_growth,'BASAL_AREA_GROWTH ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%agb_growth,'AGB_GROWTH ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%basal_area_mort,'BASAL_AREA_MORT ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%basal_area_cut,'BASAL_AREA_CUT ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%agb_mort,'AGB_MORT ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(cpoly%agb_cut,'AGB_CUT ',dsetrank,iparallel,.true.)


   return
 end subroutine fill_history_polygon
 !==========================================================================================!
 !==========================================================================================!






 !==========================================================================================!
 !==========================================================================================!
 subroutine fill_history_site(csite,sipa_index,npatches_global)

   use ed_state_vars,only: sitetype
   use grid_coms,only : nzg,nzs
   use c34constants,only:n_stoma_atts
   use ed_max_dims,only : n_pft,n_dbh
   use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
        globdims,chnkdims,chnkoffs,cnt,stride, &
        memdims,memoffs,memsize,datatype_id,setsize
   use fusion_fission_coms, only: ff_ndbh
   use hdf5

   implicit none

#if USE_INTERF
     interface
      subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required)
        use hdf5_coms,only:memsize
        real,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
        integer :: dsetrank,iparallel
        character(len=*),intent(in) :: varn
        logical,intent(in) :: required
      end subroutine hdf_getslab_r
      subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required)
        use hdf5_coms,only:memsize
        real(kind=8),dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
        integer :: dsetrank,iparallel
        character(len=*),intent(in) :: varn
        logical,intent(in) :: required
      end subroutine hdf_getslab_d
      subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required)
        use hdf5_coms,only:memsize
        integer,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
        integer :: dsetrank,iparallel
        character(len=*),intent(in) :: varn
        logical,intent(in) :: required
      end subroutine hdf_getslab_i
   end interface
#endif

   type(sitetype),target :: csite
   integer,intent(in) :: sipa_index
   integer,intent(in) :: npatches_global
   integer :: iparallel
   integer :: dsetrank
   integer :: hdferr
   integer :: ipa,ipft
   real(kind=8),allocatable, dimension(:,:) ::  buff

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

   globdims(1) = int(npatches_global,8)
   chnkdims(1) = int(csite%npatches,8)
   chnkoffs(1) = int(sipa_index - 1,8)

   memdims(1)  = int(csite%npatches,8)
   memsize(1)  = int(csite%npatches,8)
   memoffs(1)  = 0

   call hdf_getslab_i(csite%dist_type,'DIST_TYPE ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%age,'AGE ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%area,'AREA ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%fast_soil_C,'FAST_SOIL_C ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%slow_soil_C,'SLOW_SOIL_C ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%structural_soil_C,'STRUCTURAL_SOIL_C ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%structural_soil_L,'STRUCTURAL_SOIL_L ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%mineralized_soil_N,'MINERALIZED_SOIL_N ', &
        dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%fast_soil_N,'FAST_SOIL_N ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%sum_dgd,'SUM_DGD ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%sum_chd,'SUM_CHD ',dsetrank,iparallel,.true.)
   call hdf_getslab_i(csite%plantation,'PLANTATION ',dsetrank,iparallel,.true.)
   !  call hdf_getslab_i(csite%cohort_count,'COHORT_COUNT ',dsetrank,iparallel)
   call hdf_getslab_r(csite%can_enthalpy,'CAN_ENTHALPY ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%can_prss,'CAN_PRSS ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%can_theta,'CAN_THETA ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%can_temp,'CAN_TEMP ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%can_shv,'CAN_SHV ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%can_co2,'CAN_CO2 ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%can_rhos,'CAN_RHOS ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%can_depth,'CAN_DEPTH ',dsetrank,iparallel,.true.)
   !  call hdf_getslab_i(csite%pname,'PNAME ',dsetrank,iparallel)
   call hdf_getslab_r(csite%lai,'LAI_PA ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%wpa,'WPA_PA ',dsetrank,iparallel,.false.)
   call hdf_getslab_r(csite%wai,'WAI_PA ',dsetrank,iparallel,.false.)
   call hdf_getslab_i(csite%nlev_sfcwater,'NLEV_SFCWATER ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ground_shv,'GROUND_SHV ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%surface_ssh,'SURFACE_SSH ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%rough,'ROUGH ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%avg_daily_temp,'AVG_DAILY_TEMP ',dsetrank,iparallel,.true.)  
   call hdf_getslab_r(csite%mean_rh,'MEAN_RH ',dsetrank,iparallel,.true.)

   if (associated(csite%dmean_rh       )) &
        call hdf_getslab_r(csite%dmean_rh,'DMEAN_RH_PA ',dsetrank,iparallel,.false.)
   if (associated(csite%mmean_rh       )) &
        call hdf_getslab_r(csite%mmean_rh,'MMEAN_RH_PA ',dsetrank,iparallel,.false.)

   call hdf_getslab_r(csite%lambda_light,'LAMBDA_LIGHT ',dsetrank,iparallel,.true.)

   if (associated(csite%dmean_lambda_light       )) &
        call hdf_getslab_r(csite%dmean_lambda_light,'DMEAN_LAMBDA_LIGHT ',dsetrank,iparallel,.false.)

   if (associated(csite%mmean_lambda_light       )) &
        call hdf_getslab_r(csite%mmean_lambda_light,'MMEAN_LAMBDA_LIGHT ',dsetrank,iparallel,.false.)
   
   call hdf_getslab_r(csite%mean_nep,'MEAN_NEP ',dsetrank,iparallel,.true.)

   call hdf_getslab_r(csite%wbudget_loss2atm,'WBUDGET_LOSS2ATM ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%wbudget_denseffect,'WBUDGET_DENSEFFECT ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%wbudget_precipgain,'WBUDGET_PRECIPGAIN ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%wbudget_loss2runoff,'WBUDGET_LOSS2RUNOFF ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%wbudget_initialstorage,'WBUDGET_INITIALSTORAGE ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ebudget_loss2atm,'EBUDGET_LOSS2ATM ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ebudget_denseffect,'EBUDGET_DENSEFFECT ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ebudget_loss2runoff,'EBUDGET_LOSS2RUNOFF ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ebudget_netrad,'EBUDGET_NETRAD ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ebudget_latent,'EBUDGET_LATENT ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ebudget_precipgain,'EBUDGET_PRECIPGAIN ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ebudget_initialstorage,'EBUDGET_INITIALSTORAGE ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%co2budget_initialstorage,'CO2BUDGET_INITIALSTORAGE ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%co2budget_loss2atm,'CO2BUDGET_LOSS2ATM ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%co2budget_denseffect,'CO2BUDGET_DENSEFFECT ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%co2budget_gpp,'CO2BUDGET_GPP ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%co2budget_plresp,'CO2BUDGET_PLRESP ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%co2budget_rh,'CO2BUDGET_RH ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%today_A_decomp,'TODAY_A_DECOMP ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%today_Af_decomp,'TODAY_AF_DECOMP ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%dmean_A_decomp,'DMEAN_A_DECOMP ',dsetrank,iparallel,.false.)
   call hdf_getslab_r(csite%dmean_Af_decomp,'DMEAN_AF_DECOMP ',dsetrank,iparallel,.false.)
   call hdf_getslab_r(csite%mmean_A_decomp,'MMEAN_A_DECOMP ',dsetrank,iparallel,.false.)
   call hdf_getslab_r(csite%mmean_Af_decomp,'MMEAN_AF_DECOMP ',dsetrank,iparallel,.false.)
   call hdf_getslab_r(csite%veg_rough,'VEG_ROUGH ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%veg_height ,'VEG_HEIGHT ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%fsc_in,'FSC_IN ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ssc_in,'SSC_IN ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%ssl_in,'SSL_IN ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%fsn_in,'FSN_IN ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%total_plant_nitrogen_uptake,'TOTAL_PLANT_NITROGEN_UPTAKE ',dsetrank,iparallel,.true.)
   
   !  call hdf_getslab_r(csite%rshort_g,'RSHORT_G ',dsetrank,iparallel,.true.)
   !  call hdf_getslab_r(csite%rshort_g_beam,'RSHORT_G_BEAM ',dsetrank,iparallel,.true.)
   !  call hdf_getslab_r(csite%rshort_g_diffuse,'RSHORT_G_DIFFUSE ',dsetrank,iparallel,.true.)
   !  call hdf_getslab_r(csite%rlong_g,'RLONG_G ',dsetrank,iparallel,.true.)
   !  call hdf_getslab_r(csite%rlong_g_surf,'RLONG_G_SURF ',dsetrank,iparallel)
   !  call hdf_getslab_r(csite%rlong_g_incid,'RLONG_G_INCID ',dsetrank,iparallel)
   !  call hdf_getslab_r(csite%rlong_s,'RLONG_S ',dsetrank,iparallel)
   !  call hdf_getslab_r(csite%rlong_s_surf,'RLONG_S_SURF ',dsetrank,iparallel)
   !  call hdf_getslab_r(csite%rlong_s_incid,'RLONG_S_INCID ',dsetrank,iparallel)
   
   !  call hdf_getslab_r(csite%albedt,'ALBEDT ',dsetrank,iparallel)
   !  call hdf_getslab_r(csite%albedo_beam,'ALBEDO_BEAM ',dsetrank,iparallel)
   !  call hdf_getslab_r(csite%albedo_diffuse,'ALBEDO_DIFFUSE ',dsetrank,iparallel)
   !  call hdf_getslab_r(csite%rlongup,'RLONGUP ',dsetrank,iparallel)
   !  call hdf_getslab_r(csite%rlong_albedo,'RLONGUP_ALBEDO ',dsetrank,iparallel)
   call hdf_getslab_r(csite%total_snow_depth,'TOTAL_SNOW_DEPTH ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%snowfac,'SNOWFAC ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%A_decomp,'A_DECOMP ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%f_decomp,'F_DECOMP ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%rh,'RH ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%cwd_rh,'CWD_RH ',dsetrank,iparallel,.true.)
   call hdf_getslab_i(csite%fuse_flag,'FUSE_FLAG ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%plant_ag_biomass,'PLANT_AG_BIOMASS ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%mean_wflux,'MEAN_WFLUX ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%mean_latflux,'MEAN_LATFLUX ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%mean_hflux,'MEAN_HFLUX ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%mean_runoff,'MEAN_RUNOFF ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%mean_qrunoff,'MEAN_QRUNOFF ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%htry,'HTRY ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%dmean_rk4step,'DMEAN_RK4STEP ',dsetrank,iparallel,.false.)
   call hdf_getslab_r(csite%mmean_rk4step,'MMEAN_RK4STEP ',dsetrank,iparallel,.false.)
   
   dsetrank    = 2
   globdims(1) = int(nzs,8)
   chnkdims(1) = int(nzs,8)
   memdims(1)  = int(nzs,8)
   memsize(1)  = int(nzs,8)
   chnkoffs(1) = 0
   memoffs(1)  = 0
   globdims(2) = int(npatches_global,8)
   chnkdims(2) = int(csite%npatches,8)
   chnkoffs(2) = int(sipa_index - 1,8)
   memdims(2)  = int(csite%npatches,8)
   memsize(2)  = int(csite%npatches,8)
   memoffs(2)  = 0
   
   call hdf_getslab_r(csite%sfcwater_mass,'SFCWATER_MASS ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%sfcwater_energy,'SFCWATER_ENERGY ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%sfcwater_depth,'SFCWATER_DEPTH ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%rshort_s,'RSHORT_S ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%rshort_s_beam,'RSHORT_S_BEAM ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%rshort_s_diffuse,'RSHORT_S_DIFFUSE ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%sfcwater_tempk,'SFCWATER_TEMPK ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%sfcwater_fracliq,'SFCWATER_FRACLIQ ',dsetrank,iparallel,.true.)
   
   dsetrank    = 2
   globdims(1) = int(nzg,8)
   chnkdims(1) = int(nzg,8)
   memdims(1)  = int(nzg,8)
   memsize(1)  = int(nzg,8)
   chnkoffs(1) = 0
   memoffs(1)  = 0
   globdims(2) = int(npatches_global,8)
   chnkdims(2) = int(csite%npatches,8)
   chnkoffs(2) = int(sipa_index - 1,8)
   memdims(2)  = int(csite%npatches,8)
   memsize(2)  = int(csite%npatches,8)
   memoffs(2)  = 0
   
   call hdf_getslab_i(csite%ntext_soil,'NTEXT_SOIL_PA ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%soil_energy,'SOIL_ENERGY_PA ',dsetrank,iparallel,.true.)
   call hdf_getslab_r(csite%soil_water,'SOIL_WATER_PA ',dsetrank,iparallel,.true.)

   !-----------------------------------------------------------------------------------!
   !  Soil water is double precision, although it may not be DP in the dataset
   !  The following lines make provisions for this by testing the dataset.
   
!   call h5dopen_f(file_id,'SOIL_WATER_PA ', dset_id, hdferr)
!   if (hdferr /= 0 ) then
!      call fatal_error('Dataset did not have soil water?' &
!           ,'fill_history_site','ed_history_io.f90')
!   endif
   ! ---------------------------------------------------------------------------------!
   ! THESE LINES ARE USEFULL FOR DETERMINING DATA SIZE OF ANY GIVEN OBJECT IN A SET   !
   ! ---------------------------------------------------------------------------------!
   
   !call h5dget_type_f(dset_id,datatype_id,hdferr)
   !call h5tget_size_f(datatype_id,setsize,hdferr)
   !call h5dclose_f(dset_id  , hdferr)
 
! =============================================================================================
! KEEP THIS CODE AS A TEMPLATE IN CASE WE NEED TO DO SOMETHING LIKE THIS IN THE FUTURE
!  HELPFUL IF WE CHANGE DATA TYPES
! ---------------------------------------------------------------------------------------------
!   if (setsize==4_8) then  !Old precision
!      call hdf_getslab_r(csite%soil_water,'SOIL_WATER_PA ',dsetrank,iparallel,.true.)
!   else if (setsize==8_8) then ! Newer precision
!      allocate(buff(nzg,csite%npatches))
!      write (unit=*,fmt='(a)') '-------------------------------------------------------------------'
!      write (unit=*,fmt='(a)') '  Loading 8-byte precision soil water and converting to 4-byte'
!      write (unit=*,fmt='(a)') '-------------------------------------------------------------------'
!      call hdf_getslab_d(buff,'SOIL_WATER_PA ',dsetrank,iparallel,.true.)
!      csite%soil_water(1:nzg,1:csite%npatches) = sngl(buff(1:nzg,1:csite%npatches))
!      deallocate(buff)
!  else
!     call fatal_error('Soil water dataset is not real nor double?'                         &
!                     ,'fill_history_site','ed_history_io.f90')
!  end if

  !--------------------------------------------------------------------------------------------  
  
  call hdf_getslab_r(csite%soil_tempk,'SOIL_TEMPK_PA ',dsetrank,iparallel,.true.)
  call hdf_getslab_r(csite%soil_fracliq,'SOIL_FRACLIQ_PA ',dsetrank,iparallel,.true.)

  dsetrank    = 2
  globdims(1) = int(n_pft,8)
  chnkdims(1) = int(n_pft,8)
  memdims(1)  = int(n_pft,8)
  memsize(1)  = int(n_pft,8)
  chnkoffs(1) = 0_8
  memoffs(1)  = 0_8
  globdims(2) = int(npatches_global,8)
  chnkdims(2) = int(csite%npatches,8)
  chnkoffs(2) = int(sipa_index - 1,8)
  memdims(2)  = int(csite%npatches,8)
  memsize(2)  = int(csite%npatches,8)
  memoffs(2)  = 0_8
  
  call hdf_getslab_r(csite%A_o_max,'A_O_MAX ',dsetrank,iparallel,.true.) 
  call hdf_getslab_r(csite%A_c_max,'A_C_MAX ',dsetrank,iparallel,.true.) 
  call hdf_getslab_r(csite%repro,'REPRO_PA ',dsetrank,iparallel,.true.)


  dsetrank    = 2
  globdims(1) = int(n_dbh,8)
  chnkdims(1) = int(n_dbh,8)
  memdims(1)  = int(n_dbh,8)
  memsize(1)  = int(n_dbh,8)
  chnkoffs(1) = 0_8
  memoffs(1)  = 0_8
  globdims(2) = int(npatches_global,8)
  chnkdims(2) = int(csite%npatches,8)
  chnkoffs(2) = int(sipa_index - 1,8)
  memdims(2)  = int(csite%npatches,8)
  memsize(2)  = int(csite%npatches,8)
  memoffs(2)  = 0_8
  call hdf_getslab_r(csite%co2budget_gpp_dbh,'CO2BUDGET_GPP_DBH ',dsetrank,iparallel,.true.)

!!!! MAY NEED TO ADD THIS ONE
!  call hdf_getslab_r(csite%old_stoma_data_max,'OLD ',dsetrank,iparallel,.true.)

  dsetrank    = 3
  globdims(3) = int(npatches_global,8)
  chnkdims(3) = int(csite%npatches,8)
  chnkoffs(3) = int(sipa_index - 1,8)

  memdims(3)  = int(csite%npatches,8)
  memsize(3)  = int(csite%npatches,8)
  memoffs(3)  = 0_8
  
  globdims(2) = int(ff_ndbh,8)
  chnkdims(2) = int(ff_ndbh,8)
  memdims(2)  = int(ff_ndbh,8)
  memsize(2)  = int(ff_ndbh,8)
  chnkoffs(2) = 0_8
  memoffs(2)  = 0_8

  globdims(1) = int(n_pft,8)
  chnkdims(1) = int(n_pft,8)
  memdims(1)  = int(n_pft,8)
  memsize(1)  = int(n_pft,8)
  chnkoffs(1) = 0_8
  memoffs(1)  = 0_8

  call hdf_getslab_r(csite%pft_density_profile,'PFT_DENSITY_PROFILE ',dsetrank,iparallel,.true.)


  dsetrank    = 3
  globdims(3) = int(npatches_global,8)
  chnkdims(3) = int(csite%npatches,8)
  chnkoffs(3) = int(sipa_index - 1,8)

  memdims(3)  = int(csite%npatches,8)
  memsize(3)  = int(csite%npatches,8)
  memoffs(3)  = 0_8
  
  globdims(2) = int(n_pft,8)
  chnkdims(2) = int(n_pft,8)
  memdims(2)  = int(n_pft,8)
  memsize(2)  = int(n_pft,8)
  chnkoffs(2) = 0_8
  memoffs(2)  = 0_8

  globdims(1) = int(n_stoma_atts,8)
  chnkdims(1) = int(n_stoma_atts,8)
  memdims(1)  = int(n_stoma_atts,8)
  memsize(1)  = int(n_stoma_atts,8)
  chnkoffs(1) = 0_8
  memoffs(1)  = 0_8

  call hdf_getslab_r(csite%old_stoma_vector_max,'OLD_STOMA_VECTOR_MAX ',dsetrank,iparallel,.true.)

  patchloop: do ipa=1,csite%npatches
     pftloop: do ipft = 1,n_pft
        csite%old_stoma_data_max(ipft,ipa)%recalc = int(csite%old_stoma_vector_max(1,ipft,ipa))
        csite%old_stoma_data_max(ipft,ipa)%T_L    = csite%old_stoma_vector_max(2,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%e_A    = csite%old_stoma_vector_max(3,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%PAR    = csite%old_stoma_vector_max(4,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%rb_factor = csite%old_stoma_vector_max(5,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%prss   = csite%old_stoma_vector_max(6,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%phenology_factor = csite%old_stoma_vector_max(7,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%gsw_open = csite%old_stoma_vector_max(8,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%ilimit = int(csite%old_stoma_vector_max(9,ipft,ipa))
        
        csite%old_stoma_data_max(ipft,ipa)%T_L_residual = csite%old_stoma_vector_max(10,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%e_a_residual = csite%old_stoma_vector_max(11,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%par_residual = csite%old_stoma_vector_max(12,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%rb_residual  = csite%old_stoma_vector_max(13,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%prss_residual= csite%old_stoma_vector_max(14,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%leaf_residual= csite%old_stoma_vector_max(15,ipft,ipa)
        csite%old_stoma_data_max(ipft,ipa)%gsw_residual = csite%old_stoma_vector_max(16,ipft,ipa)
     end do pftloop
  end do patchloop

  return
end subroutine fill_history_site

!==========================================================================================!
!==========================================================================================!

subroutine fill_history_patch(cpatch,paco_index,ncohorts_global,green_leaf_factor)
  
  use ed_state_vars,only: patchtype
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       filespace,memspace, &
       globdims,chnkdims,chnkoffs,cnt,stride, &
       memdims,memoffs,memsize
  use consts_coms, only: cliq,cice,t3ple,tsupercool
  use c34constants,only: n_stoma_atts
  use ed_max_dims,only: n_pft, n_mort
  use ed_therm_lib, only : calc_hcapveg
  use allometry, only : area_indices
  use therm_lib, only : qwtk
  implicit none

#if USE_INTERF
  interface
     subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required)
       use hdf5_coms,only:memsize
       real,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
       integer :: dsetrank,iparallel
       character(len=*),intent(in) :: varn
       logical,intent(in) :: required
     end subroutine hdf_getslab_r
     subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required)
       use hdf5_coms,only:memsize
       real(kind=8),dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
       integer :: dsetrank,iparallel
       character(len=*),intent(in) :: varn
       logical,intent(in) :: required
     end subroutine hdf_getslab_d
     subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required)
       use hdf5_coms,only:memsize
       integer,dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff
       integer :: dsetrank,iparallel
       character(len=*),intent(in) :: varn
       logical,intent(in) :: required
     end subroutine hdf_getslab_i
  end interface
#endif

  type(patchtype),target :: cpatch
  integer,intent(in) :: paco_index
  integer,intent(in) :: ncohorts_global
  real, dimension(n_pft), intent(in) :: green_leaf_factor
  integer :: iparallel,dsetrank
  
  ! Needed for reconstructing veg_energy if using an old restart
  ! ------------------------------------------------------------
  real :: plai
  integer :: ico
  ! ------------------------------------------------------------

  iparallel = 0
  
  dsetrank = 1
  globdims = 0_8
  chnkdims = 0_8
  chnkoffs = 0_8
  memoffs  = 0_8
  memdims  = 0_8
  memsize  = 1_8
  
  globdims(1) = int(ncohorts_global,8)
  chnkdims(1) = int(cpatch%ncohorts,8)
  chnkoffs(1) = int(paco_index - 1,8)

  memdims(1)  = int(cpatch%ncohorts,8)
  memsize(1)  = int(cpatch%ncohorts,8)
  memoffs(1)  = 0_8

  if(cpatch%ncohorts>0) then

     call hdf_getslab_i(cpatch%pft,'PFT ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%nplant,'NPLANT ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%hite,'HITE ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%dbh,'DBH ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%agb,'AGB_CO ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%basarea,'BA_CO',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%dagb_dt,'DAGB_DT ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%dba_dt,'DBA_DT',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%ddbh_dt,'DDBH_DT',dsetrank,iparallel,.true.)

     call hdf_getslab_r(cpatch%bdead,'BDEAD ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%bleaf,'BLEAF ',dsetrank,iparallel,.true.)
     call hdf_getslab_i(cpatch%phenology_status,'PHENOLOGY_STATUS ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%balive,'BALIVE ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%broot,'BROOT  ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%bsapwood,'BSAPWOOD ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%lai,'LAI_CO ',dsetrank,iparallel,.true.)
     
     call hdf_getslab_r(cpatch%llspan,'LLSPAN ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%turnover_amp,'TURNOVER_AMP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%vm_bar,'VM_BAR ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%sla,'SLA ',dsetrank,iparallel,.true.)

     call hdf_getslab_r(cpatch%bstorage,'BSTORAGE ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%cbr_bar,'CBR_BAR ',dsetrank,iparallel,.true.)
     
     call hdf_getslab_r(cpatch%veg_temp,'VEG_TEMP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%veg_water,'VEG_WATER ',dsetrank,iparallel,.true.)
     
     call hdf_getslab_r(cpatch%wpa,'WPA_CO ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%wai,'WAI_CO ',dsetrank,iparallel,.true.)

     call hdf_getslab_r(cpatch%veg_energy,'VEG_ENERGY ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%hcapveg,'HCAPVEG ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%veg_fliq,'VEG_FLIQ ',dsetrank,iparallel,.true.)
     
     call hdf_getslab_r(cpatch%mean_gpp,'MEAN_GPP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%mean_leaf_resp,'MEAN_LEAF_RESP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%mean_root_resp,'MEAN_ROOT_RESP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%today_leaf_resp,'TODAY_LEAF_RESP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%today_root_resp,'TODAY_ROOT_RESP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%today_gpp,'TODAY_GPP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%today_gpp_pot,'TODAY_GPP_POT ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%today_gpp_max,'TODAY_GPP_MAX ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%growth_respiration,'GROWTH_RESPIRATION ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%storage_respiration,'STORAGE_RESPIRATION ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%vleaf_respiration,'VLEAF_RESPIRATION ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%fsn,'FSN ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%monthly_dndt,'MONTHLY_DNDT ',dsetrank,iparallel,.true.)

     if (associated(cpatch%mmean_gpp       )) &
          call hdf_getslab_r(cpatch%mmean_gpp,'MMEAN_GPP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_leaf_resp       )) &
          call hdf_getslab_r(cpatch%mmean_leaf_resp,'MMEAN_LEAF_RESP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_root_resp       )) &
          call hdf_getslab_r(cpatch%mmean_root_resp,'MMEAN_ROOT_RESP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_growth_resp       )) &
          call hdf_getslab_r(cpatch%mmean_growth_resp,'MMEAN_GROWTH_RESP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_storage_resp       )) &
          call hdf_getslab_r(cpatch%mmean_storage_resp,'MMEAN_STORAGE_RESP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_vleaf_resp       )) &
          call hdf_getslab_r(cpatch%mmean_vleaf_resp,'MMEAN_VLEAF_RESP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_leaf_resp       )) &
          call hdf_getslab_r(cpatch%dmean_leaf_resp,'DMEAN_LEAF_RESP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_root_resp       )) &
          call hdf_getslab_r(cpatch%dmean_root_resp,'DMEAN_ROOT_RESP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_gpp       )) &
          call hdf_getslab_r(cpatch%dmean_gpp,'DMEAN_GPP_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_fs_open       )) &
     call hdf_getslab_r(cpatch%dmean_fs_open,'DMEAN_FS_OPEN_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_fs_open       )) &
     call hdf_getslab_r(cpatch%mmean_fs_open,'MMEAN_FS_OPEN_CO ',dsetrank,iparallel,.false.) 
     if (associated(cpatch%dmean_fsw       )) &
     call hdf_getslab_r(cpatch%dmean_fsw,'DMEAN_FSW_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_fsw       )) &
     call hdf_getslab_r(cpatch%mmean_fsw,'MMEAN_FSW_CO ',dsetrank,iparallel,.false.) 
     if (associated(cpatch%dmean_fsn       )) &
     call hdf_getslab_r(cpatch%dmean_fsn,'DMEAN_FSN_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_fsn       )) &
     call hdf_getslab_r(cpatch%mmean_fsn,'MMEAN_FSN_CO ',dsetrank,iparallel,.false.) 
     if (associated(cpatch%mmean_leaf_maintenance )) &
     call hdf_getslab_r(cpatch%mmean_leaf_maintenance,'MMEAN_LEAF_MAINTENANCE ',dsetrank,iparallel,.false.)   
     if (associated(cpatch%mmean_root_maintenance )) &
     call hdf_getslab_r(cpatch%mmean_root_maintenance,'MMEAN_ROOT_MAINTENANCE ',dsetrank,iparallel,.false.)   
     if (associated(cpatch%mmean_leaf_drop       )) &
     call hdf_getslab_r(cpatch%mmean_leaf_drop,'MMEAN_LEAF_DROP_CO ',dsetrank,iparallel,.false.)   
     if (associated(cpatch%mmean_cb       )) &
     call hdf_getslab_r(cpatch%mmean_cb,'MMEAN_CB ',dsetrank,iparallel,.false.)   
     if (associated(cpatch%dmean_light_level       )) &
     call hdf_getslab_r(cpatch%dmean_light_level,'DMEAN_LIGHT_LEVEL ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_light_level       )) &
     call hdf_getslab_r(cpatch%mmean_light_level,'MMEAN_LIGHT_LEVEL ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_light_level       )) &
     call hdf_getslab_r(cpatch%dmean_light_level_beam,'DMEAN_LIGHT_LEVEL_BEAM ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_light_level_beam       )) &
     call hdf_getslab_r(cpatch%mmean_light_level_beam,'MMEAN_LIGHT_LEVEL_BEAM ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_light_level_diff       )) &
     call hdf_getslab_r(cpatch%dmean_light_level_diff,'DMEAN_LIGHT_LEVEL_DIFF ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_light_level_diff       )) &
     call hdf_getslab_r(cpatch%mmean_light_level_diff,'MMEAN_LIGHT_LEVEL_DIFF ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_par_v       )) &
     call hdf_getslab_r(cpatch%dmean_par_v,'DMEAN_PAR_V ',dsetrank,iparallel,.false.)

     if (associated(cpatch%dmean_par_v_beam       )) &
     call hdf_getslab_r(cpatch%dmean_par_v_beam,'DMEAN_PAR_V_BEAM ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_par_v_diff       )) &
     call hdf_getslab_r(cpatch%dmean_par_v_diff,'DMEAN_PAR_V_DIFF ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_par_v       )) &
     call hdf_getslab_r(cpatch%mmean_par_v,'MMEAN_PAR_V ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_par_v_beam       )) &
     call hdf_getslab_r(cpatch%mmean_par_v_beam,'MMEAN_PAR_V_BEAM ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_par_v_diff       )) &
     call hdf_getslab_r(cpatch%mmean_par_v_diff,'MMEAN_PAR_V_DIFF ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_beamext_level       )) &
     call hdf_getslab_r(cpatch%dmean_beamext_level,'DMEAN_BEAMEXT_LEVEL ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_beamext_level      )) &
     call hdf_getslab_r(cpatch%mmean_beamext_level,'MMEAN_BEAMEXT_LEVEL ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_diffext_level       )) &
     call hdf_getslab_r(cpatch%dmean_diffext_level,'DMEAN_DIFFEXT_LEVEL ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_diffext_level       )) &
     call hdf_getslab_r(cpatch%mmean_diffext_level,'MMEAN_DIFFEXT_LEVEL ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_diffext_level       )) &
     call hdf_getslab_r(cpatch%dmean_norm_par_beam,'DMEAN_NORM_PAR_BEAM ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_norm_par_beam       )) &
     call hdf_getslab_r(cpatch%mmean_norm_par_beam,'MMEAN_NORM_PAR_BEAM ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_norm_par_diff       )) &
     call hdf_getslab_r(cpatch%dmean_norm_par_diff,'DMEAN_NORM_PAR_DIFF ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_norm_par_diff       )) &
     call hdf_getslab_r(cpatch%mmean_norm_par_diff,'MMEAN_NORM_PAR_DIFF ',dsetrank,iparallel,.false.)
     if (associated(cpatch%dmean_lambda_light       )) &
          call hdf_getslab_r(cpatch%dmean_lambda_light,'DMEAN_LAMBDA_LIGHT_CO ',dsetrank,iparallel,.false.)
     if (associated(cpatch%mmean_lambda_light       )) &
     call hdf_getslab_r(cpatch%mmean_lambda_light,'MMEAN_LAMBDA_LIGHT_CO ',dsetrank,iparallel,.false.)

     call hdf_getslab_r(cpatch%Psi_open,'PSI_OPEN ',dsetrank,iparallel,.true.)
     call hdf_getslab_i(cpatch%krdepth,'KRDEPTH ',dsetrank,iparallel,.true.)
     call hdf_getslab_i(cpatch%first_census,'FIRST_CENSUS ',dsetrank,iparallel,.true.)
     call hdf_getslab_i(cpatch%new_recruit_flag,'NEW_RECRUIT_FLAG ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%light_level,'LIGHT_LEVEL ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%light_level_beam,'LIGHT_LEVEL_BEAM ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%light_level_diff,'LIGHT_LEVEL_DIFF ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%beamext_level,'BEAMEXT_LEVEL ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%diffext_level,'DIFFEXT_LEVEL ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%norm_par_beam,'NORM_PAR_BEAM ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%norm_par_diff,'NORM_PAR_DIFF ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%lambda_light,'LAMBDA_LIGHT_CO ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%par_v,'PAR_V ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%par_v_beam,'PAR_V_BEAM ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%par_v_diffuse,'PAR_V_DIFFUSE ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rshort_v,'RSHORT_V ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rshort_v_beam,'RSHORT_V_BEAM ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rshort_v_diffuse,'RSHORT_V_DIFFUSE ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rlong_v,'RLONG_V ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rlong_v_surf,'RLONG_V_SURF ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rlong_v_incid,'RLONG_V_INCID ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rb,'RB ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%A_open,'A_OPEN ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%A_closed,'A_CLOSED ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%Psi_closed,'PSI_CLOSED ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rsw_open,'RSW_OPEN ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%rsw_closed,'RSW_CLOSED ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%fsw,'FSW ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%fs_open,'FS_OPEN ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%stomatal_resistance,'STOMATAL_RESISTANCE ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%leaf_maintenance,'LEAF_MAINTENANCE ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%root_maintenance,'ROOT_MAINTENANCE ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%leaf_drop,'LEAF_DROP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%bseeds,'BSEEDS_CO ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%leaf_respiration,'LEAF_RESPIRATION ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%root_respiration,'ROOT_RESPIRATION ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%gpp,'GPP ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%paw_avg,'PAW_AVG ',dsetrank,iparallel,.true.)
     
     !----- 13-month dimension (12 previous months + current month). ----------------------!
     dsetrank    = 2
     globdims(1) = 13_8
     chnkdims(1) = 13_8
     chnkoffs(1) = 0_8
     memdims(1)  = 13_8
     memsize(1)  = 13_8
     memoffs(2)  = 0_8
     
     globdims(2) = int(ncohorts_global,8)
     chnkdims(2) = int(cpatch%ncohorts,8)
     chnkoffs(2) = int(paco_index - 1,8)
     
     memdims(2)  = int(cpatch%ncohorts,8)
     memsize(2)  = int(cpatch%ncohorts,8)
     memoffs(2)  = 0_8
     
     call hdf_getslab_r(cpatch%cb,'CB ',dsetrank,iparallel,.true.)
     call hdf_getslab_r(cpatch%cb_max,'CB_MAX ',dsetrank,iparallel,.true.)

     !----- 2-D, dimensioned by the number of mortality rates. ----------------------------!
     dsetrank    = 2
     globdims(1) = int(n_mort,8)
     chnkdims(1) = int(n_mort,8)
     chnkoffs(1) = 0_8
     memdims(1)  = int(n_mort,8)
     memsize(1)  = int(n_mort,8)
     memoffs(2)  = 0_8
     
     globdims(2) = int(ncohorts_global,8)
     chnkdims(2) = int(cpatch%ncohorts,8)
     chnkoffs(2) = int(paco_index - 1,8)
     
     memdims(2)  = int(cpatch%ncohorts,8)
     memsize(2)  = int(cpatch%ncohorts,8)
     memoffs(2)  = 0_8

     call hdf_getslab_r(cpatch%mort_rate,'MORT_RATE_CO ',dsetrank,iparallel,.true.)


     dsetrank    = 2
     globdims(1) = int(n_stoma_atts,8)
     chnkdims(1) = int(n_stoma_atts,8)
     chnkoffs(1) = 0_8
     memdims(1)  = int(n_stoma_atts,8)
     memsize(1)  = int(n_stoma_atts,8)
     memoffs(2)  = 0_8
     
     globdims(2) = int(ncohorts_global,8)
     chnkdims(2) = int(cpatch%ncohorts,8)
     chnkoffs(2) = int(paco_index - 1,8)
     
     memdims(2)  = int(cpatch%ncohorts,8)
     memsize(2)  = int(cpatch%ncohorts,8)
     memoffs(2)  = 0_8
     

     call hdf_getslab_r(cpatch%old_stoma_vector,'OLD_STOMA_VECTOR', &
          dsetrank,iparallel,.true.)

  cohortloop: do ico=1,cpatch%ncohorts
     cpatch%old_stoma_data(ico)%recalc = int(cpatch%old_stoma_vector(1,ico))
     cpatch%old_stoma_data(ico)%T_L    = cpatch%old_stoma_vector(2,ico)
     cpatch%old_stoma_data(ico)%e_A    = cpatch%old_stoma_vector(3,ico)
     cpatch%old_stoma_data(ico)%PAR    = cpatch%old_stoma_vector(4,ico)
     cpatch%old_stoma_data(ico)%rb_factor = cpatch%old_stoma_vector(5,ico)
     cpatch%old_stoma_data(ico)%prss = cpatch%old_stoma_vector(6,ico) 
     cpatch%old_stoma_data(ico)%phenology_factor = cpatch%old_stoma_vector(7,ico)
     cpatch%old_stoma_data(ico)%gsw_open = cpatch%old_stoma_vector(8,ico)
     cpatch%old_stoma_data(ico)%ilimit   = int(cpatch%old_stoma_vector(9,ico))
     cpatch%old_stoma_data(ico)%T_L_residual = cpatch%old_stoma_vector(10,ico)
     cpatch%old_stoma_data(ico)%e_a_residual = cpatch%old_stoma_vector(11,ico)
     cpatch%old_stoma_data(ico)%par_residual = cpatch%old_stoma_vector(12,ico)
     cpatch%old_stoma_data(ico)%rb_residual  = cpatch%old_stoma_vector(13,ico)
     cpatch%old_stoma_data(ico)%prss_residual= cpatch%old_stoma_vector(14,ico) 
     cpatch%old_stoma_data(ico)%leaf_residual= cpatch%old_stoma_vector(15,ico)
     cpatch%old_stoma_data(ico)%gsw_residual = cpatch%old_stoma_vector(16,ico)
  enddo cohortloop

endif


  return
end subroutine fill_history_patch
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine hdf_getslab_r(buff,varn,dsetrank,iparallel,required)
  
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
  
  ! Some datasets are required during model initialization, such as the state
  ! variables like soil, canopy and vegetation energy and mass, although some variables
  ! are involved only in diagnostics.  If they are missing from the restart dataset,
  ! that may compromise diagnostic variables that are being averaged over the current time
  ! period, but they will not compromise the transition of the prognostic model state from one
  ! simulation to the next.
  ! If the dataset is not required, pass it in as a .true. argument

  logical,intent(in) :: required

  ! If the the optional argument is not present, take the conservative stance
  ! and make sure that it is "not-not-not-not-required", "not-not-required", or "required".

  call h5dopen_f(file_id,trim(varn), dset_id, hdferr)
  if (hdferr /= 0 .and. required ) then

     write(unit=*,fmt=*) 'File_ID = ',file_id
     write(unit=*,fmt=*) 'Dset_ID = ',dset_id
     call fatal_error('Could not get the dataset for '//trim(varn)//'!!!' &
          ,'hdf_getslab_r','ed_history_io.f90')
     
  else if (hdferr /= 0 .and. .not.required ) then

     write(unit=*,fmt=*) 'File_ID = ',file_id
     write(unit=*,fmt=*) 'Dset_ID = ',dset_id
     write (unit=*,fmt='(a)') '-----------------------------------------------------------'
     write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
     write (unit=*,fmt='(a)') '                                                           '
     write (unit=*,fmt='(a)') ' + Variable '//trim(varn)//' not found in your history.'
     write (unit=*,fmt='(a)') ' + Initializing this variable with zero. '
     write (unit=*,fmt='(a)') ' + his may cause some of your diagnostic output related'
     write (unit=*,fmt='(a)') '   to this variable to be incorrect the current period.'
     write (unit=*,fmt='(a)') ''
     write (unit=*,fmt='(a)') '   This variable has been specified as:'
     write (unit=*,fmt='(a)') '   NOT ABSOUTELY NECESSARY TO RESTART THE PROGNOSTIC STATE'
     write (unit=*,fmt='(a)') '-----------------------------------------------------------'
     write (unit=*,fmt='(a)') ''
     
     buff=0.
     return
     
  else
     
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
           call fatal_error('Could not read in the hyperslab dataset for '//trim(varn)//'!!!' &
                ,'hdf_getslab_r','ed_history_io.f90')
        end if

     else

        call h5dread_f(dset_id, H5T_NATIVE_REAL,buff,globdims, hdferr, &
             mem_space_id = memspace, file_space_id = filespace )

        if (hdferr /= 0) then
           call fatal_error('Could not read in the hyperslab dataset for '//trim(varn)//'!!!' &
                ,'hdf_getslab_r','ed_history_io.f90')
        end if

     endif
     
     !  write(unit=*,fmt='(a)') 'History start: Loading '//trim(varn)//'...'
     
     call h5sclose_f(filespace, hdferr)
     call h5sclose_f(memspace , hdferr)
     call h5dclose_f(dset_id  , hdferr)
     
  endif
  
  
  return
end subroutine hdf_getslab_r

!==========================================================================================!
!==========================================================================================!

subroutine hdf_getslab_d(buff,varn,dsetrank,iparallel,required)
  
  use hdf5
  use hdf5_coms,only:file_id,dset_id,dspace_id,plist_id, &
       filespace,memspace, &
       globdims,chnkdims,chnkoffs,cnt,stride, &
       memdims,memoffs,memsize
  
  
  implicit none
  
  real(kind=8),dimension(memsize(1),memsize(2),memsize(3),memsize(4)) :: buff

  integer :: hdferr,dsetrank
  integer :: iparallel
  character(len=*),intent(in) :: varn

  ! Some datasets are required during model initialization, such as the state
  ! variables like soil, canopy and vegetation energy and mass, although some variables
  ! are involved only in diagnostics.  If they are missing from the restart dataset,
  ! that may compromise diagnostic variables that are being averaged over the current time
  ! period, but they will not compromise the transition of the prognostic model state from one
  ! simulation to the next.
  ! If the dataset is not required, pass it in as a .true. argument

  logical,intent(in)  :: required


  call h5dopen_f(file_id,trim(varn), dset_id, hdferr)
  if (hdferr /= 0 .and. required ) then

     write(unit=*,fmt=*) 'File_ID = ',file_id
     write(unit=*,fmt=*) 'Dset_ID = ',dset_id
     call fatal_error('Could not get the dataset for '//trim(varn)//'!!!' &
          ,'hdf_getslab_d','ed_history_io.f90')
     
  else if (hdferr /= 0 .and. .not.required ) then

     write(unit=*,fmt=*) 'File_ID = ',file_id
     write(unit=*,fmt=*) 'Dset_ID = ',dset_id
     write (unit=*,fmt='(a)') '-----------------------------------------------------------'
     write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
     write (unit=*,fmt='(a)') '                                                           '
     write (unit=*,fmt='(a)') ' + Variable '//trim(varn)//' not found in your history.'
     write (unit=*,fmt='(a)') ' + Initializing this variable with zero. '
     write (unit=*,fmt='(a)') ' + This may cause some of your diagnostic output related'
     write (unit=*,fmt='(a)') '   to this variable to be incorrect the current period.'
     write (unit=*,fmt='(a)') ''
     write (unit=*,fmt='(a)') '   This variable has been specified as:'
     write (unit=*,fmt='(a)') '   NOT ABSOLUTELY NECESSARY TO RESTART THE PROGNOSTIC STATE'
     write (unit=*,fmt='(a)') '-----------------------------------------------------------'
     write (unit=*,fmt='(a)') ''
     
     buff=0.d0
     return
     
  else
  
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
        
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,buff,globdims, hdferr, &
             mem_space_id = memspace, file_space_id = filespace, &
             xfer_prp = plist_id)
        if (hdferr /= 0) then
           call fatal_error('Could not read in the hyperslab dataset for '//trim(varn)//'!!!' &
                ,'hdf_getslab_r','ed_history_io.f90')
        end if
        
     else
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,buff,globdims, hdferr, &
             mem_space_id = memspace, file_space_id = filespace )
        
        if (hdferr /= 0) then
           call fatal_error('Could not read in the hyperslab dataset for '//trim(varn)//'!!!' &
                ,'hdf_getslab_r','ed_history_io.f90')
        end if
     endif
     
     !  write(unit=*,fmt='(a)') 'History start: Loading '//trim(varn)//'...'
     
     call h5sclose_f(filespace, hdferr)
     call h5sclose_f(memspace , hdferr)
     call h5dclose_f(dset_id  , hdferr)

  endif
  

  return
end subroutine hdf_getslab_d
!==========================================================================================!
!==========================================================================================!




!==========================================================================================!
!==========================================================================================!
subroutine hdf_getslab_i(buff,varn,dsetrank,iparallel,required)
  
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
  ! Some datasets are required during model initialization, such as the state
  ! variables like soil, canopy and vegetation energy and mass, although some variables
  ! are involved only in diagnostics.  If they are missing from the restart dataset,
  ! that may compromise diagnostic variables that are being averaged over the current time
  ! period, but they will not compromise the transition of the prognostic model state from one
  ! simulation to the next.
  ! If the dataset is not required, pass it in as a .true. argument

  logical,intent(in) :: required
  
  call h5dopen_f(file_id,trim(varn), dset_id, hdferr)
  if (hdferr /= 0 .and. required ) then
     
     write(unit=*,fmt=*) 'File_ID = ',file_id
     write(unit=*,fmt=*) 'Dset_ID = ',dset_id
     call fatal_error('Could not get the dataset for '//trim(varn)//'!!!' &
          ,'hdf_getslab_i','ed_history_io.f90')
     
  else if (hdferr /= 0 .and. .not.required ) then
     
     write(unit=*,fmt=*) 'File_ID = ',file_id
     write(unit=*,fmt=*) 'Dset_ID = ',dset_id
     write (unit=*,fmt='(a)') '-----------------------------------------------------------'
     write (unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   '
     write (unit=*,fmt='(a)') '                                                           '
     write (unit=*,fmt='(a)') ' + Variable '//trim(varn)//' not found in your history.'
     write (unit=*,fmt='(a)') ' + Initializing this variable with zero. '
     write (unit=*,fmt='(a)') ' + This may cause some of your diagnostic output related'
     write (unit=*,fmt='(a)') '   to this variable to be incorrect the current period.'
     write (unit=*,fmt='(a)') ''
     write (unit=*,fmt='(a)') '   This variable has been specified as:'
     write (unit=*,fmt='(a)') '   NOT ABSOLUTELY NECESSARY TO RESTART THE PROGNOSTIC STATE'
     write (unit=*,fmt='(a)') '-----------------------------------------------------------'
     write (unit=*,fmt='(a)') ''
     
     buff=0
     return
     
  else
  
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

  endif
     
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
     flat = ed_res * real(int(lat / ed_res)) + 0.5 * ed_res 
  else
     flat = - ed_res * real(int(-lat / ed_res)) - 0.5 * ed_res
  endif
  
  if(lon >= 0.0)then
     flon = ed_res * real(int(lon / ed_res)) + 0.5 * ed_res 
  else
     flon = - ed_res * real(int(-lon / ed_res)) - 0.5 * ed_res
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
  elseif(ed_res > 0.0999 .and. ed_res < 0.1001)then

     if(flat <= -10.0)then
        write(ed_fname,'(a,f6.2)')trim(sfilin)//'lat',flat
     elseif(flat < 0.0 .or. flat >= 10.0)then
        write(ed_fname,'(a,f5.2)')trim(sfilin)//'lat',flat
     else
        write(ed_fname,'(a,f4.2)')trim(sfilin)//'lat',flat
     endif
     if(flon <= -100.0)then
        write(ed_fname,'(a,f7.2)')trim(ed_fname)//'lon',flon
     elseif(flon <= -10.0 .or. flon >= 100.0)then
        write(ed_fname,'(a,f6.2)')trim(ed_fname)//'lon',flon
     elseif(flon < 0.0)then
        write(ed_fname,'(a,f5.2)')trim(ed_fname)//'lon',flon
     elseif(flon < 10.0)then
        write(ed_fname,'(a,f4.2)')trim(ed_fname)//'lon',flon
     else
        write(ed_fname,'(a,f5.2)')trim(ed_fname)//'lon',flon
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
