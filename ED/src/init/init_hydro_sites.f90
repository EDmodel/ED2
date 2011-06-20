!!! translate from C to fortran
! ===============================================================

subroutine read_site_file(cgrid,igr)
   ! function for loading sites within polygons and initializing polygon parms
   ! call prior to loading pss/css files but after basic polygon established
  
   use soil_coms, only: soil,slz
   use grid_coms, only: nzg
   use ed_misc_coms, only: ied_init_mode,sfilin, vary_elev, vary_rad, vary_hyd
   use mem_polygons, only: edres
   use ed_state_vars, only: edtype,polygontype,sitetype,allocate_polygontype
   use ed_max_dims, only: max_site,n_pft,str_len

   implicit none
   integer :: igr
   logical :: no_rad   = .false.  !! true turns effect OFF
   logical :: no_lapse = .false.
   logical :: no_hyd   = .false.

   character(len=str_len) :: site_name,pss_name,css_name
  
   type(edtype) :: cgrid
   type(polygontype),pointer :: cpoly

   character(len=200) :: cdummy,cdummy2
   integer :: i,nsc,nsites=1
   real(kind=8) :: area_sum = 0.0d+0
   integer :: sc             ! soil classa
   integer :: fformat = 0    ! file format
   logical :: fexist         ! file exists
   real :: Te,T0,K0
   real :: fa,fb,zmin=-2.0
   integer :: sitenum,get_site_line,get_mat_line,get_header,found_mat_header=0,lcount=0,mcount=0
   real :: area,TCI,elevation,slope,aspect
   integer,allocatable :: soilclass(:)
   integer :: ipy,isi,k
   integer :: ierr

   if(vary_elev == 0)  no_lapse = .true.
   if(vary_rad == 0)   no_rad = .true.
   if(vary_hyd == 0) no_hyd = .true.




   ! ASSUMING FOR NOW THAT THERE IS NO WATER SITE

   ! init values
   if (associated(cgrid%load_adjacency)) cgrid%load_adjacency = 0

  
   do ipy = 1,cgrid%npolygons

      cpoly => cgrid%polygon(ipy)

      call create_ed10_ed20_fname(cgrid%lat(ipy), edres, cgrid%lon(ipy) &
                                 ,trim(sfilin(igr)),pss_name,css_name,site_name)
      !! check if site file exists
      inquire(file=trim(site_name),exist=fexist)
      if(.not.fexist) then
         print*,"error opening site file ",site_name
         print*,"setting ied_init_mode to 2 and loading as single site"
      end if

      
      ! If there is no terrestrial site
      if(.not. fexist) then  !! have no site file or ied_init_mode is not 3

         ! Allocate single site vector information to
         ! the swap polygon
         
         call allocate_polygontype(cpoly,1)

         cpoly%lsl(1)  = cgrid%lsl(ipy)  ! Initialize lowest soil layer
         cpoly%area(1) = 1.0             ! Initialize the area to all

         ! Copy the soils information from the polygon to the site
         do k=1,nzg
            cpoly%ntext_soil(k,1) = cgrid%ntext_soil(k,ipy)
         enddo

         ! Set soil moisture decay function, based on second layer's K value
         ! use the second layer instead of the top in case top is organic/peat
         sc = cpoly%ntext_soil(nzg-1,1)
         cpoly%moist_f(1) = -log(soil(sc)%slcons / soil(sc)%slcons0) / 2.0

         !! derive adjustments to f
         zmin = slz(cpoly%lsl(1))
         fa = -1.0/zmin !! should be 1/(depth to bedrock)
         if(cpoly%moist_f(1)*zmin < 0.0) then
            fb = cpoly%moist_f(1)/(1.0-exp(cpoly%moist_f(1)*zmin))
            cpoly%moist_f(1) = max(fa,fb)
         else
            cpoly%moist_f(1) = fa
         endif

         cpoly%sitenum(1)   = 1
         cpoly%elevation(1) = 0.0
         cpoly%slope(1)     = 0.0
         cpoly%aspect(1)    = 0.0
         cpoly%TCI(1)       = 0.0

      else

         !! Read data from site file
         open(unit=12,file=trim(site_name),form='formatted',status='old')

         !/* read file format line */
         read(unit=12,fmt=*)cdummy,nsites,cdummy2,fformat
         !/*line format: "nsites <nsites> format <fformat>" */   
         
!         print*,"reading",nsites,"sites using file format",fformat

         call allocate_polygontype(cpoly,nsites)
     
         if(fformat <=0 .or. fformat > 3) then
            print*,""
            print*,"ERROR :: unrecognized file format specifier in",site_name
            stop
         endif
         if(fformat == 2) cgrid%load_adjacency(ipy) = 1
         nsc = 1
         if(fformat == 3) nsc = nzg

         !/* discard file header line */
         read(unit=12,fmt=*)
         
         !/* read data*/
         allocate(soilclass(nzg))
         lcount = 0 
         mcount = 0
         found_mat_header = 0
         isi = 0
         area_sum = 0.0d+0
         count_sites: do

            isi = isi + 1  ! The site counter
            
            !           read(12,*,iostat=ierr)time,pname,trk,age,area,fsc,stsc,  &
            !                stsl,ssc,psc,msn,fsn,water(1:nwater)
            !           if(ierr /= 0)exit count_patches
            
            lcount = lcount + 1 !! line counter
            get_site_line=0
            get_mat_line=0
            get_header=0 !! line flags
            
            !/********* decide what type of line to read ***************/
            select case (fformat)
            case (1,3)
               if(lcount <= nsites) then !we're reading data
                  get_site_line=1
                  get_mat_line=0
                  get_header=0
               else                      !we're past data, discard
                  get_site_line=0
                  get_mat_line=0
                  get_header=1             
               endif
            case (2)
               if(lcount <= nsites) then  !know we're in the site section
                  get_site_line=1
                  get_mat_line=0

               elseif(found_mat_header.eq.1)then !know we're in the matrix section
                  if(mcount < nsites)then
                     get_site_line = 0
                     get_mat_line=1
                     get_header = 0
                  else  !//already got all the matrix lines, discard what's left
                     get_site_line = 0
                     get_mat_line = 0
                     get_header = 1
                  endif
               else 
                  
                  !!assume line is header, remove and then assume header found
                  
                  get_site_line = 0
                  get_mat_line = 0
                  get_header = 1
                  found_mat_header=1
                  
               endif
               
            end select

!            print*,"line indicators",get_site_line,get_mat_line,get_header
    
            if(get_site_line == 1)then   !/********** READ SITE LINE ***************/
               
               !/* line format: sitenum, area, TCI, elevation, slope, aspect,ntext_soil(s) */              
               read(unit=12,fmt=*,iostat=ierr)sitenum,area,TCI,elevation,slope,aspect,soilclass(1:nsc)
               if(ierr == 0) then
                  !/*create data object for each new site */
!                  print*,sitenum, area, TCI, elevation,slope,aspect,soilclass(1:nsc)
                  
                  cpoly%lsl(isi)  = cgrid%lsl(ipy)  ! Initialize lowest soil layer
                  cpoly%area(isi) = area            ! Initialize the area to all
                  
                  ! Copy the soils information from the polygon to the site
                  do k=1,nzg
                     cpoly%ntext_soil(k,isi) = cgrid%ntext_soil(k,ipy)
                  enddo

                  area_sum = area_sum + dble(area)
                  cpoly%sitenum(isi)      = sitenum
                  cpoly%elevation(isi)    = elevation
                  cpoly%slope(isi)        = slope
                  cpoly%aspect(isi)       = aspect
                  cpoly%TCI(isi)          = TCI+13.96962

                  !! flags to turn effects off
                  if(no_lapse) cpoly%elevation(isi) = 0.0
                  if(no_hyd)   cpoly%TCI(isi)       = 8.0
                  if(no_rad) then
                               cpoly%slope(isi)     = 0.0
                               cpoly%aspect(isi)    = 0.0
                  end if

!print*,"SITE",cpoly%elevation(isi),cpoly%TCI(isi),cpoly%slope(isi),cpoly%aspect(isi)

                  if(fformat == 3) then
                     do i=1,nzg
                        cpoly%ntext_soil(i,isi) = soilclass(i)
                     end do
                  else
                     do i=1,nzg
                        cpoly%ntext_soil(i,isi) = soilclass(1)
                     end do
                  end if
                  !//Currently do nothing with setting site-level soils

                  sc = cpoly%ntext_soil(nzg-1,1)
                  cpoly%moist_f(isi) = -log(soil(sc)%slcons / soil(sc)%slcons0) / 2.0
                  !! derive adjustments to f
                  zmin = slz(cpoly%lsl(isi))
                  fa = -1.0/zmin !! should be 1/(depth to bedrock)
                  if(cpoly%moist_f(isi)*zmin < 0.0) then
                     fb = cpoly%moist_f(isi)/(1.0-exp(cpoly%moist_f(isi)*zmin))
                     cpoly%moist_f(isi) = max(fa,fb)
                  else
                     cpoly%moist_f(isi) = fa
                  endif
               end if  ! end valid line
            end if        !//********** END READ SITE LINE**************

            ! RGK
            if(get_header == 1) read(unit=12,fmt=*,iostat=ierr) !/* discard file line */
            ! RGK
            if(get_mat_line == 1) then

               !/*read line into site adjacency matrix*/
               
               read(unit=12,fmt=*,iostat=ierr) cgrid%site_adjacency(mcount,1:(nsites+1),ipy)

               mcount = mcount+1
               
            endif
            if(ierr /= 0) exit count_sites
        
         end do count_sites
         deallocate(soilclass)
      end if

      !adjust areas
      !assume that if area_sum ~ 1 need to renormalize terrestrial
      if( area_sum > 0.995d+0) then
         cpoly%area(:) = real(dble(cpoly%area(:))/area_sum)
      end if


      if(cgrid%load_adjacency(ipy) /= 0) then  
         call calc_flow_routing(cgrid,ipy)
      endif

      ! calculate summary stats - pass 1: Te
      ! On terrestrial site still
      Te = 0.0 
      area_sum = 0.0d+0
      do isi = 1,cpoly%nsites
         sc = cpoly%ntext_soil(nzg-1,isi)
         K0 = soil(sc)%slcons0
         T0 = K0/cpoly%moist_f(isi)
         Te = Te + T0*cpoly%area(isi)
         area_sum = area_sum + dble(cpoly%area(isi))
      end do

      Te = Te/real(area_sum)
      cgrid%Te(ipy) = Te

      !pass 2: W, Wbar

      cgrid%wbar(ipy) = 0.0

      do isi = 1,cpoly%nsites
         sc = cpoly%ntext_soil(nzg-1,isi)
         K0 = soil(sc)%slcons0
         T0 = K0/cpoly%moist_f(isi)
         cpoly%moist_W(isi) = cpoly%TCI(isi) + log(Te) - log(T0)
         cgrid%wbar(ipy) = cgrid%wbar(ipy) + real(dble(cpoly%moist_W(isi))*dble(cpoly%area(isi))/area_sum)
      end do

      !call dump_ed(myPolygon)

   end do


end subroutine read_site_file
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine calc_flow_routing(cgrid,ipy)
   use ed_state_vars, only: polygontype, edtype
   type(edtype), target :: cgrid ! Alias for current grid
   integer, intent(in) :: ipy    ! Current polygon Polygon ID
   integer :: i,ihys,ines ! ihys -> current hydro site ID. ines-> next hydro site ID
   type(polygontype),pointer :: cpoly ! Alias for current polygon
   real :: side = 10.0  !cell side length (m) used for calculating adjacency
   real(kind=8) :: row_sum
   integer, external :: find_rank
   integer, allocatable, dimension(:) :: hyrank
   integer :: nsites

   !if we loaded an adjacency matrix from file, we don't need to calculate flow routing
   if(cgrid%load_adjacency(ipy) == 1) return   

   cpoly => cgrid%polygon(ipy)
   nsites=cpoly%nsites
   
   allocate (hyrank(nsites))

   ! set flow routing -> sort by TCI, it will give the sites the hydro order
   call rank_up(nsites,cpoly%TCI,hyrank)


   !init structures and recalculate sitenums list based on hydro order
   do ihys=1,nsites
      do ines=1,nsites
         cgrid%site_adjacency(ihys,ines,ipy) = 0.0
      end do
   end do
   !routing established, now calc approximate adjacency matrix
   do i=1,nsites-1
      ihys=find_rank(i,nsites,hyrank)
      ines=find_rank(i+1,nsites,hyrank)
      cgrid%site_adjacency(ihys,ihys,ipy) = 2.0*side*(side-1.0)*cpoly%area(ihys) &
               + side/2.0*cpoly%area(ihys)/(cpoly%area(ihys)+cpoly%area(ines))
      cgrid%site_adjacency(ihys,ines,ipy) = side*cpoly%area(ines)/(cpoly%area(ihys)+cpoly%area(ihys))
   end do
   !Now find the adjacency for the last rank site:
   ihys=find_rank(nsites,nsites,hyrank)
   ines=find_rank(1,nsites,hyrank)
   cgrid%site_adjacency(ihys,ihys,ipy) = 2.0*side*(side-1.0)*cpoly%area(ihys) &
            + side/2.0*cpoly%area(ihys)/(cpoly%area(ihys)+cpoly%area(ines))
   cgrid%site_adjacency(ihys,ines,ipy) = side*cpoly%area(ines)/(cpoly%area(ihys)+cpoly%area(ihys))

   !normalize routing (rows sum to 1)  
   do ihys=1,nsites
      row_sum=sum(dble(cgrid%site_adjacency(ihys,:,ipy)))
      cgrid%site_adjacency(ihys,:,ipy)=real(dble(cgrid%site_adjacency(ihys,:,ipy))/row_sum)
   end do

   !check routing
!   print*,"CHECK ROUTING" 
!   do i=1,myPolygon%nsites
!      print*,cgrid%site_adjacency(i,:,ipy)
!   end do

   deallocate(hyrank)
   return
end subroutine calc_flow_routing
