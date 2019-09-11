module landuse_init_module
  contains

!==========================================================================================!
!==========================================================================================!
!     This file will read the landuse files and assign the anthropogenic disturbance       !
! matrices.                                                                                !
!------------------------------------------------------------------------------------------!
subroutine landuse_init

   use ed_state_vars , only : edtype          & ! structure
                            , polygontype     & ! structure
                            , sitetype        & ! structure
                            , edgrid_g        ! ! structure
   use consts_coms   , only : erad            & ! intent(in)
                            , pio180          ! ! intent(in)
   use disturb_coms  , only : lutime          & ! intent(in)
                            , max_lu_years    & ! intent(in)
                            , num_lu_trans    & ! intent(in)
                            , ianth_disturb   & ! intent(in)
                            , lu_database     ! ! intent(in)
!                             , sl_pft          & ! intent(in)
!                             , sl_scale        & ! intent(in)
!                             , sl_nyrs         & ! intent(in)
!                             , sl_yr_first     & ! intent(in)
!                             , sl_mindbh_harvest & ! intent(in)
!                             , sl_prob_harvest ! ! intent(in)
   use ed_misc_coms  , only : iyeara          & ! intent(in)
                            , iyearz          ! ! intent(in)
   use grid_coms     , only : ngrids          ! ! intent(in)
   use ed_max_dims   , only : str_len         & ! intent(in)
                            , huge_lu         & ! intent(in)
                            , n_pft           & ! intent(in)
                            , maxlist         ! ! intent(in)
  ! use detailed_coms , only : idetailed         ! ! intent(in) ! Maybe not needed
                            

   implicit none
   !----- Local variables -----------------------------------------------------------------!
   type(edtype)          , pointer            :: cgrid
   type(polygontype)     , pointer            :: cpoly
   type(sitetype)        , pointer            :: csite
   type(lutime)          , pointer            :: clutime
   type(lutime)          , pointer            :: onelutime
   character(len=str_len), dimension(maxlist) :: full_list
   character(len=str_len), dimension(maxlist) :: lu_list
   real                  , dimension(maxlist) :: llon_list
   real                  , dimension(maxlist) :: llat_list
   real                  , dimension(maxlist) :: file_ldist
   character(len=6)                           :: hform
   character(len=str_len)                     :: lu_name
   character(len=str_len)                     :: cdum
!    character(len=13)                          :: hifmt ! if writing LU table
!    character(len=15)                          :: hffmt ! if writing LU table  
   integer                                    :: nharvest
   integer               , dimension(n_pft)   :: harvest_pft
   integer                                    :: nf
   integer                                    :: nflist
   integer                                    :: nfllu
   integer                                    :: ncl
   integer                                    :: iyear
   integer                                    :: igr
   integer                                    :: ipy
   integer                                    :: isi
   integer                                    :: ipft
   integer                                    :: h
   integer                                    :: sim_years
   integer                                    :: yd_1st
   integer                                    :: yd_this
   integer                                    :: yd_last
   logical                                    :: inside
   real                  , dimension(n_pft)   :: mindbh_slog
   real                  , dimension(n_pft)   :: harvprob_slog_g
   real                  , dimension(n_pft)   :: harvprob_slog_l
   real                  , dimension(n_pft)   :: mindbh_fplt
   real                  , dimension(n_pft)   :: harvprob_fplt_g
   real                  , dimension(n_pft)   :: harvprob_fplt_l
   real                  , dimension(num_lu_trans) :: landuse_now
   real                                       :: lu_area
   real                                       :: lu_area_i
   real                                       :: wlon
   real                                       :: elon
   real                                       :: slat
   real                                       :: nlat
   !----- Local constants. ----------------------------------------------------------------!
   character(len=12)    , parameter           :: fffmt    = '(a,1x,f12.5)'
   character(len=13)    , parameter           :: esfmt    = '(a,1x,es12.5)'
   integer              , parameter           :: hoff     = 18
   real                 , parameter           :: huge_dbh = huge(1.)
   !----- External function. --------------------------------------------------------------!
   real                 , external            :: dist_gc
!    real                 , external            :: solid_area ! Need only for writing LU table
   !---------------------------------------------------------------------------------------!


   !----- Finding number of simulation years ----------------------------------------------!
   sim_years = iyearz-iyeara+1

   !------ Decide whether to write lu settings. ----------------------------------------!
   ! write_lu_settings = btest(idetailed,6) .and. ianth_disturb /= 0
   !------------------------------------------------------------------------------------!

   !----- Crashing the run if the user set up a very long run... --------------------------!
   if (ianth_disturb /= 0 .and. sim_years > max_lu_years) then
      write (unit=*,fmt='(a,1x,i5)') 'IYEARA       (From namelist)        :',iyeara
      write (unit=*,fmt='(a,1x,i5)') 'IYEARZ       (From namelist)        : ',iyearz
      write (unit=*,fmt='(a,1x,i5)') 'MAX_LU_YEARS (From disturb_coms.f90): ',max_lu_years
      write (unit=*,fmt='(a)') ' Your run is too long.  Try increasing max_lu_years,'
      write (unit=*,fmt='(a)') ' so MAX_LU_YEARS >= IYEARZ-IYEARA+1, then recompile'
      write (unit=*,fmt='(a)') ' your model.'
      call fatal_error ('Simulation is too long for anthropogenic disturbance.'            &
                       ,'landuse_init','landuse_init.f90')
   end if


   gridloop: do igr = 1,ngrids

      !------------------------------------------------------------------------------------!
      !     Find the list of disturbance rate files.                                       !
      !------------------------------------------------------------------------------------!
      if (ianth_disturb >= 1) then
         call ed_filelist(full_list,lu_database(igr),nflist)
         call ed1_fileinfo('.lu',nflist,full_list,nfllu,lu_list,llon_list,llat_list)
      end if
      !------------------------------------------------------------------------------------!

      cgrid=>edgrid_g(igr)

      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         select case (ianth_disturb)
         case (0)
            !------------------------------------------------------------------------------!
            !      Anthropogenic disturbance is not used this time, allocate only a single !
            ! landuse year.                                                                !
            !------------------------------------------------------------------------------!
            allocate(cpoly%clutimes(1,cpoly%nsites))
            !------------------------------------------------------------------------------!

            !----- Set the parameters in a way that no logging/ploughing will happen. -----!
            do isi = 1,cpoly%nsites
               cpoly%num_landuse_years(isi)                  = 1
               cpoly%mindbh_harvest   (1:n_pft,isi)          = huge_dbh
               cpoly%prob_havest_g    (1:n_pft,isi)          = 0.
               cpoly%prob_havest_l    (1:n_pft,isi)          = 0.
               cpoly%clutimes(1,isi)%landuse_year            = iyeara
               cpoly%clutimes(1,isi)%landuse(1:num_lu_trans) = 0.0
            end do
            !------------------------------------------------------------------------------!


            !----- No plantations. --------------------------------------------------------!
            cpoly%plantation(:) = 0
            !------------------------------------------------------------------------------!
            
         case (1,2)

            !------------------------------------------------------------------------------!
            !     Comput the distance between the current polygon and all the files.       !
            !------------------------------------------------------------------------------!
            do nf=1,nfllu
               file_ldist(nf) = dist_gc(cgrid%lon(ipy),llon_list(nf)                       &
                                       ,cgrid%lat(ipy),llat_list(nf) )
            end do

            !------------------------------------------------------------------------------!
            !    Pick the closest file.  This is not a guarantee that it will be used      !
            ! because the closest polygon area must contain the point associated to the    !
            ! current polygon.                                                             !
            !------------------------------------------------------------------------------!
            ncl     = minloc(file_ldist(1:nfllu),dim=1)
            lu_name = lu_list(ncl)
            write (unit=*,fmt='(2a)') 'Using land use file: ',trim(lu_name)

            !------------------------------------------------------------------------------!
            !    Open the patch file and read in all patches.                              !
            !------------------------------------------------------------------------------!
            open(unit=12,file=trim(lu_name),form='formatted',status='old',action='read')

            !------------------------------------------------------------------------------!
            !     Initialise the temporary arrays with data that will not cause any        !
            ! harvesting.  The actual variables will be read in the following block.       !
            !------------------------------------------------------------------------------!
            harvest_pft(1:n_pft)   = -1
            mindbh_slog(1:n_pft)   = huge_dbh
            harvprob_slog_g(1:n_pft) = 0.
            harvprob_slog_l(1:n_pft) = 0.
            mindbh_fplt(1:n_pft)   = huge_dbh
            harvprob_fplt_g(1:n_pft) = 0.
            harvprob_fplt_l(1:n_pft) = 0.
            
            !----- Initialise plantation patches if plantation information is available. --!
            cpoly%plantation(:) = 0
            call read_plantation_fractions(cpoly,cgrid%lon(ipy),cgrid%lat(ipy),igr)            


            !----- Define the format for the header. --------------------------------------!
            write(hform,fmt='(a,i3.3,a)') '(a',str_len,')'

            !----- Read the header. -------------------------------------------------------!
            read (unit=12,fmt=hform) cdum
            cdum = cdum(hoff:)
            read (cdum, fmt=*) wlon

            read (unit=12,fmt=hform) cdum
            cdum = cdum(hoff:)
            read (cdum, fmt=*) elon

            read (unit=12,fmt=hform) cdum
            cdum = cdum(hoff:)
            read (cdum, fmt=*) slat

            read (unit=12,fmt=hform) cdum
            cdum = cdum(hoff:)
            read (cdum, fmt=*) nlat

            read (unit=12,fmt=hform) cdum
            cdum = cdum(hoff:)
            read (cdum, fmt=*) lu_area

            read (unit=12,fmt=hform) cdum
            cdum = cdum(hoff:)
            read (cdum, fmt=*) yd_1st

            read (unit=12,fmt=hform) cdum
            cdum = cdum(hoff:)
            read (cdum, fmt=*) yd_last

            read (unit=12,fmt=hform) cdum
            cdum = cdum(hoff:)
            read (cdum, fmt=*) nharvest

            if (nharvest > 0 ) then
               !---------------------------------------------------------------------------!
               !     If nharvest is not 0, then the harvesting is not PFT-blind, read the  !
               ! PFT information in the file.                                              !
               !---------------------------------------------------------------------------!
               read (unit=12,fmt=hform)  cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) (harvest_pft(h)  ,h=1,nharvest)

               read (unit=12,fmt=hform)  cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) (mindbh_slog(h)  ,h=1,nharvest)

               read (unit=12,fmt=hform)  cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) (harvprob_slog_g(h)  ,h=1,nharvest)

               read (unit=12,fmt=hform)  cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) (harvprob_slog_l(h)  ,h=1,nharvest)

               read (unit=12,fmt=hform)  cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) (mindbh_fplt(h)      ,h=1,nharvest)

               read (unit=12,fmt=hform)  cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) (harvprob_fplt_g(h)  ,h=1,nharvest)

               read (unit=12,fmt=hform)  cdum
               cdum = cdum(hoff:)
               read (cdum, fmt=*) (harvprob_fplt_l(h)  ,h=1,nharvest)
            else
               !---------------------------------------------------------------------------!
               !     No specific PFT information was given, this is likely to be a case in !
               ! which the logging is based on absolute target biomass (PFT- and           !
               ! DBH-blind).                                                               !
               !---------------------------------------------------------------------------!
               h = 0
               nopftloop: do ipft=1,n_pft
                  h=h+1
                  harvest_pft(h)       = ipft
                  mindbh_slog(1:n_pft) = 0.
                  mindbh_fplt(1:n_pft) = 0.
               end do nopftloop
            end if
            read (unit=12,fmt=*) 

            !----- Use file_lat to compute the physical area sampled by the file. ---------!
            if (lu_area == 0.) then
               write (unit=*,fmt='(a)')           '---------------------------------------'
               write (unit=*,fmt='(2(a,1x))')     ' - File:    ',trim(lu_name)
               write (unit=*,fmt='(a,1x,es12.5)') ' - Wlon:    ',wlon
               write (unit=*,fmt='(a,1x,es12.5)') ' - Elon:    ',elon
               write (unit=*,fmt='(a,1x,es12.5)') ' - Slat:    ',slat
               write (unit=*,fmt='(a,1x,es12.5)') ' - Nlat:    ',nlat
               write (unit=*,fmt='(a,1x,es12.5)') ' - Lu_area: ',lu_area
               write (unit=*,fmt='(a,1x,i6)')     ' - Yd_1st:  ',yd_1st
               write (unit=*,fmt='(a,1x,i6)')     ' - Yd_last: ',yd_last
               write (unit=*,fmt='(a,1x,i6)')     ' - Nharvest:',nharvest
               write (unit=*,fmt='(a)')           '---------------------------------------'
               call fatal_error('Land use area is zero, it doesn''t make any sense!'       &
                               ,'landuse_init','landuse_init.f90')
            else
               lu_area_i = 1. / lu_area
            end if

            !----- Determine whether this block contains the current polygon. -------------!
            inside = cgrid%lon(ipy) >= wlon .and. cgrid%lon(ipy) <= elon .and.             &
                     cgrid%lat(ipy) >= slat .and. cgrid%lat(ipy) <= nlat


            !------------------------------------------------------------------------------!
            !     Here we will only use the land information if the polygon centre is      !
            ! inside the block of land use disturbance we are about to read.               !
            !------------------------------------------------------------------------------!
            if (inside) then
               !----- File exists, allocate the maximum number of years. ------------------!
               allocate(cpoly%clutimes(max_lu_years,cpoly%nsites))


               !----- Copy the file information to the first site. ------------------------!
               isi = 1
               csite => cpoly%site(isi)


               !----- Determine the number of disturbance years. --------------------------!
               cpoly%num_landuse_years(isi) = max(yd_last,iyearz)-min(yd_1st,iyeara) + 1


               !----- Initialise the PFT-dependent arrays. --------------------------------!
               cpoly%mindbh_harvest(1:n_pft,isi) = huge_dbh
               cpoly%prob_harvest_g(1:n_pft,isi) = 0.
               cpoly%prob_harvest_l(1:n_pft,isi) = 0.

			  ! From Marcos
			  !----- Fill the arrays with the appropriate PFT. ------------------------!
			  select case(cpoly%plantation(isi))
			  case (0)
				 harvloop_slog: do h=1,nharvest
					ipft = harvest_pft(h)
					if (ipft >= 1 .and. ipft <= n_pft) then
					   cpoly%mindbh_harvest(ipft,isi) = mindbh_slog  (h)
					   cpoly%prob_harvest_g  (ipft,isi) = harvprob_slog_g(h)
					   cpoly%prob_harvest_l  (ipft,isi) = harvprob_slog_l(h)
					end if
				 end do harvloop_slog
			  case (1)
				 harvloop_fplt: do h=1,nharvest
					ipft = harvest_pft(h)
					if (ipft >= 1 .and. ipft <= n_pft) then
					   cpoly%mindbh_harvest(ipft,isi) = mindbh_fplt  (h)
					   cpoly%prob_harvest_g  (ipft,isi) = harvprob_fplt_g(h)
					   cpoly%prob_harvest_l  (ipft,isi) = harvprob_fplt_l(h)
					end if
				 end do harvloop_fplt
			  end select
			  !------------------------------------------------------------------------!

! 
!                ! Original
!                !----- Fill the arrays with the appropriate PFT. ---------------------------!
!                do ipft=1,n_pft
!                   harvloop: do h=1,nharvest
!                      if (harvest_pft(h) == ipft) then
!                         cpoly%mindbh_primary    (ipft,isi) = mindbh_slog  (h)
!                         cpoly%probharv_primary  (ipft,isi) = harvprob_slog(h)
!                         cpoly%mindbh_secondary  (ipft,isi) = mindbh_fplt  (h)
!                         cpoly%probharv_secondary(ipft,isi) = harvprob_fplt(h)
!                         exit harvloop
!                      end if
!                   end do harvloop
!                end do


               !----- Padding disturbances with zero before first available lu year. ------!
               iyear = 0
               do yd_this = iyeara,(yd_1st-1)
                  iyear = iyear + 1
                  clutime => cpoly%clutimes(iyear,isi)

                  clutime%landuse_year            = yd_this
                  clutime%landuse(1:num_lu_trans) = 0.0
               end do

               !---- Reading the years that have data -------------------------------------!
               do yd_this = yd_1st,yd_last
                  iyear = iyear + 1
                  clutime => cpoly%clutimes(iyear,isi)

                  read(unit=12,fmt=*) clutime%landuse_year, clutime%landuse(1:num_lu_trans)
                  
                  !------------------------------------------------------------------------!
                  !    Here we normalise by the area, except when landuse(12) and/or       !
                  ! landuse(14) are negative (a special flag to use the selective logging. !
                  !------------------------------------------------------------------------!
                  if ( clutime%landuse(12) > 0. ) then
                     clutime%landuse(12) = lu_area_i * clutime%landuse(12)
                  end if
                  if ( clutime%landuse(14) > 0. ) then
                     clutime%landuse(14) = lu_area_i * clutime%landuse(14) 
                  end if
                  ! Note: These two columns aren't really used right now
                  clutime%landuse(16) = lu_area_i * clutime%landuse(16)
                  clutime%landuse(18) = lu_area_i * clutime%landuse(18)
               end do

               !----- Padding disturbances with zero after last available lu year. --------!
               do yd_this = (yd_last+1),iyearz
                  iyear = iyear + 1
                  clutime => cpoly%clutimes(iyear,isi)
                  clutime%landuse_year            = yd_this
                  clutime%landuse(1:num_lu_trans) = 0.0
               end do


               !---------------------------------------------------------------------------!
               !      Copy the information from the first site to the other, if they       !
               ! exist.                                                                    !
               !---------------------------------------------------------------------------!
               siteloop: do isi = 2,cpoly%nsites
                  csite => cpoly%site(isi)

                  !----- Determine the number of disturbance years. -----------------------!
                  cpoly%num_landuse_years(isi) = cpoly%num_landuse_years(1)

                  !----- PFT-dependent harvest characteristics. ------------------------!
                  cpoly%mindbh_harvest(:,isi) = cpoly%mindbh_harvest(:,1)
                  cpoly%prob_harvest_g(:,isi) = cpoly%prob_harvest_g(:,1)
                  cpoly%prob_harvest_l(:,isi) = cpoly%prob_harvest_l(:,1)

                  !----- Disturbances. ----------------------------------------------------!
                  do iyear = 1,cpoly%num_landuse_years(isi)
                     clutime   => cpoly%clutimes(iyear,isi)
                     onelutime => cpoly%clutimes(iyear,1)
                     clutime%landuse_year            = onelutime%landuse_year 
                     clutime%landuse(1:num_lu_trans) = onelutime%landuse(1:num_lu_trans)
                  end do

               end do siteloop
            else
               !---------------------------------------------------------------------------!
               !      No GLU data for this site.  Probably water.                          !
               !---------------------------------------------------------------------------!
               write (unit=*,fmt='(a)') '-------------------------------------------------'
               write (unit=*,fmt='(a)') ' The closest land use point is too far away...'
               write (unit=*,fmt='(a)') ' - File:                        ',trim(lu_name)
               write (unit=*,fmt=fffmt) ' - Polygon longitude:           ',cgrid%lon(ipy)
               write (unit=*,fmt=fffmt) ' - Polygon latitude:            ',cgrid%lat(ipy)
               write (unit=*,fmt=fffmt) ' - Closest LU central longitude:',llon_list(ncl)
               write (unit=*,fmt=fffmt) ' - Closest LU central latitude: ',llat_list(ncl)
               write (unit=*,fmt=fffmt) ' - Closest LU western edge:     ',wlon
               write (unit=*,fmt=fffmt) ' - Closest LU eastern edge:     ',elon
               write (unit=*,fmt=fffmt) ' - Closest LU southern edge:    ',slat
               write (unit=*,fmt=fffmt) ' - Closest LU northern edge:    ',nlat
               write (unit=*,fmt=esfmt) ' - Distance:                    ',file_ldist(ncl)
               write (unit=*,fmt=*)     ' '
               write (unit=*,fmt='(a)') ' We will assign no land use disturbance rate.'
               write (unit=*,fmt='(a)') '-------------------------------------------------'

               !----- Allocate just 1 landuse year. ---------------------------------------!
               allocate(cpoly%clutimes(1,cpoly%nsites))

               !----- Set the parameters in a way that no logging/ploughing will happen. --!
               do isi = 1,cpoly%nsites
                  cpoly%num_landuse_years(isi)                  = 1
                  cpoly%mindbh_harvest(1:n_pft,isi)             = huge_dbh
                  cpoly%prob_harvest_g(1:n_pft,isi)             = 0.
                  cpoly%prob_harvest_l(1:n_pft,isi)             = 0.
                  cpoly%clutimes(1,isi)%landuse_year            = iyeara
                  cpoly%clutimes(1,isi)%landuse(1:num_lu_trans) = 0.0
               end do
            end if
            !------------------------------------------------------------------------------!


            !----- Close the land use file, outside the if statement. ---------------------!
            close(unit=12,status='keep')
            !------------------------------------------------------------------------------!

!           case (2)
!              !---------------------------------------------------------------------------!
!              !      Make the land use data based on ED2IN.                               !
!              !      Work with the first site, then copy the data to the others.          !
!              !---------------------------------------------------------------------------!
!              isi   = 1
!              csite => cpoly%site(isi)
!              !---------------------------------------------------------------------------!
! 
! 
!              !----- No plantations. -----------------------------------------------------!
!              cpoly%plantation(:) = 0
!              !---------------------------------------------------------------------------!
! 
! 
! 
!              !----- Determine the number of disturbance years. --------------------------!
!              cpoly%num_landuse_years(isi) = sim_years
!              !---------------------------------------------------------------------------!
! 
! 
! 
!              !----- File exists, allocate the maximum number of years. ------------------!
!              allocate(cpoly%clutimes(sim_years,cpoly%nsites))
!              !---------------------------------------------------------------------------!
! 
! 
!              !----- Initialise the PFT-dependent arrays. --------------------------------!
!              cpoly%mindbh_primary    (1:n_pft,isi) = huge_dbh
!              cpoly%probharv_primary  (1:n_pft,isi) = 0.
!              cpoly%mindbh_secondary  (1:n_pft,isi) = huge_dbh
!              cpoly%probharv_secondary(1:n_pft,isi) = 0.
!              !---------------------------------------------------------------------------!
! 
! 
! 
!              !------ Find the number of PFT that can be harvested. ----------------------!
!              nharvest = count(sl_pft >= 1 .and. sl_pft <= n_pft)
!              !---------------------------------------------------------------------------!
! 
!              !----- Fill the arrays with the appropriate PFT. ---------------------------!
!              harvloop_two: do h=1,nharvest
!                 ipft = sl_pft(h)
!                 if (ipft >= 1 .and. ipft <= n_pft) then
!                    cpoly%mindbh_harvest(ipft,isi) = sl_mindbh_harvest(h)
!                    cpoly%prob_harvest  (ipft,isi) = sl_prob_harvest  (h)
!                 end if
!              end do harvloop_two
!              !---------------------------------------------------------------------------!
! 
! 
! 
!              !---------------------------------------------------------------------------!
!              !      Fill in the disturbance matrices and biomass target.                 !
!              !---------------------------------------------------------------------------!
!              iyear = 0
!              do yd_this = iyeara,iyearz
!                 iyear = iyear + 1
!                 clutime => cpoly%clutimes(iyear,isi)
! 
!                 clutime%landuse_year            = yd_this
!                 clutime%landuse(1:num_lu_trans) = 0.
! 
!                 !------------------------------------------------------------------------!
!                 !     Decide whether to include logging disturbance in this year.        !
!                 !------------------------------------------------------------------------!
!                 if (yd_this >= sl_yr_first) then
!                    if ( (sl_scale == 1) .or. (mod(yd_this-sl_yr_first,sl_nyrs) == 0) )   &
!                    then
!                       clutime%landuse(12) = -1.0
!                       clutime%landuse(14) = -1.0
!                    end if
!                 end if
!                 !------------------------------------------------------------------------!
!              end do
!              !---------------------------------------------------------------------------!
! 
! 
!              !---------------------------------------------------------------------------!
!              !      Copy the information from the first site to the others.              !
!              !---------------------------------------------------------------------------!
!              siteloop_two: do isi = 2,cpoly%nsites
!                 csite => cpoly%site(isi)
! 
!                 !----- Determine the number of disturbance years. -----------------------!
!                 cpoly%num_landuse_years(isi) = cpoly%num_landuse_years(1)
!                 !------------------------------------------------------------------------!
! 
! 
!                 !----- PFT-dependent harvest characteristics. ---------------------------!
!                 cpoly%mindbh_harvest(:,isi) = cpoly%mindbh_harvest(:,1)
!                 cpoly%prob_harvest  (:,isi) = cpoly%prob_harvest  (:,1)
!                 !------------------------------------------------------------------------!
! 
! 
! 
!                 !----- Disturbances. ----------------------------------------------------!
!                 do iyear = 1,cpoly%num_landuse_years(isi)
!                    clutime   => cpoly%clutimes(iyear,isi)
!                    onelutime => cpoly%clutimes(iyear,1)
!                    clutime%landuse_year            = onelutime%landuse_year 
!                    clutime%landuse(1:num_lu_trans) = onelutime%landuse(1:num_lu_trans)
!                 end do
!                 !------------------------------------------------------------------------!
! 
!              end do siteloop_two
! 
         end select
         
      end do polyloop
   end do gridloop

   return
end subroutine landuse_init
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine reads the plantation fraction.                                       !
!------------------------------------------------------------------------------------------!
subroutine read_plantation_fractions(cpoly,polylon,polylat,igr)
   use ed_state_vars , only : polygontype         ! ! structure
   use ed_max_dims   , only : str_len             ! ! intent(in)
   use disturb_coms  , only : plantation_file     & ! intent(in)
                            , max_plantation_dist & ! intent(in)
                            , min_plantation_frac ! ! intent(in)

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(polygontype)     , target      :: cpoly
   real                  , intent(in)  :: polylon
   real                  , intent(in)  :: polylat
   integer               , intent(in)  :: igr
   !----- Local variables -----------------------------------------------------------------!
   logical                             :: exans
   integer                             :: ierr
   integer                             :: isi
   real, dimension(:)    , allocatable :: plantlon
   real, dimension(:)    , allocatable :: plantlat
   real, dimension(:)    , allocatable :: plantdist
   real, dimension(:)    , allocatable :: fracplant
   real                                :: rdum
   integer                             :: ndat
   integer                             :: n
   !----- Local constants. ----------------------------------------------------------------!
   character(len=12)    , parameter    :: fffmt='(a,1x,f12.5)'
   character(len=13)    , parameter    :: esfmt='(a,1x,es12.5)'
   !----- External functions. -------------------------------------------------------------!
   real                  , external    :: dist_gc
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     If the user left the plantation file empty, it means that they don't want to      !
   ! use plantation files, skip the subroutine and don't print any warnings.               !
   !---------------------------------------------------------------------------------------!
   if (len_trim(plantation_file(igr)) == 0) then
      return
   else
      !----- Check whether plantation file exists. ----------------------------------------!
      inquire(file=trim(plantation_file(igr)),exist=exans)
      if (.not.exans)then
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') 'File :'//trim(plantation_file(igr))//' not found...'
         write (unit=*,fmt='(a)') 'Assuming that there are no plantations'
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         return
      end if
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     If we reach this point, then there is a plantation file. Open it.                 !
   !---------------------------------------------------------------------------------------!
   open (unit=12, file=trim(plantation_file(igr)), form='formatted', status='old')
   !----- The first loop will determine how many points are available. --------------------!
   ndat = 0
   countloop: do
      read (unit=12,fmt=*,iostat=ierr) rdum
      if (ierr /= 0) exit countloop
      
      ndat = ndat + 1
   end do countloop
   rewind(unit=12)

   !----- Now that we know the number of lines, allocate the temporary vectors. -----------!
   allocate (plantlon(ndat),plantlat(ndat),plantdist(ndat),fracplant(ndat))

   !----- Now we read the data. -----------------------------------------------------------!
   read_plantation: do n=1,ndat
      read (unit=12,fmt=*) plantlat(n), plantlon(n), fracplant(n)
      plantdist(n) = dist_gc(polylon,plantlon(n),polylat,plantlat(n))
   end do read_plantation
   close(unit=12,status='keep')

   !----- Determine which point was the closest one. --------------------------------------!
   n = minloc(plantdist,dim=1)

   !----- Check whether the distance is not too large.  If it is, don't use it. -----------!
   if (plantdist(n) > max_plantation_dist) then
      write (unit=*,fmt='(a)') '----------------------------------------------------------'
      write (unit=*,fmt=*)     ' The closest plantation point is too far away...'
      write (unit=*,fmt=fffmt) ' - Polygon longitude:                       ',polylon
      write (unit=*,fmt=fffmt) ' - Polygon latitude:                        ',polylat
      write (unit=*,fmt=fffmt) ' - Plantation closest longitude:            ',plantlon(n)
      write (unit=*,fmt=fffmt) ' - Plantation closest latitude:             ',plantlat(n)
      write (unit=*,fmt=esfmt) ' - Distance between polygon and plantation: ',plantdist(n)
      write (unit=*,fmt=esfmt) ' - Maximum accepted distance:               '              &
                                  ,max_plantation_dist
      write (unit=*,fmt=*)     ' '
      write (unit=*,fmt='(a)') ' We will assume no plantations in this polygon.'
      write (unit=*,fmt='(a)') '----------------------------------------------------------'

   elseif (fracplant(n) > min_plantation_frac) then
      !----- Assume that all sites are plantation sites. ----------------------------------!
      do isi = 1,cpoly%nsites
         cpoly%plantation(isi) = 1
      end do
   end if
   !---------------------------------------------------------------------------------------!

   deallocate(plantlon,plantlat,fracplant,plantdist)

   return
end subroutine read_plantation_fractions
!==========================================================================================!
!==========================================================================================!

end module landuse_init_module