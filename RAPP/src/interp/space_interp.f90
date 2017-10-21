!==========================================================================================!
!==========================================================================================!
!     space_interp.f90                                                                     !
!                                                                                          !
!     This file contains several functions and subroutines that are used during the        !
! spatial interpolation (aka objective analysis).                                          !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will map the position of all lon/lat grid points relative to the     !
! Gaussian grid, and leave the points outside the domain undefined.                        !
!------------------------------------------------------------------------------------------!
subroutine map_lolaxgauss()
   use mod_grid   , only : grid_g        ! ! intent(in)
   use mod_interp , only : interp_buffer & ! intent(inout)
                         , mxgauss       & ! intent(in)
                         , mygauss       & ! intent(in)
                         , mxlola        & ! intent(in)
                         , mylola        & ! intent(in)
                         , xla           & ! intent(out)
                         , xlz           & ! intent(out)
                         , yla           & ! intent(out)
                         , ylz           ! ! intent(out)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer        :: xg
   integer        :: yg
   integer        :: xl
   integer        :: yl
   real           :: dlongauss
   real           :: dlatgauss
   !----- External functions. -------------------------------------------------------------!
   real, external :: dist_gc
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     First we find the edges in which both grids overlap. Outside this region, we will !
   ! not bother looking for data.                                                          !
   !---------------------------------------------------------------------------------------!
   !------ Western boundary. --------------------------------------------------------------!
   fwestloop: do xla=1,mxlola
      if (grid_g(3)%lon(xla,1) >= grid_g(1)%lon(1,1)) exit fwestloop
   end do fwestloop
   !------ Eastern boundary. --------------------------------------------------------------!
   feastloop: do xlz=xla,mxlola
      if (grid_g(3)%lon(xlz,1) >= grid_g(1)%lon(mxgauss,1)) exit feastloop
   end do feastloop
   xlz = xlz - 1
   !------ Northern boundary. -------------------------------------------------------------!
   fnorthloop: do yla=1,mylola
      if (grid_g(3)%lat(1,yla) <= grid_g(1)%lat(1,1)) exit fnorthloop
   end do fnorthloop
   !------ Southern boundary. -------------------------------------------------------------!
   fsouthloop: do ylz=yla,mylola
      if (grid_g(3)%lat(1,ylz) <= grid_g(1)%lat(1,mygauss)) exit fsouthloop
   end do fsouthloop
   ylz = ylz - 1
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Now we store the x and y point immediately to the northwest of each Gaussian      !
   ! grid point.  Then we compute the relative position of the lon/lat grid in the box     !
   ! defined by ixe..ixe+1 and jno..jno+1, which we call dix and diy.                      !
   !---------------------------------------------------------------------------------------!
   ylloop: do yl=yla,ylz
      xlloop: do xl=xla,xlz
         !----- Mapping the longitude. ----------------------------------------------------!
         xgloop: do xg=2,mxgauss-1
            if (grid_g(1)%lon(xg,1) > grid_g(3)%lon(xl,yl)) exit xgloop
         end do xgloop
         xg = xg - 1
         interp_buffer%iwe(xl,yl) = xg
         !---------------------------------------------------------------------------------!


         !----- Mapping the latitude. -----------------------------------------------------!
         ygloop: do yg=2,mygauss-1
            if (grid_g(1)%lat(1,yg) < grid_g(3)%lat(xl,yl)) exit ygloop
         end do ygloop
         yg = yg - 1
         interp_buffer%jno(xl,yl) = yg
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Finding the inverse of distance between two consecutive Gaussian grid       !
         ! points at around the lon/lat coordinate.                                        !
         !---------------------------------------------------------------------------------!
         dlongauss = dist_gc(grid_g(1)%lon(xg,yg),grid_g(1)%lon(xg+1,yg)                   &
                            ,grid_g(3)%lat(xl,yl),grid_g(3)%lat(xl,yl) )
         dlatgauss = dist_gc(grid_g(3)%lon(xl,yl),grid_g(3)%lon(xl,yl)                     &
                            ,grid_g(1)%lat(xg,yg),grid_g(1)%lat(xg,yg+1) )
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Finding the relative distance in longitude and latitude between the        !
         ! Lon/Lat grid point and the "wall" at the north and the west.  This will be used !
         ! to refine the first pass objective analysis.                                    !
         !---------------------------------------------------------------------------------!
         interp_buffer%dix(xl,yl) = dist_gc(grid_g(1)%lon(xg,yg),grid_g(3)%lon(xl,yl)      &
                                           ,grid_g(3)%lat(xl,yl),grid_g(3)%lat(xl,yl) )    &
                                  / dlongauss

         interp_buffer%diy(xl,yl) = dist_gc(grid_g(3)%lon(xl,yl),grid_g(3)%lon(xl,yl)      &
                                           ,grid_g(1)%lat(xg,yg),grid_g(3)%lat(xl,yl) )    &
                                  / dlatgauss
         !---------------------------------------------------------------------------------!

      end do xlloop
   end do ylloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine map_lolaxgauss
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will compute some dimension dependent variables that are represent-  !
! ative of the regular Lon/Lat grid around a given point in the Gaussian grid.  After this !
! is found, then we estimate kappa0, which will control the decaying rate of weights as a  !
! function of distance.  In the beginning, we determine the best D0(2 Delta_n) based on    !
! the value of gamma provided by the user.                                                 !
!------------------------------------------------------------------------------------------!
subroutine assign_interp_dims()
   use mod_interp , only : interp_buffer & ! intent(inout)
                         , gamma0        & ! intent(in)
                         , mxgauss       & ! intent(in)
                         , mygauss       & ! intent(in)
                         , xla           & ! intent(out)
                         , xlz           & ! intent(out)
                         , yla           & ! intent(out)
                         , ylz           ! ! intent(out)
   use mod_grid   , only : grid_g        ! ! intent(in)
   use rconstants , only : pii           ! ! intent(in)
   
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer        :: x            !
   integer        :: y            !
   real           :: deltax       ! Delta x
   real           :: deltay       ! Delta y
   real           :: halfdlon     ! Half the average delta_lon of the Lon/Lat grid
   real           :: halfdlat     ! Half the average delta_lat of the Lon/Lat grid
   real           :: west         ! Longitude half delta_lon to the west
   real           :: east         ! Longitude half delta_lon to the east
   real           :: north        ! Latitude half delta_lon to the north
   real           :: south        ! Latitude half delta_lon to the south
   real           :: d02dtn       ! D0(2Delta_n) that makes D1* = exp(-1)
   real           :: lnd02dtn     ! ln(D0_2Delta_n)
   !----- External functions. -------------------------------------------------------------!
   real, external :: oa_firstpass ! 
   real, external :: dist_gc      ! Great circle distance.
   real, external :: d0barnes     ! Function that finds D0 based on Barnes (1973).
   !---------------------------------------------------------------------------------------!



   !------ Initialising some variables. ---------------------------------------------------!
   d02dtn   = d0barnes(gamma0)
   lnd02dtn = log(d02dtn)
   !---------------------------------------------------------------------------------------!

   !------ Checking the average grid spacing. ---------------------------------------------!
   if (xla < xlz .and. yla < ylz) then
      halfdlon = 0.5 * (grid_g(3)%lon(xlz,1)-grid_g(3)%lon(xla,1))/real(xlz-xla)
      halfdlat = 0.5 * (grid_g(3)%lat(1,yla)-grid_g(3)%lat(1,ylz))/real(ylz-yla)
   else
      call fatal_error('Your domain area is too small.   '                                 &
                     //'Try setting lonw/lone and latn/lats further apart'                 &
                      ,'assign_deltan_kappa0','interp_utils.f90')
   end if
   !---------------------------------------------------------------------------------------!


   !------ Internal area of the domain. ---------------------------------------------------!
   myloop: do y = 1, mygauss
      mxloop: do x = 1, mxgauss
         west  = grid_g(1)%lon(x,y) - halfdlon
         east  = grid_g(1)%lon(x,y) + halfdlon
         north = grid_g(1)%lat(x,y) + halfdlat
         south = grid_g(1)%lat(x,y) - halfdlat
         !----- Find a typical distance of a regular lon/lat grid at around point x,y. ---!
         deltax = dist_gc(west,east,grid_g(1)%lat(x,y),grid_g(1)%lat(x,y))
         deltay = dist_gc(grid_g(1)%lon(x,y),grid_g(1)%lon(x,y),south,north)
         interp_buffer%deltan(x,y) = 0.50 * (deltax+deltay)

         !----- Now that delta_n is defined, we assign kappa0. ----------------------------!
         interp_buffer%kappa0(x,y) = - (2 * interp_buffer%deltan(x,y) * pii)               &
                                     * (2 * interp_buffer%deltan(x,y) * pii)               &
                                     * lnd02dtn
      end do mxloop
   end do myloop
   !---------------------------------------------------------------------------------------!


   return
end subroutine assign_interp_dims
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine lola_2_gauss()
   use mod_interp , only : interp_buffer & ! intent(inout)
                         , minweight     & ! intent(in)
                         , gamma0        & ! intent(in)
                         , mxgauss       & ! intent(in)
                         , mygauss       & ! intent(in)
                         , mxlola        & ! intent(in)
                         , mylola        ! ! intent(in)
   use mod_grid   , only : grid_g        & ! intent(in)
                         , ssxp          & ! intent(in)
                         , ssyp          & ! intent(in)
                         , sstp          ! ! intent(in)
   use mod_ncep   , only : ncep_g        ! ! intent(inout)
   use therm_lib  , only : rslif         ! ! function
   use rconstants , only : toodry        ! ! intent(in)

   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer        :: xx         ! Counter
   integer        :: yy         ! Counter
   integer        :: tt         ! Counter
   real           :: rsat       ! Saturation mixing ratio.
   real           :: qsat       ! Saturation specific humidity.
   !---------------------------------------------------------------------------------------!


   write (unit=*,fmt='(a)') '     - Running the objective analysis...'


   !---------------------------------------------------------------------------------------!
   !    For all variables we need to interpolate, we will perform the objective analysis   !
   ! in two steps: 
   !   STEP 1.  Find all weights and mask, then compute the first pass of our objective    !
   !            analysis, a simple weighted average.                                       !
   !   STEP 2.  Assess the residual of the first pass and run the next iterative step.     !
   !---------------------------------------------------------------------------------------!
   !----- Mixing ratio. -------------------------------------------------------------------!
   call oa_1stpass(sstp(3),'Specific humidity',grid_g(3)%lon,grid_g(3)%lat,ncep_g(3)%shum  &
                                              ,grid_g(1)%lon,grid_g(1)%lat,ncep_g(1)%shum )
   call oa_2ndpass(sstp(3),'Specific humidity',grid_g(3)%lon,grid_g(3)%lat,ncep_g(3)%shum  &
                                              ,grid_g(1)%lon,grid_g(1)%lat,ncep_g(1)%shum )
   !----- Ice-liquid potential temperature. -----------------------------------------------!
   call oa_1stpass(sstp(3),'Temperature      ',grid_g(3)%lon,grid_g(3)%lat,ncep_g(3)%temp  &
                                              ,grid_g(1)%lon,grid_g(1)%lat,ncep_g(1)%temp )
   call oa_2ndpass(sstp(3),'Temperature      ',grid_g(3)%lon,grid_g(3)%lat,ncep_g(3)%temp  &
                                              ,grid_g(1)%lon,grid_g(1)%lat,ncep_g(1)%temp )
   !----- Zonal wind (for wind direction only). -------------------------------------------!
   call oa_1stpass(sstp(3),'Zonal wind       ',grid_g(3)%lon,grid_g(3)%lat,ncep_g(3)%uwnd  &
                                              ,grid_g(1)%lon,grid_g(1)%lat,ncep_g(1)%uwnd )
   call oa_2ndpass(sstp(3),'Zonal wind       ',grid_g(3)%lon,grid_g(3)%lat,ncep_g(3)%uwnd  &
                                              ,grid_g(1)%lon,grid_g(1)%lat,ncep_g(1)%uwnd )
   !----- Meridional wind (for wind direction only). --------------------------------------!
   call oa_1stpass(sstp(3),'Meridional wind  ',grid_g(3)%lon,grid_g(3)%lat,ncep_g(3)%vwnd  &
                                              ,grid_g(1)%lon,grid_g(1)%lat,ncep_g(1)%vwnd )
   call oa_2ndpass(sstp(3),'Meridional wind  ',grid_g(3)%lon,grid_g(3)%lat,ncep_g(3)%vwnd  &
                                              ,grid_g(1)%lon,grid_g(1)%lat,ncep_g(1)%vwnd )
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     We will now constrain specific humidity to avoid unrealistic values.  The object- !
   ! ive analysis is not positively defined, so some spurious negative values may appear.  !
   ! Likewise, we may have some spurious super-saturation because we interpolated both     !
   ! temperature and specific humidity.                                                    !
   !---------------------------------------------------------------------------------------!
   ttloop: do tt=1,sstp(1)
      yyloop: do yy=1,mygauss
         xxloop: do xx=1,mxgauss
            rsat = rslif(ncep_g(1)%pres(xx,yy,tt),ncep_g(1)%temp(xx,yy,tt))
            qsat = rsat / (1. + rsat)
            ncep_g(1)%shum(xx,yy,tt) = min(qsat,max(toodry,ncep_g(1)%shum(xx,yy,tt)))
         end do xxloop
      end do yyloop
   end do ttloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine lola_2_gauss
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will perform the first pass of the objective analysis.                !
!------------------------------------------------------------------------------------------!
subroutine oa_1stpass(mtp,vname,lonlola,latlola,vlola,longauss,latgauss,vgauss)
   use mod_interp , only : interp_buffer & ! intent(inout)
                         , minweight8    & ! intent(in)
                         , mxgauss       & ! intent(in)
                         , mygauss       & ! intent(in)
                         , mxlola        & ! intent(in)
                         , mylola        & ! intent(in)
                         , xla           & ! intent(in)
                         , xlz           & ! intent(in)
                         , yla           & ! intent(in)
                         , ylz           ! ! intent(in)
   use mod_ioopts , only : missflg_real  & ! intent(in)
                         , edgeoff       ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                                , intent(in)    :: mtp      ! # of time steps
   character(len=*)                       , intent(in)    :: vname    ! Variable name
   real   , dimension(mxlola ,mylola )    , intent(in)    :: lonlola  ! Long.    (lon/lat)
   real   , dimension(mxlola ,mylola )    , intent(in)    :: latlola  ! Latitude (lon/lat)
   real   , dimension(mxlola ,mylola ,mtp), intent(in)    :: vlola    ! Variable (lon/lat)
   real   , dimension(mxgauss,mygauss)    , intent(in)    :: longauss ! Long.    (Gaussian)
   real   , dimension(mxgauss,mygauss)    , intent(in)    :: latgauss ! Latitude (Gaussian)
   real   , dimension(mxgauss,mygauss,mtp), intent(inout) :: vgauss   ! Variable (Gaussian)
   !----- Local variables. ----------------------------------------------------------------!
   integer              :: xl       ! Counter
   integer              :: yl       ! Counter
   integer              :: xg       ! Counter
   integer              :: yg       ! Counter
   integer              :: tt       ! Counter
   integer              :: xw       ! Westernmost point for output
   integer              :: xe       ! Easternmost point for output
   integer              :: yn       ! Northernmost point for output
   integer              :: ys       ! Southernmost point for output
   integer              :: misscnt  ! Number of missing points in this domain.
   real(kind=8)         :: weight8  ! Weight in double precision
   real                 :: distance ! Distance between Lon/Lat and Gaussian grid points.
   real                 :: varmin   ! Minimum value after the objective analysis
   real                 :: varmax   ! Maximum value after the objective analysis
   !----- External functions. -------------------------------------------------------------!
   real   , external    :: dist_gc  ! Great circle distance
   real   , external    :: wei_ave  ! Weighted average
   !---------------------------------------------------------------------------------------!


   !------ Gaussian grid loop. ------------------------------------------------------------!
   ygloop: do yg=1,mygauss
      xgloop: do xg=1,mxgauss

         !---------------------------------------------------------------------------------!
         !   Lon/Lat grid loop.  Here we will set up some parameters used in the objective !
         ! analysis, namely the distance-dependent weight and the mask to skip points that !
         ! are too far from the grid point.                                                !
         !---------------------------------------------------------------------------------!
         ylloop: do yl=yla,ylz
            xlloop: do xl=xla,xlz
               !----- 1. Find the square of distance. -------------------------------------!
               distance = dist_gc(longauss(xg,yg),lonlola(xl,yl)                           &
                                 ,latgauss(xg,yg),latlola(xl,yl)                           )
               interp_buffer%rm2(xl,yl) = distance * distance


               !----- 2. Find the mask so we can skip some points... ----------------------!
               weight8 = dexp(- dble(interp_buffer%rm2(xl,yl))                             &
                              / dble(interp_buffer%kappa0(xg,yg)))
               interp_buffer%mask(xl,yl)   = weight8 > minweight8

               !----- 3. Find the weights of objective analysis. --------------------------!
               if (interp_buffer%mask(xl,yl)) then
                  interp_buffer%weight(xl,yl) = sngl(weight8)
               else
                  interp_buffer%weight(xl,yl) = 0.
               end if
           end do xlloop
         end do ylloop

         !---------------------------------------------------------------------------------!
         !    Now we find the first pass of the objective analysis, which is just an       !
         ! weighted average.                                                               !
         !---------------------------------------------------------------------------------!
         ttloop: do tt=1,mtp
            vgauss(xg,yg,tt) = wei_ave(mxlola,mylola,xla,xlz,yla,ylz                       &
                                      ,vlola(1:mxlola,1:mylola,tt)                         &
                                      ,interp_buffer%weight,interp_buffer%mask)
         end do ttloop

      end do xgloop
   end do ygloop

   !----- Quick statistics to entretain the user and warn about possible problems. --------!
   varmin  = minval(vlola(xla:xlz,yla:ylz,:),mask= vlola(xla:xlz,yla:ylz,:) /= missflg_real)
   varmax  = maxval(vlola(xla:xlz,yla:ylz,:),mask= vlola(xla:xlz,yla:ylz,:) /= missflg_real)
   misscnt = count(vlola(xla:xlz,yla:ylz,:) == missflg_real)
   write (unit=*,fmt='(2(a,1x),2(a,1x,es14.7,1x),a,1x,i6,a)')                              &
                               '         [|] Objective analysis, reference  :',vname       &
                              ,'. Range: [',varmin,':',varmax,'] . # missing: ',misscnt,'.'
   !----- Quick statistics to entretain the user and warn about possible problems. --------!
   xw = 1 + edgeoff
   xe = mxgauss - edgeoff
   yn = 1 + edgeoff
   ys = mygauss - edgeoff
   varmin = minval(vgauss(xw:xe,yn:ys,:),mask=vgauss(xw:xe,yn:ys,:) /= missflg_real)
   varmax = maxval(vgauss(xw:xe,yn:ys,:),mask=vgauss(xw:xe,yn:ys,:) /= missflg_real)
   misscnt = count(vgauss(xw:xe,yn:ys,:) == missflg_real)

   write (unit=*,fmt='(2(a,1x),2(a,1x,es14.7,1x),a,1x,i6,a)')                              &
                               '         [|] Objective analysis, 1st pass   :',vname       &
                              ,'. Range: [',varmin,':',varmax,'] . # missing: ',misscnt,'.'

   return
end subroutine oa_1stpass
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will assess the residual error of the previous pass, then compute    !
! the objective analysis incorporating the residual.                                       !
!------------------------------------------------------------------------------------------!
subroutine oa_2ndpass(mtp,vname,lonlola,latlola,vlola,longauss,latgauss,vgauss)
   use mod_interp , only : interp_buffer & ! intent(inout)
                         , minweight8    & ! intent(in)
                         , gamma0        & ! intent(in)
                         , xla           & ! intent(in)
                         , xlz           & ! intent(in)
                         , yla           & ! intent(in)
                         , ylz           & ! intent(in)
                         , mxgauss       & ! intent(in)
                         , mygauss       & ! intent(in)
                         , mxlola        & ! intent(in)
                         , mylola        ! ! intent(in)
   use mod_ioopts , only : missflg_real  & ! intent(in)
                         , edgeoff       ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                                , intent(in)    :: mtp      ! # of time steps
   character(len=*)                       , intent(in)    :: vname    ! Variable name
   real   , dimension(mxlola ,mylola )    , intent(in)    :: lonlola  ! Long.    (lon/lat)
   real   , dimension(mxlola ,mylola )    , intent(in)    :: latlola  ! Latitude (lon/lat)
   real   , dimension(mxlola ,mylola ,mtp), intent(in)    :: vlola    ! Variable (lon/lat)
   real   , dimension(mxgauss,mygauss)    , intent(in)    :: longauss ! Long.    (Gaussian)
   real   , dimension(mxgauss,mygauss)    , intent(in)    :: latgauss ! Latitude (Gaussian)
   real   , dimension(mxgauss,mygauss,mtp), intent(inout) :: vgauss   ! Variable (Gaussian)
   !----- Local variables. ----------------------------------------------------------------!
   integer              :: xl       ! Counter for lon/lat grid
   integer              :: yl       ! Counter for lon/lat grid
   integer              :: xg       ! Counter for Gaussian grid
   integer              :: yg       ! Counter for Gaussian grid
   integer              :: tt       ! Counter for time
   integer              :: xw       ! Mapping - point immediately to the west
   integer              :: xe       ! Mapping - point immediately to the east
   integer              :: yn       ! Mapping - point immediately to the north
   integer              :: ys       ! Mapping - point immediately to the south
   integer              :: misscnt  ! Number of missing points in this domain.
   real(kind=8)         :: weight8  ! Weight in double precision
   real(kind=8)         :: kappan   ! Kappa of the nth iteration
   real                 :: voan     ! Interpolated variable just north of the point
   real                 :: voas     ! Interpolated variable just south of the point
   real                 :: vest     ! Estimated value using a bilinear interpolation.
   real                 :: distance ! Distance between Lon/Lat and Gaussian grid points.
   real                 :: varmin   ! Minimum value after the objective analysis
   real                 :: varmax   ! Maximum value after the objective analysis
   !----- External functions. -------------------------------------------------------------!
   real   , external    :: dist_gc  ! Great circle distance
   real   , external    :: wei_ave  ! Weighted average
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! A.  We now estimate the residual of the previous pass, by comparing the Lon/Lat grid  !
   !     with a simple bi-linear interpolation of the previous pass result.                !
   !---------------------------------------------------------------------------------------!
   interp_buffer%residu = 0.
   ryloop: do yl=yla,ylz
      rxloop: do xl=xla,xlz

         !---------------------------------------------------------------------------------!
         ! STEP 1. Find aliases for mapping indices.                                       !
         !---------------------------------------------------------------------------------!
         xw = interp_buffer%iwe(xl,yl)
         xe = xw + 1
         yn = interp_buffer%jno(xl,yl) 
         ys = yn + 1

         !---------------------------------------------------------------------------------!
         ! STEP 2.  We now find the residual for each variable, every time.  The way we do !
         !          this is by making a simple bi-linear interpolation using the four      !
         !          points around the Longitude/Latitude point.                            !
         !---------------------------------------------------------------------------------!
         rtloop: do tt=1,mtp
            !------------------------------------------------------------------------------!
            ! 2.a. Compute the X-direction interpolation.  This is done by interpolating   !
            !      the Gaussian points northwest and northeast of the current Lon/Lat      !
            !      point to the same longitude as the Lon/Lat point.  A similar procedure  !
            !      is done for the points to the southwest and southeast.                  !
            !------------------------------------------------------------------------------!
            voan =       interp_buffer%dix(xl,yl)  * vgauss(xw,yn,tt)                      &
                 + (1. - interp_buffer%dix(xl,yl)) * vgauss(xe,yn,tt)
            voas =       interp_buffer%dix(xl,yl)  * vgauss(xw,ys,tt)                      &
                 + (1. - interp_buffer%dix(xl,yl)) * vgauss(xe,ys,tt)

            !------------------------------------------------------------------------------!
            ! 2.b. Now that we reduced to one dimension, make a simple linear interpol-    !
            !      ation to find the estimated value at the Lon/Lat point.                 !
            !------------------------------------------------------------------------------!
            vest =       interp_buffer%diy(xl,yl)  * voan                                  &
                 + (1. - interp_buffer%diy(xl,yl)) * voas

            !------------------------------------------------------------------------------!
            ! 2.c. Compute the "residual".  It is actually what is missing to reach the    !
            !      expected value, and it is simply the difference between the expected    !
            !      and the estimated values.                                               ! 
            !------------------------------------------------------------------------------!
            if (vlola(xl,yl,tt) == missflg_real) then
               interp_buffer%residu(xl,yl,tt) = missflg_real
            else
               interp_buffer%residu(xl,yl,tt) = vlola(xl,yl,tt) - vest
            end if
         end do rtloop

      end do rxloop
   end do ryloop
   
   !---------------------------------------------------------------------------------------!
   ! B.  We will now add this nth pass to the Gaussian grid point, which will include      !
   !     the part not accounted by the previous passes.  These steps are very similar to   !
   !     the first pass, except that the kappa is scaled down, and the weighted average is !
   !     applied to the "residual" rather than the actual values.                          !
   !---------------------------------------------------------------------------------------!
   !------ Gaussian grid loop. ------------------------------------------------------------!
   ygloop: do yg=1,mygauss
      xgloop: do xg=1,mxgauss

         !---------------------------------------------------------------------------------!
         ! STEP 3.  Lon/Lat grid loop.  Here we will set up some parameters used in the    !
         !          objective analysis, namely the distance-dependent weight and the mask  !
         !          to skip points that are too far from the grid point.                   !
         !---------------------------------------------------------------------------------!
         ylloop: do yl=yla,ylz
            xlloop: do xl=xla,xlz
               !---------------------------------------------------------------------------!
               ! 3.a. Find the square of distance.                                         !
               !---------------------------------------------------------------------------!
               distance = dist_gc(longauss(xg,yg),lonlola(xl,yl)                           &
                                 ,latgauss(xg,yg),latlola(xl,yl)                           )
               interp_buffer%rm2(xl,yl) = distance * distance


               !---------------------------------------------------------------------------!
               ! 3.b. Find the mask so we can skip some points. The mask is always done    !
               !      using the coarsest weight ..                                         !
               !---------------------------------------------------------------------------!
               kappan  = dble(gamma0) * dble(interp_buffer%kappa0(xg,yg))
               weight8 = dexp(- dble(interp_buffer%rm2(xl,yl)) / kappan)
               interp_buffer%mask(xl,yl)   = weight8 > minweight8


               !---------------------------------------------------------------------------!
               ! 3.c. Find the weights of the response function.                           !
               !---------------------------------------------------------------------------!
               if (interp_buffer%mask(xl,yl)) then
                  interp_buffer%weight(xl,yl) = sngl(weight8)
               else
                  !----- Too far away, don't consider this point. -------------------------!
                  interp_buffer%weight(xl,yl) = 0.
               end if

            end do xlloop
         end do ylloop

         !---------------------------------------------------------------------------------!
         ! STEP 4.  Now we find the first pass of the objective analysis, which is just an !
         !          weighted average.                                                      !
         !---------------------------------------------------------------------------------!
         ttloop: do tt=1,mtp
            vgauss(xg,yg,tt) = vgauss(xg,yg,tt)                                            &
                             + wei_ave(mxlola,mylola,xla,xlz,yla,ylz                       &
                                      ,interp_buffer%residu(:,:,tt)                        &
                                      ,interp_buffer%weight,interp_buffer%mask)
         end do ttloop
      end do xgloop
   end do ygloop

   !---------------------------------------------------------------------------------------!
   ! C. Quick statistics to entretain the user and warn about possible problems.           !
   !---------------------------------------------------------------------------------------!
   xw = 1 + edgeoff
   xe = mxgauss - edgeoff
   yn = 1 + edgeoff
   ys = mygauss - edgeoff
   varmin = minval(vgauss(xw:xe,yn:ys,:),mask=vgauss(xw:xe,yn:ys,:) /= missflg_real)
   varmax = maxval(vgauss(xw:xe,yn:ys,:),mask=vgauss(xw:xe,yn:ys,:) /= missflg_real)
   misscnt = count(vgauss(xw:xe,yn:ys,:) == missflg_real)

   write (unit=*,fmt='(2(a,1x),2(a,1x,es14.7,1x),a,1x,i6,a)')                              &
                               '         [|] Objective analysis, second pass:',vname       &
                              ,'. Range: [',varmin,':',varmax,'] . # missing: ',misscnt,'.'


   return
end subroutine oa_2ndpass
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function simply computes the weighted average of an array, skipping the points  !
! we do not want to include.                                                               !
!------------------------------------------------------------------------------------------!
real function wei_ave(mxp,myp,xa,xz,ya,yz,var,weight,mask)
    use mod_ioopts , only : missflg_real  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                    , intent(in) :: mxp
   integer                    , intent(in) :: myp
   integer                    , intent(in) :: xa
   integer                    , intent(in) :: xz
   integer                    , intent(in) :: ya
   integer                    , intent(in) :: yz
   real   , dimension(mxp,myp), intent(in) :: var
   real   , dimension(mxp,myp), intent(in) :: weight
   logical, dimension(mxp,myp), intent(in) :: mask
   !----- Local variables. ----------------------------------------------------------------!
   integer           :: x
   integer           :: y
   real              :: wsum
   real              :: fwsum
   !---------------------------------------------------------------------------------------!


   !----- Initialise the weighted sum. ----------------------------------------------------!
   wsum  = 0.
   fwsum = 0.

   !----- Loop over every point and add the points that contribute to the Obj. Analysis. --!
   yloop: do y=ya,yz
      xloop: do x=xa,xz
         !----- Skip this point if it has a tiny influence on the total weight. -----------!
         if ((.not. mask(x,y)) .or. var(x,y) == missflg_real) cycle xloop

         wsum  = wsum  + weight(x,y)
         fwsum = fwsum + weight(x,y) * var(x,y)
      end do xloop
   end do yloop

   !----- Find the first pass, a simple weighted average. ---------------------------------!
   if (wsum > 0.) then
      wei_ave = fwsum / wsum
   else
      wei_ave = 0.
   end if

   return
end function wei_ave
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function will define the best D0(2 Delta_n for the Barnes objective analysis as !
! in Koch et al. (1983).  The value of D0 is chosen to give D1*=exp(-1) for a given gamma. !
!        D1* = D0 * (1 + D0^(gamma-1) - D0^gamma) (Koch et al, 1983, equation 11.)         !
!------------------------------------------------------------------------------------------!
real function d0barnes(gamma0)
   use therm_lib, only : toler    & ! intent(in)
                       , maxfpo   ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real   , intent(in) :: gamma0 ! Parameter gamma
   !----- Local variables. ----------------------------------------------------------------!
   real                :: deriv      ! Function derivative                         [   ---]
   real                :: fun        ! Function for which we seek a root.          [   ---]
   real                :: funa       ! Smallest  guess function                    [   ---]
   real                :: funz       ! Largest   guess function                    [   ---]
   real                :: d01st      ! First guess.                                [   ---]
   real                :: d0a        ! Smallest guess (or previous guess)          [   ---]
   real                :: d0z        ! Largest   guess (or new guess in Newton)    [   ---]
   real                :: delta      ! Aux. var --- 2nd guess for bisection        [   ---]
   real                :: gammam1    ! Gamma - 1.                                  [   ---]
   real                :: gammap1    ! Gamma + 1.                                  [   ---]
   real                :: expm1      ! exp(-1.)
   integer             :: itn,itb    ! Iteration counter                           [   ---]
   logical             :: converged  ! Convergence handle                          [   ---]
   logical             :: zside      ! Flag to check for one-sided approach...     [   ---]
   !---------------------------------------------------------------------------------------!


   !----- Assigning some constants. -------------------------------------------------------!
   expm1   = exp(-1.)
   gammam1 = gamma0 - 1.
   gammap1 = gamma0 + 1.

   !---------------------------------------------------------------------------------------! 
   !     Here we need to be careful when we choose the first guess. This function will     !
   ! fail if the guess becomes negative, so we must ensure that we chose a good first      !
   ! guess. If not, we make the first guess smaller until Newton's method moves away from  !
   ! 0.  This is a very specific characteristic with this function.                        !
   !---------------------------------------------------------------------------------------! 
   d01st = min(gamma0,0.99)
   d0a   = -1.
   deriv = 0. 
   do while (d0a < d01st .or. abs(deriv) < toler)
      !----- First guess and function/derivative evaluation. ------------------------------!
      d01st = 0.1 * d01st 
      funa  = d01st + d01st**gamma0 - d01st**gammap1 - expm1
      deriv =  1. + gamma0 * d01st**gammam1 - gammap1 * d01st**gamma0
      d0a   = d01st - funa / deriv

   end do
   !----- Copying just in case it fails at the first iteration. ------------------------!
   funa  = d0a + d0a**gamma0 - d0a**gammap1 - expm1
   deriv =  1. + gamma0 * d0a**gammam1 - gammap1 * d0a**gamma0
   d0z   = d0a
   fun   = funa
   !---------------------------------------------------------------------------------------! 



   !----- Enter Newton's method loop: -----------------------------------------------------!
   converged = .false.
   newloop: do itn = 1,maxfpo/6
      if (abs(deriv) < toler) exit newloop !----- Too dangerous, go with regula falsi -----!
      !----- Copying the previous guess ---------------------------------------------------!
      d0a   = d0z
      funa  = fun
      !----- New guess, its function and derivative evaluation ----------------------------!
      d0z   = d0a - fun/deriv
      fun   = d0z + d0z**gamma0 - d0z**gammap1 - expm1
      deriv = 1. + gamma0 * d0z**gammam1 - gammap1 * d0z**gamma0
      
      converged = abs(d0a-d0z) < toler * d0z
      if (converged) then
         d0barnes = 0.5 * (d0a + d0z)
         return
      elseif (fun == 0.) then
         d0barnes = d0z
         return
      end if
   end do newloop
   !---------------------------------------------------------------------------------------!
   !     If I reached this point then it's because Newton's method failed. Using regula    !
   ! falsi instead.                                                                        !
   !---------------------------------------------------------------------------------------!
   if (funa * fun < 0.) then
      funz  = fun
      zside = .true.
   else
      if (abs(fun-funa) < 100.*toler*d0a) then
         delta = 100.*toler*delta
      else
         delta = max(abs(funa * (d0z-d0a)/(fun-funa)),100.*toler*d0a)
      end if
      d0z = d0a + delta
 
      zside = .false.
      zgssloop: do itb=1,maxfpo
         d0z   = d0a + real((-1)**itb * (itb+3)/2) * delta
         funz  = d0z + d0z**gamma0 - d0z**gammap1 - expm1
         zside = funa*funz < 0
         if (zside) exit zgssloop
      end do zgssloop
      if (.not. zside) then
         write (unit=*,fmt='(a)') ' No second guess for you...'
         write (unit=*,fmt='(2(a,1x,es14.7))') 'd0a=',d0a,'funa=',funa
         write (unit=*,fmt='(2(a,1x,es14.7))') 'd0z=',d0z,'func=',funz
         call abort_run('Failed finding the second guess for regula falsi'                 &
                       ,'d0barnes','interp_utils.f90')
      end if
   end if

   fpoloop: do itb=itn,maxfpo
      d0barnes =  (funz*d0a-funa*d0z)/(funz-funa)

      !------------------------------------------------------------------------------------!
      !     Now that we updated the guess, check whether they are really close. If so, it  !
      ! converged, I can use this as my guess.                                             !
      !------------------------------------------------------------------------------------!
      converged = abs(d0barnes-d0a) < toler * d0barnes
      if (converged) exit fpoloop

      !------ Finding the new function ----------------------------------------------------!
      fun       =  d0barnes + d0barnes**gamma0 - d0barnes**gammap1 - expm1

      !------ Defining my new interval based on the intermediate value theorem. -----------!
      if (fun*funa < 0. ) then
         d0z   = d0barnes
         funz  = fun
         !----- If we are updating zside again, modify aside (Illinois method) ------------!
         if (zside) funa = funa * 0.5
         !----- We just updated zside, setting zside to true. -----------------------------!
         zside = .true.
      else
         d0a    = d0barnes
         funa   = fun
         !----- If we are updating aside again, modify aside (Illinois method) ------------!
         if (.not. zside) funz = funz * 0.5
         !----- We just updated aside, setting aside to true. -----------------------------!
         zside = .false.
      end if
   end do fpoloop

   if (.not. converged) then
      call abort_run('D0(2 Delta_n) didn''t converge, giving up!!!'                        &
                    ,'d0barnes','interp_utils.f90')
   end if

   return
end function d0barnes
!==========================================================================================!
!==========================================================================================!
