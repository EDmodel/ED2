!==========================================================================================!
!==========================================================================================!
!  RAPP. module model: contains variables from the input model.                            !
!------------------------------------------------------------------------------------------!
module mod_model

   use mod_maxdims , only : maxstr    & ! intent(in)
                          , maxgrds   & ! intent(in)
                          , nzpmax    & ! intent(in)
                          , nxpmax    & ! intent(in)
                          , nypmax    & ! intent(in)
                           ,maxtimes  ! ! intent(in)
   use mod_time    , only : time_stt  ! ! structure
   implicit none

   
   !----- Variables from the header -------------------------------------------------------!
   integer                            :: ihtran  ! Flag for horizontal coordinate:        
                                                 !   0 - Lon/Lat;
                                                 !   1 - Polar-stereographic;
                                                 !   2 - Lambert;
   !
   real                               :: polelon ! It depends on ihtran:
                                                 !   0 - Lon. @ NW corner;
                                                 !   1 - Pole longitude
                                                 !   2 - Longitude parallel to columns
   real                               :: polelat ! It depends on ihtran:
                                                 !   0 - Lat. @ NW corner;
                                                 !   1 - Pole latitude
                                                 !   2 - Lambert's first std. latitude
   real                               :: l2ndlat ! Used only if ihtran=2, it's Lambert's
                                                 !   second standard latitude

   real   , dimension(maxgrds)        :: centlon ! Longitude of the grid centre
   real   , dimension(maxgrds)        :: centlat ! Latitude of the grid centre
   !
   real   , dimension(maxgrds)        :: deltaxn ! Grid-dependent delta x
   real   , dimension(maxgrds)        :: deltayn ! Grid-dependent delta y

   !----- Full dataset information. -------------------------------------------------------!
   integer, dimension(maxgrds)        :: nnxp    ! Number of points in the X direction;
   integer, dimension(maxgrds)        :: nnyp    ! Number of points in the Y direction;
   integer, dimension(maxgrds)        :: nnzp    ! Number of points in the Z direction;
   integer, dimension(maxgrds)        :: nntp    ! Number of points in time;

   !----- Number of input grids -----------------------------------------------------------!
   integer                            :: ngrids
  
   !----- Some other potential dimensions (not in use). -----------------------------------!
   integer                            :: nzg     ! # of soil layers
   integer                            :: nzs     ! # of snow layers
   integer                            :: npatch  ! # of patches 
   integer                            :: nclouds ! # of clouds
   integer                            :: nwave   ! # of wave numbers
  
   !----- Vertical levels associated to Polar-Stereographic and vertical coordinates-------!
   real   , dimension(nxpmax,maxgrds) :: xtn     ! x coordinate of cell centre on PS
   real   , dimension(nxpmax,maxgrds) :: xmn     ! x coordinate of higher cell boundary 
   real   , dimension(nypmax,maxgrds) :: ytn     ! y coordinate of cell centre on PS
   real   , dimension(nypmax,maxgrds) :: ymn     ! y coordinate of higher cell boundary 
   real   , dimension(nzpmax,maxgrds) :: ztn     ! z coordinate of cell centre on PS
   real   , dimension(nzpmax,maxgrds) :: zmn     ! z coordinate of higher cell boundary 

   !----- Scratch Time strucures ----------------------------------------------------------!
   type(time_stt), dimension(maxtimes) :: this_time   ! Scratch to hold file time
   type(time_stt)                      :: zero_time   ! Scratch to hold init. time

end module mod_model
!==========================================================================================!
!==========================================================================================!
