!==========================================================================================!
!==========================================================================================!
!     This module contains several PFT-indepenendent parameters that control leaf          !
! phenology.  The PFT-dependent paramaters are all in pft_coms.                            !
!                                                                                          !
!    Kind reminder: DO NOT INITIALIZE NON-PARAMETERS IN THEIR MODULES - not all compilers  !
!                   will actually initialize them.  The subroutine "init_phen_coms" in     !
!                   ed_params.f90 will assign initial values.                              !
!------------------------------------------------------------------------------------------!
module phenology_coms

   use ed_max_dims, only: str_len ! ! intent(in)

   implicit none


   !=======================================================================================!
   !=======================================================================================!
   !    The following variables are initialised by the namelist, not by ed_params.f90.     !
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !<IPHEN_SCHEME -- It controls how phenology is updated for non-evergreen PFTs.         \n
   !<                This variable NO LONGER controls which leaf phenology strategy each  \n
   !<                PFT adopts. This is defined directly by PFT-parameter phenology.     \n
   !<                                                                                     \n
   !< 0: Predicted internally for all plant functional types.                             \n
   !< 1: Phenology is prescribed for cold-deciduous broadleaf trees and predicted for     \n
   !<    all other PFTs.                                                                  \n
   !---------------------------------------------------------------------------------------!
   integer                 :: iphen_scheme
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! REPRO_SCHEME -- This controls plant reproduction and dispersal.                       !
   !                 0.  Reproduction off.  Useful for very short runs only.               !
   !                 1.  Original reproduction scheme.  Seeds are exchanged between        !
   !                     patches belonging to the same site, but they can't go outside     !
   !                     their original site.                                              !
   !                 2.  Similar to 1, but seeds are exchanged between patches belonging   !
   !                     to the same polygon, even if they are in different sites.  They   !
   !                     can't go outside their original polygon, though.  This is the     !
   !                     same as option 1 if there is only one site per polygon.           !
   !                 3.  Similar to 2, but reproduction allocation for tropical trees is   !
   !                     assumed to be an asymptote as a function of height, following one !
   !                     of the functional forms proposed by Wenk and Falster (2015). This !
   !                     is experimental, use it at your own risk.                         !
   !---------------------------------------------------------------------------------------!
   integer                 :: repro_scheme
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! THETACRIT -- Leaf drought phenology threshold.  The sign matters here:                !
   !              >= 0. -- This is the relative soil moisture above the wilting point      !
   !                       below which the drought-deciduous plants will start shedding    !
   !                       their leaves                                                    !
   !              <  0. -- This is the soil potential in MPa below which the drought-      !
   !                       -deciduous plants will start shedding their leaves.  The wilt-  !
   !                       ing point is by definition -1.5MPa, so make sure that the value !
   !                       is above -1.5.                                                  !
   !---------------------------------------------------------------------------------------!
   real                   :: thetacrit
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     The following variables control the phenology prescribed from observations:       !
   !                                                                                       !
   ! IPHENYS1 -- First year for spring phenology                                           !
   ! IPHENYSF -- Final year for spring phenology                                           !
   ! IPHENYF1 -- First year for fall/autumn phenology                                      !
   ! IPHENYFF -- Final year for fall/autumn phenology                                      !
   ! PHENPATH -- path and prefix of the prescribed phenology data.                         !
   !                                                                                       !
   ! If the years don't cover the entire simulation period, they will be recycled.         !
   !---------------------------------------------------------------------------------------!
   integer                :: iphenys1
   integer                :: iphenysf
   integer                :: iphenyf1
   integer                :: iphenyff 
   character(len=str_len) :: phenpath
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Parameters that control the phenology response to radiation.                      !
   !                                                                                       !
   ! RADINT -- Intercept                                                                   !
   ! RADSLP -- Slope.                                                                      !
   !---------------------------------------------------------------------------------------!
   real                   :: radint
   real                   :: radslp
   !---------------------------------------------------------------------------------------!


   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !      Variables that are to be initialised in init_phen_coms (ed_params.f90).          !
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Before plants drop their leaves, they retain this fraction of their leaf carbon   !
   ! and nitrogen and put it into storage.                                                 !
   !---------------------------------------------------------------------------------------!
   real    :: retained_carbon_fraction
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Factor that controls the fine-root "elongation factor" relative to leaf           !
   ! elongation factor.  This is currently applied only for PFTs with phenology = 5.       !
   !                                                                                       !
   ! e_root = (e_leaf + root_phen_factor - 1) / root_phen_factor.                          !
   !                                                                                       !
   ! root_phen_factor > 1.  Fine roots will senesce more slowly than leaf shedding.        !
   ! root_phen_factor = 1.  Fine root elongation factor will be the same as for leaves.    !
   ! root_phen_factor < 1.  Fine roots will senesce more rapidly than leaf shedding.       !
   ! root_phen_factor = 0.  Special flag to disable fine-root phenology.                   !
   ! root_phen_factor < 0.  Non-sensical, currently assume the same as 0.                  !
   !---------------------------------------------------------------------------------------!
   real    :: root_phen_factor
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Flag that checks whether to Use soil potential rather than soil moisture to drive !
   ! phenology.                                                                            !
   !---------------------------------------------------------------------------------------!
   logical :: spot_phen
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Threshold for shedding all leaves when leaf water potential is very low. .       !
   !---------------------------------------------------------------------------------------!
   real    :: f_psi_xdry
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Minimum elongation factor before plants give up completely and shed all remain-  !
   ! ing leaves.                                                                           !
   !---------------------------------------------------------------------------------------!
   real    :: elongf_min
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Minimum elongation factor that allows plants to start flushing out new leaves if !
   ! they are drought deciduous and have been losing leaves.                               !
   !---------------------------------------------------------------------------------------!
   real    :: elongf_flush
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Leaf offset parameters are from:                                                 !
   !      White et al. 1997, Global Biogeochemical Cycles 11(2) 217-234                    !
   !---------------------------------------------------------------------------------------!
   real    :: dl_tr
   real    :: st_tr1
   real    :: st_tr2
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Phenology parameters for cold deciduous trees:                                   !
   !      Botta et al. 2000, Global Change Biology, 6, 709--725                            !
   !---------------------------------------------------------------------------------------!
   real    :: phen_a
   real    :: phen_b
   real    :: phen_c
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     This variable is the maximum distance between the coordinates of a prescribed     !
   ! phenology file and the actual polygon that we will still consider close enough to be  !
   ! representative.  If the user wants to run with prescribed phenology and the closest   !
   ! file is farther away from the polygon than the number below, the simulation will      !
   ! stop.                                                                                 !
   !---------------------------------------------------------------------------------------!
   real :: max_phenology_dist
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Derived type describing prescribed phenology.                                    !
   !---------------------------------------------------------------------------------------!
   type prescribed_phen
      !----- Number of years for which prescribed phenology is available. -----------------!
      integer :: nyears

      !----- The years for which prescribed phenology is available. -----------------------!
      integer, dimension(:), pointer :: years

      !----- Two parameters of the springtime logistic function describing leaf flush. ----!
      real, dimension(:), pointer :: flush_a
      real, dimension(:), pointer :: flush_b

      !----- Two parameters of the autumn logistic function describing leaf color. --------!
      real, dimension(:), pointer :: color_a
      real, dimension(:), pointer :: color_b
      
   end type prescribed_phen
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     These parameters control the specific leaf area as a function of the turnover     !
   ! rate                                                                                  !
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Variables controlling the light phenology as in Kim et al. (2012)                !
   !---------------------------------------------------------------------------------------!
   !----- Radiation window for running average [days] -------------------------------------!
   real :: radavg_window
   !----- Turnover weight, the inverse of the window. -------------------------------------!
   real :: radavg_wgt
   !----- Turnover window for running average [days] --------------------------------------!
   real :: turnamp_window
   !----- Turnover weight, the inverse of the window. -------------------------------------!
   real :: turnamp_wgt
   !----- Minimum instantaneous turnover rate amplitude [n/d]. ----------------------------!
   real :: turnamp_min
   !----- Maximum instantaneous turnover rate amplitude [n/d]. ----------------------------!
   real :: turnamp_max
   !----- Minimum radiation [W/m2], below which the turnover no longer responds. ----------!
   real :: radto_min
   !----- Maximum radiation [W/m2], above which the turnover no longer responds. ----------!
   real :: radto_max
   !----- Lifespan window for running average [days]. -------------------------------------!
   real :: llspan_window
   !----- Lifespan weight, the inverse of the window. -------------------------------------!
   real :: llspan_wgt
   !----- Minimum instantaneous life span [months]. ---------------------------------------!
   real :: llspan_min
   !----- Maximum instantaneous life span [months]. ---------------------------------------!
   real :: llspan_max
   !----- Instantaneous life span in case the turnover rate is 0. -------------------------!
   real :: llspan_inf
   !----- Vm0 window for running average [days]. ------------------------------------------!
   real :: vm0_window
   !----- Vm0 weight, the inverse of the window. ------------------------------------------!
   real :: vm0_wgt
   !----- Parameters that define the instantaneous Vm0 as a function of leaf life span. ---!
   real :: vm0_tran
   real :: vm0_slope
   real :: vm0_amp
   real :: vm0_min
   !----- SLA window for running average [days]. ------------------------------------------!
   real :: sla_window
   !----- SLA weight, the inverse of the window. ------------------------------------------!
   real :: sla_wgt
   !---------------------------------------------------------------------------------------!

end module phenology_coms
!==========================================================================================!
!==========================================================================================!
