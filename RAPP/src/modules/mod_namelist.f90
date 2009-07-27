!==========================================================================================!
!==========================================================================================!
!  Module namelist: contains the namelist structure, as well as initial values.            !
!------------------------------------------------------------------------------------------!
module mod_namelist
   use mod_maxdims, only: maxstr
   implicit none

   type namelist_struct

      !----- Paths, prefixes, and types of inputs and outputs -----------------------------!
      character(len=maxstr) :: intype    ! Type of input file you want
      character(len=maxstr) :: inpath    ! Path of the input dataset
      character(len=maxstr) :: outpref   ! Prefix of output file
      
      !----- Time information -------------------------------------------------------------!
      integer               :: iyeara    ! First year to process
      integer               :: iyearz    ! Last year to process
      real                  :: inpfrq    ! Input variable frequency
      real                  :: radfrq    ! Output frequency for radiation fluxes
      real                  :: dtinc     ! Time increment for radiation interp. integration
      
      !----- Height information. ----------------------------------------------------------!
      real                  :: ref_hgt   ! Reference height.

      !----- Output edge information. -----------------------------------------------------!
      real                  :: lonw      ! Westernmost longitude
      real                  :: lone      ! Easternmost longitude
      real                  :: lats      ! Southernmost latitude
      real                  :: latn      ! Northernmost latitude
      integer               :: edgeoff   ! # of edge points not to be included in output.
      
      !----- Interportalion parameters (see Koch et al., 1983). ---------------------------!
      real                  :: gamma0    ! Parameter gamma
      real                  :: minweight ! Minimum weight to include the point in the OA.

   end type namelist_struct
   !---------------------------------------------------------------------------------------!

   type(namelist_struct) :: nl ! Structure containing all data read by the namelist


   contains
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine initialise_namelist()
      !------------------------------------------------------------------------------------!
      !   Subroutine that initialises the namelist structures with some default variables, !
      ! in case the user doesn't provide some structures. This is going to initialise all  !
      ! structures with non-sense numbers so it will fail in case the user doesn't fill    !
      ! the namelist properly and attempts to use it.                                      !
      !------------------------------------------------------------------------------------!
      implicit none
      
      !----- Initialising IO options. -----------------------------------------------------!
      nl%intype       = 'whatever'
      nl%inpath       = '/dev/null'
      nl%outpref      = '/dev/null'
      
      !----- Initialise time options. -----------------------------------------------------!
      nl%iyeara       = +huge(1)
      nl%iyearz       = -huge(1)
      nl%inpfrq       = -huge(1.) 
      nl%radfrq       = -huge(1.)
      nl%dtinc        = -huge(1.)
      
      !----- Initialise height options. ---------------------------------------------------!
      nl%ref_hgt      = -huge(1.)
      
      !----- Initialise the output domain options. ----------------------------------------!
      nl%lonw         = +huge(1.)
      nl%lone         = -huge(1.)
      nl%lats         = +huge(1.)
      nl%latn         = -huge(1.)
      nl%edgeoff      = -huge(1)
      
      !----- Initialise the interpolation parameters. -------------------------------------!
      nl%gamma0       = -huge(1.)
      nl%minweight    = -huge(1.)
      return
   end subroutine initialise_namelist
   !=======================================================================================!
   !=======================================================================================!
end module mod_namelist
!==========================================================================================!
!==========================================================================================!
