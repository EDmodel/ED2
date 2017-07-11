!==========================================================================================!
!==========================================================================================!
!      This subroutine assigns the right number of vertical and "environmental" points     !
! associated with a given variable type.                                                   !
!------------------------------------------------------------------------------------------!
subroutine ze_dims(ng,idim,ismaster,nzpts,nepts)
   use mem_grid  , only : nnzp    & ! intent(in)
                        , npatch  & ! intent(in)
                        , nzg     & ! intent(in)
                        , nzs     ! ! intent(in)
   use node_mod  , only : mmzp    ! ! intent(in)
   use mem_aerad , only : nwave   ! ! intent(in)
   use mem_cuparm, only : nclouds ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer, intent(in)  :: ng       ! Grid number
   integer, intent(in)  :: idim     ! Dimension type
   logical, intent(in)  :: ismaster ! Flag to tell whether this is the master node or not.
   integer, intent(out) :: nzpts    ! Points in the vertical axis
   integer, intent(out) :: nepts    ! Points in the environmental "axis"
   !----- Local variables. ----------------------------------------------------------------!
   integer              :: nzatm    ! Number of vertical atmospheric levels.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------! 
   !    Check whether this is a master node or a computing node and pick the variable that !
   ! makes the most sense.  As of version 4.0.6, both variables should be the same, but    !
   ! that could change in the future.                                                      !
   !---------------------------------------------------------------------------------------! 
   if (ismaster) then
      nzatm = nnzp(ng)
   else
      nzatm = mmzp(ng)
   end if
   !---------------------------------------------------------------------------------------! 


   !----- Find the dimensions based on the variable type. ---------------------------------!
   select case(idim)
   case (2)
      !----- Two dimensions (longitude, latitude). ----------------------------------------!
      nzpts = 1
      nepts = 1

   case (3)
      !----- Three dimensions (altitude, longitude, latitude). ----------------------------!
      nzpts = nzatm
      nepts = 1
 
   case (4)
      !----- Four dimensions (soil depth, longitude, latitude, patch). --------------------!
      nzpts = nzg
      nepts = npatch
 
   case (5)
      !----- Four dimensions (snow depth, longitude, latitude, patch). --------------------!
      nzpts = nzs
      nepts = npatch
 
   case (6)
      !----- Three dimensions (longitude, latitude, patch). -------------------------------!
      nzpts = 1
      nepts = npatch
 
   case (7)
      !----- Three dimensions (longitude, latitude, wave lenght). -------------------------!
      nzpts = 1
      nepts = nwave
 
   case (8)
      !----- Four dimensions (altitude, longitude, latitude, cloud type). -----------------!
      nzpts = nzatm
      nepts = nclouds
 
   case (9)
      !----- Three dimensions (longitude, latitude, cloud type). --------------------------!
      nzpts = 1
      nepts = nclouds
 
   case default
      !----- Invalid dimension, kill the run and warn the user. ---------------------------!
      write (unit=*,fmt='(a)')       '----------------------------------------------------'
      write (unit=*,fmt='(a,1x,i6)') ' Dimension type : ',idim
      write (unit=*,fmt='(a)')       '----------------------------------------------------'
      call abort_run('Invalid dimension type!','ze_dims','varutils.f90')
   end select
   !---------------------------------------------------------------------------------------!  
   return
end subroutine ze_dims
!==========================================================================================!
!==========================================================================================!
