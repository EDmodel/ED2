!==========================================================================================!
!==========================================================================================!
!    Module detailed_coms: this module contains variables used to control some ED detailed !
!  output, which can be used for debugging.                                                !
!------------------------------------------------------------------------------------------!
module detailed_coms

   implicit none

   !---------------------------------------------------------------------------------------!
   !  IDETAILED -- This flag controls the possible detailed outputs, mostly used for       !
   !               debugging purposes.  Notice that this doesn't replace the normal debug- !
   !               ger options, the idea is to provide detailed output to check bad        !
   !               assumptions.  The options are additive, and the indices below represent !
   !               the different types of output:                                          !
   !                                                                                       !
   !               1  -- Detailed budget (every DTLSM)                                     !
   !               2  -- Detailed photosynthesis (every DTLSM)                             !
   !               4  -- Detailed output from the integrator (every HDID)                  !
   !               8  -- Thermodynamic bounds for sanity check (every DTLSM)               !
   !               16 -- Daily error stats (which variable caused the time step to shrink) !
   !               32 -- Allometry parameters, and minimum and maximum sizes               !
   !                     (two files, only at the beginning)                                !
   !                                                                                       !
   !               In case you don't want any detailed output (likely for most runs), set  !
   !               IDETAILED to zero.  In case you want to generate multiple outputs, add  !
   !               the number of the sought options: for example, if you want detailed     !
   !               photosynthesis and detailed output from the integrator, set IDETAILED   !
   !               to 6 (2 + 4).  Any combination of the above outputs is acceptable, al-  !
   !               though all but the last produce a sheer amount of txt files, in which   !
   !               case you may want to look at variable PATCH_KEEP.  It is also a good    !
   !               idea to set IVEGT_DYNAMICS to 0 when using the first five outputs.      !
   !---------------------------------------------------------------------------------------!
   integer :: idetailed
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! PATCH_KEEP -- This option will eliminate all patches except one from the initial-     !
   !               isation.  This is only used when one of the first five types of         !
   !               detailed output is active, otherwise it will be ignored.  Options are:  !
   !                 -2.  Keep only the patch with the lowest potential LAI                !
   !                 -1.  Keep only the patch with the highest potential LAI               !
   !                  0.  Keep all patches.                                                !
   !                > 0.  Keep the patch with the provided index.  In case the index is    !
   !                      not valid, the model will crash.                                 !
   !---------------------------------------------------------------------------------------!
   integer :: patch_keep
   !---------------------------------------------------------------------------------------!

end module detailed_coms
!==========================================================================================!
!==========================================================================================!
