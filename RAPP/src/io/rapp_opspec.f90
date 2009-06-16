!==========================================================================================!
!==========================================================================================!
!    Subroutine rapp_opspec.  This subroutine will check the options provided by the user, !
! to make the run to crash in case some set doesn´t make sense. It is called right after   !
! the namelist is loaded.                                                                  !
!------------------------------------------------------------------------------------------!
subroutine rapp_opspec()
   use mod_maxdims , only : maxstr               ! ! intent(in)
   use mod_ioopts  , only : intype               & ! intent(in)
                          , iyeara               & ! intent(in)
                          , iyearz               & ! intent(in)
                          , inpfrq               & ! intent(in)
                          , radfrq               & ! intent(in)
                          , lonw                 & ! intent(in)
                          , lone                 & ! intent(in)
                          , lats                 & ! intent(in)
                          , latn                 ! ! intent(in)
   use rconstants  , only : day_sec              ! ! intent(in)
   implicit none
   !----- Internal variables --------------------------------------------------------------!
   character(len=maxstr)             :: reason
   integer                           :: ifaterr ! Error flag
   !---------------------------------------------------------------------------------------!
   
   !----- ifaterr will count how many bad settings the user made. -------------------------!
   ifaterr = 0
   
   !----- Checking the type of input data. ------------------------------------------------!
   select case (trim(intype))
   case ('ncep')
      !----- Valid type, move on... -------------------------------------------------------!
      continue
   case default
      write (reason,fmt='(a,1x,a,a)')                                                      &
        'Invalid INTYPE, the only valid type is NCEP and yours is set to',intype,'...'
      call opspec_fatal(reason,'rapp_opspec')
      ifaterr = ifaterr +1
   end select
   !---------------------------------------------------------------------------------------! 



   
   !----- Check whether iyeara and iyearz make sense... -----------------------------------!
   if (iyearz < iyeara) then
      write (reason,fmt='(a,2(1x,i4,1x,a))')                                               &
         'IYEARZ must be after IYEARA, or at least the same. Your IYEARA =',iyeara         &
        ,'and your IYEARZ =',iyearz,'...'
      call opspec_fatal(reason,'rapp_opspec')
      ifaterr = ifaterr +1   
   end if
   !---------------------------------------------------------------------------------------! 




   !----- Check whether the input frequency makes sense. It must be a divisor of one day. -!
   if (inpfrq <= 0.) then
      write(reason,fmt='(a,1x,es14.7,a)')                                                  &
         'INPFRQ must be positive, and yours is set to ',inpfrq,'...'
      call opspec_fatal(reason,'rapp_opspec')
      ifaterr = ifaterr +1   
   elseif (mod(day_sec,inpfrq) /= 0.) then
      write(reason,fmt='(a,1x,es14.7,a)')                                                  &
         'INPFRQ must be a divisor of one day, and yours is set to ',inpfrq,' seconds...'
      call opspec_fatal(reason,'rapp_opspec')
      ifaterr = ifaterr +1   
   end if
   !---------------------------------------------------------------------------------------! 




   !---------------------------------------------------------------------------------------! 
   !    Check whether the output frequency for radiation makes sense. It must be a divisor !
   ! of the input frequency.                                                               !
   !---------------------------------------------------------------------------------------!
   if (radfrq <= 0.) then
      write(reason,fmt='(a,1x,es14.7,a)')                                                  &
         'RADFRQ must be positive, and yours is set to ',radfrq,'...'
      call opspec_fatal(reason,'rapp_opspec')
      ifaterr = ifaterr +1   
   elseif (mod(day_sec,inpfrq) /= 0.) then
      write(reason,fmt='(a,1x,es14.7,a)')                                                  &
          'RADFRQ must be a divisor of INPFRQ. Currently RADFRQ =',radfrq,' and '          &
         ,'INPFRQ =',inpfrq,'...'
      call opspec_fatal(reason,'rapp_opspec')
      ifaterr = ifaterr +1   
   end if
   !---------------------------------------------------------------------------------------! 



   !---------------------------------------------------------------------------------------! 
   !     Check whether the longitude and latidute make sense.                              !
   !---------------------------------------------------------------------------------------! 
   if (lonw < -180. .or. lonw > 180.) then
      write (reason,fmt='(a,1x,f9.3,a)')                                                   &
          'Invalid LONW (between -180 and 180).  Currently LONW is set to ',lonw,'...'
   end if
   if (lone < -180. .or. lone > 180.) then
      write (reason,fmt='(a,1x,f9.3,a)')                                                   &
          'Invalid LONE (between -180 and 180).  Currently LONE is set to ',lonw,'...'
   end if
   if (lonw >= lone) then
      write (reason,fmt='(a,2(1x,f9.3,1x,a))')                                             &
          'LONE must be greater than LONW.  Currently LONW =',lonw,'and LONE =',lone,'...'
   end if
   if (lats < -90. .or. lats > 90.) then
      write (reason,fmt='(a,1x,f9.3,a)')                                                   &
          'Invalid LATS (between -90 and 90).  Currently LATS is set to ',lats,'...'
   end if
   if (latn < -90. .or. latn > 90.) then
      write (reason,fmt='(a,1x,f9.3,a)')                                                   &
          'Invalid LATN (between -90 and 90).  Currently LATN is set to ',latn,'...'
   end if
   if (lats >= latn) then
      write (reason,fmt='(a,2(1x,f9.3,1x,a))')                                             &
          'LATN must be greater than LATS.  Currently LATS =',lats,'and LATN =',latn,'...'
   end if


   !----- Stop the run if there are any fatal errors. -------------------------------------!
   if (ifaterr > 0) then
      write (unit=*,fmt='(a)')       ' -----------RAPP_OPSPEC -----------------------------'
      write (unit=*,fmt='(a,1x,i5)') ' Fatal errors:',ifaterr
      write (unit=*,fmt='(a)')       ' ----------------------------------------------------'
      call fatal_error('Fatal errors at namelist','rapp_opspec','rapp_opspec.f90')
   end if
   
   return
end subroutine rapp_opspec
!==========================================================================================!
!==========================================================================================!
