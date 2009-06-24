!==========================================================================================!
!==========================================================================================!
!    Subroutine load_namelist. This subroutine will fill the appropriate modules with the  !
! information gathered by the namelist.                                                    !
!------------------------------------------------------------------------------------------!
subroutine load_namelist(rapp_in)
   use mod_maxdims,  only : maxstr               ! ! intent(in)
   use mod_namelist, only : nl                   & ! structure
                          , initialise_namelist  ! ! subroutine
   use mod_ioopts  , only : intype               & ! intent(out)
                          , inpath               & ! intent(out)
                          , outpref              & ! intent(out)
                          , iyeara               & ! intent(out)
                          , iyearz               & ! intent(out)
                          , inpfrq               & ! intent(out)
                          , radfrq               & ! intent(out)
                          , radratio             & ! intent(out)
                          , nsteps               & ! intent(out)
                          , nrads                & ! intent(out)
                          , nstepsi              & ! intent(out)
                          , nradsi               & ! intent(out)
                          , ref_hgt              & ! intent(out)
                          , lonw                 & ! intent(out)
                          , lone                 & ! intent(out)
                          , lats                 & ! intent(out)
                          , latn                 & ! intent(out)
                          , edgeoff              ! ! intent(out)
   use mod_interp  , only : dtinc                & ! intent(out)
                          , dtinc_fday           & ! intent(out)
                          , gamma0               & ! intent(out)
                          , minweight            & ! intent(out)
                          , minweight8           ! ! intent(out)
   use therm_lib   , only : level                & ! intent(out)
                          , vapour_on            & ! intent(out)
                          , cloud_on             & ! intent(out)
                          , bulk_on              ! ! intent(out)
   use rconstants  , only : day_sec              ! ! intent(in)
   implicit none
   !---------------------------------------------------------------------------------------!
   !    Variable declaration
   !---------------------------------------------------------------------------------------!
   !----- Arguments -----------------------------------------------------------------------!
   character(len=maxstr)              , intent(in) :: rapp_in   ! Namelist file name
   !----- Internal variables --------------------------------------------------------------!
   character(len=maxstr), dimension(1)             :: intypeaux ! Aux variable
   integer                                         :: ierr      ! Error flag
   !---------------------------------------------------------------------------------------!

   namelist /rapp_options/ nl

   !---------------------------------------------------------------------------------------!
   ! 1. First step is to initialise the namelist with some default (non-sense) parameters  !
   !---------------------------------------------------------------------------------------!
   call initialise_namelist()
   
   !---------------------------------------------------------------------------------------!
   ! 2. Reading the namelist                                                               !
   !---------------------------------------------------------------------------------------!
   open (unit=10,file=trim(rapp_in),status='old',action='read',iostat=ierr)
   if (ierr /= 0) then
      call fatal_error('Problems opening namelist '//trim(rapp_in)//'!!!'                  &
                      ,'load_namelist','load_namelist.f90')
   else
      read (unit=10,nml=rapp_options,iostat=ierr)
      if (ierr /= 0) then
        call fatal_error('Problems reading namelist '//trim(rapp_in)//'!!!'                &
                        ,'load_namelist','load_namelist.f90')
      end if
   end if

   !---------------------------------------------------------------------------------------!
   ! 3. Copying the namelist to the proper modules                                         !
   !---------------------------------------------------------------------------------------!
   intypeaux(1) = nl%intype
   inpath       = nl%inpath
   outpref      = nl%outpref

   iyeara       = nl%iyeara
   iyearz       = nl%iyearz
   inpfrq       = nl%inpfrq
   radfrq       = nl%radfrq
   dtinc        = nl%dtinc
   
   ref_hgt      = nl%ref_hgt

   lonw         = nl%lonw
   lone         = nl%lone
   lats         = nl%lats
   latn         = nl%latn
   
   edgeoff      = nl%edgeoff
   
   gamma0       = nl%gamma0
   minweight    = nl%minweight

   !---------------------------------------------------------------------------------------!
   ! 3. Converting some setting variables into lower-case, so the code is case insenstive  !
   !    for options given as characters.                                                   !
   !---------------------------------------------------------------------------------------!
   call tolower (intypeaux,1)
   intype = intypeaux(1)

   !---------------------------------------------------------------------------------------!
   ! 4. Copying the namelist to the proper modules                                         !
   !---------------------------------------------------------------------------------------!
   call dump_namelist()
   
   !---------------------------------------------------------------------------------------!
   ! 5. Checking whether the setting is a valid one.                                       !
   !---------------------------------------------------------------------------------------!
   call rapp_opspec()
   
   !---------------------------------------------------------------------------------------!
   ! 6. Find the various ratios between time scales                                        !
   !---------------------------------------------------------------------------------------!
   dtinc_fday  = dtinc/day_sec
   radratio    = nint(inpfrq / radfrq)
   nsteps      = nint(inpfrq / dtinc )
   nrads       = nint(radfrq / dtinc )
   nstepsi     =      dtinc  / inpfrq
   nradsi      =      dtinc  / radfrq

   !---------------------------------------------------------------------------------------!
   ! 7. Making the longitude in line with the meteorological dataset range (either 0..360  !
   !    or -180..180).                                                                     !
   !---------------------------------------------------------------------------------------!
   select case (trim(intype))
   case ('ncep')
      if (lonw < 0.) then
         lonw = lonw + 360.
         lone = lone + 360.
      end if
   end select

   !---------------------------------------------------------------------------------------!
   ! 8. Find the double precision version of minweight.                                    !
   !---------------------------------------------------------------------------------------!
   minweight8 = dble(minweight)


   !---------------------------------------------------------------------------------------!
   ! 9. Assigning the thermodynamic level, with full ice/liquid partition.                 !
   !---------------------------------------------------------------------------------------!
   level     = 3
   vapour_on = .true.
   cloud_on  = .true.
   bulk_on   = .true.

   return
end subroutine load_namelist
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine simply dumps the data read at the namelist on screen, so the users     !
! won't get bored and will able to check that RAPP understood what they asked.             !
!------------------------------------------------------------------------------------------!
subroutine dump_namelist()
   use mod_maxdims,  only : maxstr               ! ! intent(in)
   use mod_ioopts  , only : intype               & ! intent(in)
                          , inpath               & ! intent(in)
                          , outpref              & ! intent(in)
                          , iyeara               & ! intent(in)
                          , iyearz               & ! intent(in)
                          , inpfrq               & ! intent(in)
                          , radfrq               & ! intent(in)
                          , lonw                 & ! intent(in)
                          , lone                 & ! intent(in)
                          , lats                 & ! intent(in)
                          , latn                 & ! intent(in)
                          , edgeoff              ! ! intent(in)
   use mod_interp  , only : gamma0               & ! intent(in)
                          , minweight            ! ! intent(in)
   implicit none
   !----- Internal variables --------------------------------------------------------------!
   integer    :: v
   !---------------------------------------------------------------------------------------!

   write (unit=*,fmt='(102a)'      ) ('-',v=1,102)
   write (unit=*,fmt='(a)'         ) ' RAPP - Namelist information:                       '
   write (unit=*,fmt='(102a)'      ) ('-',v=1,102)
   write (unit=*,fmt='(a,1x,a)'    ) ' -> Intype:        ',trim(intype)
   write (unit=*,fmt='(a,1x,a)'    ) ' -> Inpath:        ',trim(inpath)
   write (unit=*,fmt='(a,1x,a)'    ) ' -> Outpref:       ',trim(outpref)
   write (unit=*,fmt='(a)'         ) ' '
   
   write (unit=*,fmt='(a,1x,i5)'   ) ' -> Iyeara:        ',iyeara
   write (unit=*,fmt='(a,1x,i5)'   ) ' -> Iyearz:        ',iyearz
   write (unit=*,fmt='(a,1x,f10.3)') ' -> Inpfrq:        ',inpfrq
   write (unit=*,fmt='(a,1x,f10.3)') ' -> Radfrq:        ',radfrq
   write (unit=*,fmt='(a)'         ) ' '

   write (unit=*,fmt='(a,1x,f10.3)') ' -> Lonw:          ',lonw
   write (unit=*,fmt='(a,1x,f10.3)') ' -> Lone:          ',lone
   write (unit=*,fmt='(a,1x,f10.3)') ' -> Lats:          ',lats
   write (unit=*,fmt='(a,1x,f10.3)') ' -> Latn:          ',latn
   write (unit=*,fmt='(a)'         ) ' '

   write (unit=*,fmt='(a,1x,i5)'   ) ' -> Edgeoff:       ',edgeoff
   write (unit=*,fmt='(a)'         ) ' '

   write (unit=*,fmt='(a,1x,f10.3)') ' -> Gamma0:        ',gamma0
   write (unit=*,fmt='(a,1x,f10.3)') ' -> Minweight:     ',minweight
   write (unit=*,fmt='(a)'         ) ' '

   write (unit=*,fmt='(102a)'      ) ('-',v=1,102)

   return
end subroutine dump_namelist
!==========================================================================================!
!==========================================================================================!
