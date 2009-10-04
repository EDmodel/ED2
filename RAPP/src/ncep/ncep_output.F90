!==========================================================================================!
!==========================================================================================!
!     This subroutine will control the output of the meteorological driver to an ED-       !
! friendly format.                                                                         !
!------------------------------------------------------------------------------------------!
subroutine ncep_output(month,year)
   use mod_maxdims, only : maxstr        ! ! intent(in)
   use mod_ioopts , only : outpref       & ! intent(in)
                         , edgeoff       ! ! intent(in)
   use mod_time   , only : monnames      ! ! intent(in)
   use mod_grid   , only : grid_g        & ! intent(in)
                         , sstp          & ! intent(in)
                         , ssxp          & ! intent(in)
                         , ssyp          ! ! intent(in)
   use mod_interp , only : ndownscal     ! ! intent(in)
   use mod_model  , only : ngrids        ! ! intent(in)
   use mod_ncep   , only : ncep_g        & ! intent(in)
                         , nvars_ol1     & ! intent(in)
                         , nvars_ol2     & ! intent(in)
                         , nvars_ol4     & ! intent(in)
                         , vars_ol1      & ! intent(in)
                         , vars_ol2      & ! intent(in)
                         , vars_ol4      ! ! intent(in)
#if USE_HDF5
   use hdf5
   use hdf5_utils !, only : shdf5_open_f  & ! function
                  !       , shdf5_orec_f  & ! function
                  !       , shdf5_close_f ! ! function
#endif
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                  , intent(in)  :: month
   integer                  , intent(in)  :: year
#if USE_HDF5
   !----- Local variables. ----------------------------------------------------------------!
   character(len=maxstr)                  :: ol1outname
   character(len=maxstr)                  :: ol2outname
   character(len=maxstr)                  :: olroutname
   character(len=3     )                  :: mmm
   integer                                :: ng
   integer                                :: nv
   integer                                :: nr
   integer                                :: mxo
   integer                                :: myo
   real   , dimension(:,:)  , allocatable :: outbuffg
   real   , dimension(:,:,:), allocatable :: outbuffs
   real   , dimension(:,:,:), allocatable :: outbufff
   integer, dimension(2)                  :: griddims
   integer, dimension(3)                  :: statedims
   integer, dimension(3)                  :: fluxdims
   !---------------------------------------------------------------------------------------!

   mmm=monnames(month)
   write (unit=*,fmt='(a,1x,2a,i4.4,a)')                                                   &
                                  '     - Writing the output files for',mmm,'/',year,'...'


   !---------------------------------------------------------------------------------------!
   !     Initialise the dimensions for output variables.  The output will contain the      !
   ! state arrays plus rain and LW radiation from grid one, and the other flux arrays from !
   ! grid two.  The points not further than edgeoff from the edge won't be included.       !
   !---------------------------------------------------------------------------------------!
   mxo = ssxp(1) - 2*edgeoff
   myo = ssyp(1) - 2*edgeoff

   !---------------------------------------------------------------------------------------!
   !     We now allocate the buffers.  We need these buffers because the output in HDF5    !
   ! for ED has time as the leftmost dimension, and also because the output array may be   !
   ! smaller than the working arrays.                                                      !
   !---------------------------------------------------------------------------------------!
   allocate(outbuffg(        mxo,myo))
   allocate(outbuffs(sstp(1),mxo,myo))
   allocate(outbufff(sstp(2),mxo,myo))

   !----- Assigning the dimensions for the dimenson vectors. ------------------------------!
   griddims  = (/              mxo,     myo /)
   statedims = (/ sstp(1),     mxo,     myo /)
   fluxdims  = (/ sstp(2),     mxo,     myo /)



   !----- Building the output file names. -------------------------------------------------!
   write(ol1outname,fmt='(2a,i4.4,2a)') trim(outpref),'_OL1_',year,trim(mmm),'.h5'
   write(ol2outname,fmt='(2a,i4.4,2a)') trim(outpref),'_OL2_',year,trim(mmm),'.h5'

   !---------------------------------------------------------------------------------------!
   !    We are now going to write the "Grid 1" file (the one with no time interpolation).  !
   ! We loop over the variables, just because if we decide to change the order, everything !
   ! will be done in one place (mod_ncep.f90).                                             !
   !---------------------------------------------------------------------------------------!
   write (unit=*,fmt='(a,1x,2a)') '         [|] Writing file:',trim(ol1outname),'...'
   call shdf5_open_f(trim(ol1outname),'W',1)
   ol1loop: do nv = 1, nvars_ol1
      select case (trim(vars_ol1(nv)))
      case ('lon')   !----- Longitude. ----------------------------------------------------!
         call fillbuff_2d(ssxp(1),ssyp(1),mxo,myo,edgeoff,grid_g(1)%lon,outbuffg)
         call shdf5_orec_f( 2, griddims ,trim(vars_ol1(nv)),rvara=outbuffg)

      case ('lat')   !----- Latitude. -----------------------------------------------------!
         call fillbuff_2d(ssxp(1),ssyp(1),mxo,myo,edgeoff,grid_g(1)%lat,outbuffg)
         call shdf5_orec_f( 2, griddims ,trim(vars_ol1(nv)),rvara=outbuffg)

      case ('hgt')   !----- Reference height. ---------------------------------------------!
         call fillbuff_2d(ssxp(1),ssyp(1),mxo,myo,edgeoff,grid_g(1)%lev,outbuffg)
         call shdf5_orec_f( 2, griddims ,trim(vars_ol1(nv)),rvara=outbuffg)

      case ('tmp')   !----- Temperature. --------------------------------------------------!
         call fillbuff_3d(ssxp(1),ssyp(1),sstp(1),mxo,myo,edgeoff,ncep_g(1)%temp,outbuffs)
         call shdf5_orec_f( 3, statedims,trim(vars_ol1(nv)),rvara=outbuffs)

      case ('pres')  !----- Pressure. -----------------------------------------------------!
         call fillbuff_3d(ssxp(1),ssyp(1),sstp(1),mxo,myo,edgeoff,ncep_g(1)%pres,outbuffs)
         call shdf5_orec_f( 3, statedims,trim(vars_ol1(nv)),rvara=outbuffs)

      case ('sh')    !----- Specific humidity. --------------------------------------------!
         call fillbuff_3d(ssxp(1),ssyp(1),sstp(1),mxo,myo,edgeoff,ncep_g(1)%shum,outbuffs)
         call shdf5_orec_f( 3, statedims,trim(vars_ol1(nv)),rvara=outbuffs)

      case ('ugrd')  !----- Zonal wind. ---------------------------------------------------!
         call fillbuff_3d(ssxp(1),ssyp(1),sstp(1),mxo,myo,edgeoff,ncep_g(1)%uwnd,outbuffs)
         call shdf5_orec_f( 3, statedims,trim(vars_ol1(nv)),rvara=outbuffs)

      case ('vgrd')  !----- Meridional wind. ----------------------------------------------!
         call fillbuff_3d(ssxp(1),ssyp(1),sstp(1),mxo,myo,edgeoff,ncep_g(1)%vwnd,outbuffs)
         call shdf5_orec_f( 3, statedims,trim(vars_ol1(nv)),rvara=outbuffs)

      case ('prate') !----- Precipitation rate --------------------------------------------!
         call fillbuff_3d(ssxp(1),ssyp(1),sstp(1),mxo,myo,edgeoff,ncep_g(1)%prate,outbuffs)
         call shdf5_orec_f( 3, statedims,trim(vars_ol1(nv)),rvara=outbuffs)

      case ('dlwrf') !----- Downwelling radiation flux ------------------------------------!
         call fillbuff_3d(ssxp(1),ssyp(1),sstp(1),mxo,myo,edgeoff,ncep_g(1)%dlwrf,outbuffs)
         call shdf5_orec_f( 3, statedims,trim(vars_ol1(nv)),rvara=outbuffs)

      end select
   end do ol1loop
   call shdf5_close_f()
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    We are now going to write the "Grid 2" file (the one with time interpolation).     !
   ! We loop over the variables, just because if we decide to change the order, everything !
   ! will be done in one place (mod_ncep.f90).                                             !
   !---------------------------------------------------------------------------------------!
   write (unit=*,fmt='(a,1x,2a)') '         [|] Writing file:',trim(ol2outname),'...'
   call shdf5_open_f(trim(ol2outname),'W',1)
   ol2loop: do nv = 1, nvars_ol2
      select case (trim(vars_ol2(nv)))
      case ('lon')   !----- Longitude. ----------------------------------------------------!
         call fillbuff_2d(ssxp(2),ssyp(2),mxo,myo,edgeoff,grid_g(2)%lon,outbuffg)
         call shdf5_orec_f( 2, griddims ,trim(vars_ol2(nv)),rvara=outbuffg)

      case ('lat')   !----- Latitude. -----------------------------------------------------!
         call fillbuff_2d(ssxp(2),ssyp(2),mxo,myo,edgeoff,grid_g(2)%lat,outbuffg)
         call shdf5_orec_f( 2, griddims ,trim(vars_ol2(nv)),rvara=outbuffg)

      case ('nbdsf') !----- Near-IR beam flux. --------------------------------------------!
         call fillbuff_3d(ssxp(2),ssyp(2),sstp(2),mxo,myo,edgeoff,ncep_g(2)%nbdsf,outbufff)
         call shdf5_orec_f( 3, fluxdims,trim(vars_ol2(nv)),rvara=outbufff)

      case ('nddsf') !----- Near-IR diffuse flux. -----------------------------------------!
         call fillbuff_3d(ssxp(2),ssyp(2),sstp(2),mxo,myo,edgeoff,ncep_g(2)%nddsf,outbufff)
         call shdf5_orec_f( 3, fluxdims,trim(vars_ol2(nv)),rvara=outbufff)

      case ('vbdsf') !----- Visible beam flux. --------------------------------------------!
         call fillbuff_3d(ssxp(2),ssyp(2),sstp(2),mxo,myo,edgeoff,ncep_g(2)%vbdsf,outbufff)
         call shdf5_orec_f( 3, fluxdims,trim(vars_ol2(nv)),rvara=outbufff)

      case ('vddsf') !----- Visible diffuse flux. -----------------------------------------!
         call fillbuff_3d(ssxp(2),ssyp(2),sstp(2),mxo,myo,edgeoff,ncep_g(2)%vddsf,outbufff)
         call shdf5_orec_f( 3, fluxdims,trim(vars_ol2(nv)),rvara=outbufff)

      end select
   end do ol2loop
   call shdf5_close_f()
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    We will now write the files with the different precipitation realisations.         !
   !---------------------------------------------------------------------------------------!
   dnscloop: do nr=1,ndownscal
      write(olroutname,fmt='(2a,i2.2,a,i4.4,2a)')                                          &
                                             trim(outpref),'_R',nr,'_',year,trim(mmm),'.h5'

      ng = nr + 3

      !------------------------------------------------------------------------------------!
      !    We are now going to write the realisation files (with time interpolation).      !
      ! Each downscaled precipitation will go to a different file.  As in the previous     !
      ! cases, we loop over the variables, just because if we decide to change the order,  !
      ! everything will be done in one place (mod_ncep.f90).                               !
      !------------------------------------------------------------------------------------!
      write (unit=*,fmt='(a,1x,2a)') '         [|] Writing file:',trim(olroutname),'...'
      call shdf5_open_f(trim(olroutname),'W',1)
      ol4loop: do nv = 1, nvars_ol4
         select case (trim(vars_ol4(nv)))
         case ('lon')   !----- Longitude. -------------------------------------------------!
            call fillbuff_2d(ssxp(ng),ssyp(ng),mxo,myo,edgeoff,grid_g(ng)%lon,outbuffg)
            call shdf5_orec_f( 2, griddims ,trim(vars_ol4(nv)),rvara=outbuffg)

         case ('lat')   !----- Latitude. --------------------------------------------------!
            call fillbuff_2d(ssxp(ng),ssyp(ng),mxo,myo,edgeoff,grid_g(ng)%lat,outbuffg)
            call shdf5_orec_f( 2, griddims ,trim(vars_ol4(nv)),rvara=outbuffg)

         case ('prate') !----- Precipitation rate -----------------------------------------!
            call fillbuff_3d(ssxp(ng),ssyp(ng),sstp(ng),mxo,myo,edgeoff,ncep_g(ng)%prate   &
                            ,outbufff)
            call shdf5_orec_f( 3, fluxdims,trim(vars_ol4(nv)),rvara=outbufff)

         end select
      end do ol4loop
      call shdf5_close_f()
      !------------------------------------------------------------------------------------!
   end do dnscloop



   deallocate(outbuffg)
   deallocate(outbuffs)
   deallocate(outbufff)
#else
   call fatal_error('You can''t use hdf5 routines without compiling RAPP with HDF5.'       &
                         ,'ncep_output','ncep_output.F90')

#endif
   return
end subroutine ncep_output
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will copy the subset of work array ARIN to the output buffer AROUT,  !
! for time-invariant variables.                                                            !
!------------------------------------------------------------------------------------------!
subroutine fillbuff_2d(mxi,myi,mxo,myo,eoff,arin,arout)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                        , intent(in)  :: mxi   ! # of X points, input array
   integer                        , intent(in)  :: myi   ! # of Y points, input array
   integer                        , intent(in)  :: mxo   ! # of X points, output array
   integer                        , intent(in)  :: myo   ! # of Y points, output array
   integer                        , intent(in)  :: eoff  ! Edge offset of output array
   real   , dimension(mxi,myi)    , intent(in)  :: arin  ! Input array  (work)
   real   , dimension(mxo,myo)    , intent(out) :: arout ! Output array (buffer)
   !----- Local variables. ----------------------------------------------------------------!
   integer                                      :: xi    ! Index of input array, X
   integer                                      :: yi    ! Index of input array, Y
   integer                                      :: xo    ! Index of output array, X
   integer                                      :: yo    ! Index of output array, Y
   logical                                      :: okx   ! Flag for correct X dimensions
   logical                                      :: oky   ! Flag for correct Y dimensions
   !---------------------------------------------------------------------------------------!

   !----- Quick dimension check. ----------------------------------------------------------!
   okx = mxo == mxi - 2 * eoff
   oky = myo == myi - 2 * eoff

   !----- Run only when the input vs. output dimensions are consistent. -------------------!
   if (okx .and. oky) then
      yoloop: do yo=1,myo
         yi = yo + eoff

         xoloop: do xo=1,mxo
            xi = xo + eoff
            arout(xo,yo) = arin(xi,yi) 
         end do xoloop
      end do yoloop
   else
      write (unit=*,fmt='(a,1x,i5)') 'MXI  =',mxi
      write (unit=*,fmt='(a,1x,i5)') 'MYI  =',myi
      write (unit=*,fmt='(a,1x,i5)') 'MXO  =',mxo
      write (unit=*,fmt='(a,1x,i5)') 'MYO  =',myo
      write (unit=*,fmt='(a,1x,i5)') 'EOFF =',eoff
      write (unit=*,fmt='(a,1x,l1)') 'MXO = MXI-2EOFF?',okx
      write (unit=*,fmt='(a,1x,l1)') 'MYO = MYI-2EOFF?',oky
      call fatal_error('Mismatch between input and output grid (check values above...)'    &
                      ,'fillbuff_2d','ncep_output.F90')
   end if

   return
end subroutine fillbuff_2d
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will copy the subset of work array ARIN to the output buffer AROUT,  !
! for time-dependent variables.  The output buffer has the time index as its leftmost      !
! dimension.                                                                               !
!------------------------------------------------------------------------------------------!
subroutine fillbuff_3d(mxi,myi,mtp,mxo,myo,eoff,arin,arout)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                        , intent(in)  :: mxi   ! # of X points, input array
   integer                        , intent(in)  :: myi   ! # of Y points, input array
   integer                        , intent(in)  :: mtp   ! # of T points, both arrays
   integer                        , intent(in)  :: mxo   ! # of X points, output array
   integer                        , intent(in)  :: myo   ! # of Y points, output array
   integer                        , intent(in)  :: eoff  ! Edge offset of output array
   real   , dimension(mxi,myi,mtp), intent(in)  :: arin  ! Input array  (work)
   real   , dimension(mtp,mxo,myo), intent(out) :: arout ! Output array (buffer)
   !----- Local variables. ----------------------------------------------------------------!
   integer                                      :: xi    ! Index of input array, X
   integer                                      :: yi    ! Index of input array, Y
   integer                                      :: tt    ! Time index, both arrays
   integer                                      :: xo    ! Index of output array, X
   integer                                      :: yo    ! Index of output array, Y
   logical                                      :: okx   ! Flag for correct X dimensions
   logical                                      :: oky   ! Flag for correct Y dimensions
   !---------------------------------------------------------------------------------------!

   !----- Quick dimension check. ----------------------------------------------------------!
   okx = mxo == mxi - 2 * eoff
   oky = myo == myi - 2 * eoff

   !----- Run only when the input vs. output dimensions are consistent. -------------------!
   if (okx .and. oky) then
      yoloop: do yo=1,myo
         yi = yo + eoff

         xoloop: do xo=1,mxo
            xi = xo + eoff

            ttloop: do tt=1,mtp
               arout(tt,xo,yo) = arin(xi,yi,tt) 
            end do ttloop

         end do xoloop
      end do yoloop
   else
      write (unit=*,fmt='(a,1x,i5)') 'MXI  =',mxi
      write (unit=*,fmt='(a,1x,i5)') 'MYI  =',myi
      write (unit=*,fmt='(a,1x,i5)') 'MXO  =',mxo
      write (unit=*,fmt='(a,1x,i5)') 'MYO  =',myo
      write (unit=*,fmt='(a,1x,i5)') 'EOFF =',eoff
      write (unit=*,fmt='(a,1x,l1)') 'MXO = MXI-2EOFF?',okx
      write (unit=*,fmt='(a,1x,l1)') 'MYO = MYI-2EOFF?',oky
      call fatal_error('Mismatch between input and output grid (check values above...)'    &
                      ,'fillbuff_3d','ncep_output.F90')
   end if

   return
end subroutine fillbuff_3d
!==========================================================================================!
!==========================================================================================!
