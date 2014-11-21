!==========================================================================================!
!==========================================================================================!
!                                                                                          !
!                         RAMSPOST - RAMS Post Processor for GrADS.                        !
!                                                                                          !
!------------------------------------------------------------------------------------------!
program ramspost
   use rpost_coms
   use rout_coms , only  : rout               & ! intent(out)
                         , routgrads          & ! intent(out)
                         , alloc_rout         & ! subroutine
                         , undef_rout         & ! subroutine
                         , undefflg           ! ! intent(in)
   use rpost_dims, only  : str_len            ! ! intent(in)
   use brams_data
   use leaf_coms , only  : sfclyr_init_params ! ! sub-routine
   use somevars
   implicit none
   !---- Fixed size variables. ------------------------------------------------------------!
   character(len=str_len), dimension(maxfiles)                       :: fln
   character(len=str_len)                                            :: inp
   character(len=str_len)                                            :: fprefix
   character(len=str_len)                                            :: gprefix
   character(len=str_len)                                            :: cfln
   character(len=str_len), dimension(maxfiles)                       :: wfln
   character(len=40)     , dimension(maxvars)                        :: vpln
   character(len=40)                                                 :: tmpdesc
   character(len=20)     , dimension(maxvars)                        :: vp
   character(len=10)     , dimension(maxvars)                        :: vpun
   character(len=20)                                                 :: tmpvar
   character(len=20)                                                 :: proj
   character(len=1)                                                  :: cgrid
   character(len=2)                                                  :: ccgrid
   character(len=2)                                                  :: patchnumber
   character(len=2)                                                  :: cldnumber
   character(len=15)                                                 :: chdate
   character(len=15)                                                 :: chstep
   character(len=15)                                                 :: xchstep
   character(len=20)                                                 :: cvar
   integer                                                           :: nvp
   integer                                                           :: nfiles
   integer               , dimension(maxvars)                        :: nzvp
   integer                                                           :: nrec
   integer                                                           :: ipresslev
   real                  , dimension(nplmax)                         :: iplevs
   integer                                                           :: inplevs
   integer               , dimension(maxgrds)                        :: zlevmax
   integer               , dimension(maxvars)                        :: ndim
   integer                                                           :: fim_inp
   integer                                                           :: hunit
   integer                                                           :: nstep
   integer                                                           :: ip
   integer                                                           :: ic
   integer                                                           :: iep_ng
   integer                                                           :: iep_np
   integer                                                           :: iep_nc
   integer                                                           :: iep_ngrids
   integer                                                           :: iyear
   integer                                                           :: imonth
   integer                                                           :: idate
   integer                                                           :: ihour
   integer                                                           :: imin
   integer                                                           :: iistep
   integer                                                           :: ng
   integer                                                           :: nnvp
   integer                                                           :: iv
   integer                                                           :: nfn
   integer                                                           :: n
   integer                                                           :: imon
   integer                                                           :: nxpg
   integer                                                           :: nypg
   integer                                                           :: nzpg
   integer               , dimension(maxgrds)                        :: nxgrads
   integer               , dimension(maxgrds)                        :: nygrads
   integer               , dimension(maxgrds)                        :: nxa
   integer               , dimension(maxgrds)                        :: nxb
   integer               , dimension(maxgrds)                        :: nya
   integer               , dimension(maxgrds)                        :: nyb
   integer               , dimension(maxgrds)                        :: iep_nx
   integer               , dimension(maxgrds)                        :: iep_ny
   integer               , dimension(maxgrds)                        :: iep_nz
   real                  , dimension(maxgrds)                        :: lati
   real                  , dimension(maxgrds)                        :: latf
   real                  , dimension(maxgrds)                        :: loni
   real                  , dimension(maxgrds)                        :: lonf
   real                  , dimension(nzpmax,maxgrds)                 :: dep_zlev
   real                  , dimension(2,maxgrds)                      :: dep_glat
   real                  , dimension(2,maxgrds)                      :: dep_glon
   real                                                              :: rlatmin
   real                                                              :: rlatmax
   real                                                              :: rlonmin
   real                                                              :: rlonmax
   !----- Constants. ----------------------------------------------------------------------!
   character(len=3), dimension(12), parameter :: cmo = (/ 'jan', 'feb', 'mar', 'apr'       &
                                                        , 'may', 'jun', 'jul', 'aug'       &
                                                        , 'sep', 'oct', 'nov', 'dec' /)
   character(len=3), dimension(3) , parameter :: ctu = (/ 'sec',  'mn',  'hr' /)
   !---------------------------------------------------------------------------------------!


   namelist /rp_input/  fprefix,nvp,vp,gprefix,nstep,proj,lati,latf,loni,lonf,zlevmax      &
                       ,ipresslev,inplevs,iplevs


   !---------------------------------------------------------------------------------------!
   !    Print a banner to entertain the user.                                              !
   !---------------------------------------------------------------------------------------!
   write (unit=*,fmt='(a)') '=========================================================='
   write (unit=*,fmt='(a)') '    Ramspost for BRAMS-4.0.6 '
   write (unit=*,fmt='(a)') '=========================================================='


   !---------------------------------------------------------------------------------------!
   !    Load the namelist filename.  If none is given, assume default.                     !
   !---------------------------------------------------------------------------------------!
   call getarg(1,inp)
   fim_inp=len_trim(inp)
   if (fim_inp < 3) inp='ramspost.inp'
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Read the namelist.                                                                 !
   !---------------------------------------------------------------------------------------!
   write (unit=*,fmt='(a)' )      ' '
   write (unit=*,fmt='(3(a,1x))') ' - Readin namelist from file ',trim(inp),'...'

   open  (unit=15,file=trim(inp),status='old',action='read')
   read  (unit=15,nml=rp_input)
   close (unit=15,status='keep')
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Make some character variables lower case.                                          !
   !---------------------------------------------------------------------------------------!
   call tolower_sca(proj)
   !---------------------------------------------------------------------------------------!



   !----- Determine some dimensions. ------------------------------------------------------!
   call RAMS_anal_init(nfiles,fln,fprefix,dep_zlev,iep_nx,iep_ny,iep_nz,iep_ng,iep_np      &
                      ,iep_nc,iep_ngrids)
   write (unit=*,fmt='(92a)'    )    ('-',n=1,92)
   write (unit=*,fmt='(a,1x,i5)')    ' + NFILES     =',nfiles
   write (unit=*,fmt='(a,1x,i5)')    ' + IEP_NG     =',iep_ng
   write (unit=*,fmt='(a,1x,i5)')    ' + IEP_NP     =',iep_np
   write (unit=*,fmt='(a,1x,i5)')    ' + IEP_NC     =',iep_nc
   write (unit=*,fmt='(a,1x,i5)')    ' + IEP_NGRIDS =',iep_ngrids
   do ng=1,iep_ngrids
      !------------------------------------------------------------------------------------!
      !     Make sure that we don't exceed the maximum number of layers. allowed in this   !
      ! grid.                                                                              !
      !------------------------------------------------------------------------------------!
      if (zlevmax(ng) >= iep_nz(ng)) zlevmax(ng) = iep_nz(ng)-1
      !------------------------------------------------------------------------------------!

      write (unit=*,fmt='(a,1x,i5)') ' + GRID       =',ng
      write (unit=*,fmt='(a,1x,i5)') '   - IEP_NX   =',iep_nx(ng)
      write (unit=*,fmt='(a,1x,i5)') '   - IEP_NY   =',iep_ny(ng)
      write (unit=*,fmt='(a,1x,i5)') '   - IEP_NZ   =',iep_nz(ng)
      write (unit=*,fmt='(a,1x,i5)') '   - ZLEVMAX  =',zlevmax(ng)
   end do
   write (unit=*,fmt='(92a)'    )    ('-',n=1,92)
   write (unit=*,fmt='(a)'      )    ' '
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Determine the initial time, the time step, and the best units for time interval.  !
   !---------------------------------------------------------------------------------------!
   call RAMS_get_time_init(1,iyear,imonth,idate,ihour,imin)
   call RAMS_get_time_step(iistep,hunit,nfiles)
   write(chdate,fmt='(3(i2.2,a),i4.4)') ihour,':',imin,'z',idate,cmo(imonth),iyear
   write(chstep,fmt='(8x,i3,a)') iistep,trim(ctu(hunit))
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Allocate the output buffer structures.                                           !
   !---------------------------------------------------------------------------------------!
   allocate (rout     (iep_ngrids))
   allocate (routgrads(iep_ngrids))
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Grid loop.                                                                       !
   !---------------------------------------------------------------------------------------!
   gridloop: do ng=1,iep_ngrids
      write (unit=*,fmt='(92a)'    )    ('=',n=1,92)
      write (unit=*,fmt='(a)'      )    ' '
      write (unit=*,fmt='(a,1x,i5)')    ' + Writing Grid ',ng


      !------------------------------------------------------------------------------------!
      !    Allocate pointers from rout and routgrads.                                      !
      !------------------------------------------------------------------------------------!
      call alloc_rout(rout(ng),iep_nx(ng),iep_ny(ng),iep_nz(ng),iep_ng,iep_np,iep_nc       &
                     ,inplevs)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Initialise NNVP with the number of variables coming from the namelist.  In case !
      ! it is a multiple-class variable, the actual number of output variables will        !
      ! increase.  Conversely, if the user asked for a variable that doesn't exist, we     !
      ! take one number out and warn him/her.                                              !
      !------------------------------------------------------------------------------------!
      nnvp = nvp
      !------------------------------------------------------------------------------------!

      write(cgrid,'(i1)') ng

      iv   = 1

      !----- Get the prefix and remove the trailing -head.txt so we can append grid info. -!
      cfln = trim(fln(1))
      ip   = len_trim(cfln) - 9
      cfln = cfln(1:ip)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Read latitude and longitude of the "thermodynamic points" of BRAMS.  Save them  !
      ! to rlat and rlon, respectively.                                                    !
      !------------------------------------------------------------------------------------!
      write(unit=*,fmt='(a)')       ' '
      write(unit=*,fmt='(2(a,1x))') '   - Variable:  ','lat'

      call ep_getvar('lat',iep_nx(ng),iep_ny(ng),1,ng,cfln,vpln(iv),vpun(iv),n,iep_np      &
                    ,iep_nc,iep_ng)
      call atob(iep_nx(ng)*iep_ny(ng),rout(ng)%r2,rout(ng)%rlat)
      write(unit=*,fmt='(a,1x,i5)') '     # Output variable type:  ',n

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(2(a,1x))') '   - Variable:  ','lon'

      call ep_getvar('lon',iep_nx(ng),iep_ny(ng),1,ng,cfln ,vpln(iv),vpun(iv),n,iep_np     &
                    ,iep_nc,iep_ng)
      call atob(iep_nx(ng)*iep_ny(ng),rout(ng)%r2,rout(ng)%rlon)
      write(unit=*,fmt='(a,1x,i5)') '     # Output variable type:  ',n
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the dimensions of x and y domains depending on the sought projection, and !
      ! find the mapping to interpolate values to the regular lon-lat grid if needed.      !
      !------------------------------------------------------------------------------------!
      call geo_grid(iep_nx(ng),iep_ny(ng),rout(ng)%rlat,rout(ng)%rlon,dep_glon(1,ng)       &
                   ,dep_glon(2,ng),dep_glat(1,ng),dep_glat(2,ng),rlatmin,rlatmax,rlonmin   &
                   ,rlonmax,nxgrads(ng),nygrads(ng),proj)
      call alloc_rout(routgrads(ng),nxgrads(ng),nygrads(ng),iep_nz(ng),iep_ng,iep_np       &
                     ,iep_nc,inplevs)
      call array_interpol(ng,nxgrads(ng),nygrads(ng),iep_nx(ng),iep_ny(ng)                 &
                         ,routgrads(ng)%iinf,routgrads(ng)%jinf,routgrads(ng)%rmi,proj)
      call define_lim(ng,nxgrads(ng),nygrads(ng),dep_glat(1,ng),dep_glat(2,ng)             &
                     ,dep_glon(1,ng),dep_glon(2,ng),lati(ng),latf(ng),loni(ng),lonf(ng)    &
                     ,nxa(ng),nxb(ng),nya(ng),nyb(ng),proj,iep_nx(ng),iep_ny(ng)           &
                     ,rout(ng)%rlat,rout(ng)%rlon)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Open the binary file now.                                                     !
      !------------------------------------------------------------------------------------!
      open(unit=19,file=trim(gprefix)//'_g'//cgrid//'.gra',form='unformatted'              &
          ,access='direct',status='replace',action='write'                                 &
          ,recl=4*(nxb(ng)-nxa(ng)+1)*(nyb(ng)-nya(ng)+1))
      nrec = 0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Loop over all files that are to be used.                                        !
      !------------------------------------------------------------------------------------!
      fileloop: do nfn=1,nfiles,nstep
         write(unit=*,fmt='(a)')        ' '
         write(unit=*,fmt='(a,1x,i5)')  '   - Timestep: ',nfn

         !---------------------------------------------------------------------------------!
         !      Get the prefix and remove the trailing -head.txt so we can append grid     !
         ! info.                                                                           !
         !---------------------------------------------------------------------------------!
         cfln = trim(fln(nfn))
         ip   = len_trim(cfln) - 9
         cfln = cfln(1:ip)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     In case ipresslev is not zero (vertical intepolation), we must load         !
         ! topography and Exner function so we can interpolate variables.                  !
         !---------------------------------------------------------------------------------!
         select case (ipresslev)
         case (0)
            continue
         case default
            !----- Load topography. -------------------------------------------------------!
            write(unit=*,fmt='(a)') ' '
            write(unit=*,fmt='(2(a,1x))') '     * Variable:  ','topo'
            call ep_getvar('topo',iep_nx(ng),iep_ny(ng),iep_nz(ng),ng,cfln,vpln(iv),vpun(iv),n      &
                          ,iep_np,iep_nc,iep_ng)
            call atob(iep_nx(ng)*iep_ny(ng),rout(ng)%r2,rout(ng)%topo)
            write(unit=*,fmt='(a,1x,i5)') '       # Output variable type:  ',n
            !------------------------------------------------------------------------------!

            !----- Load Exner function. ---------------------------------------------------!
            write(unit=*,fmt='(a)') ' '
            write(unit=*,fmt='(2(a,1x))') '     * Variable:  ','pi'
            call ep_getvar('pi',iep_nx(ng),iep_ny(ng),iep_nz(ng),ng,cfln,vpln(iv),vpun(iv) &
                          ,n,iep_np,iep_nc,iep_ng)
            write(unit=*,fmt='(a)')       '       # Got variable.  Entering atob...'
            call atob(iep_nx(ng)*iep_ny(ng)*iep_nz(ng),rout(ng)%r3,rout(ng)%exner)
            write(unit=*,fmt='(a,1x,i5)') '       # Output variable type:  ',n
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Loop over all other variables.                                              !
         !---------------------------------------------------------------------------------!
         varloop: do iv=1,nvp
            call undef_rout(rout(ng)     ,.false.)
            call undef_rout(routgrads(ng),.false.)
            write(unit=*,fmt='(a)') ' '
            write(unit=*,fmt='(2(a,1x))') '     * Variable:  ',trim(vp(iv))
            call ep_getvar(trim(vp(iv)),iep_nx(ng),iep_ny(ng),iep_nz(ng),ng,cfln,vpln(iv)  &
                          ,vpun(iv),ndim(iv),iep_np,iep_nc,iep_ng)
            write(unit=*,fmt='(a,1x,i5)') '       # Output variable type:  ',ndim(iv)


            !------------------------------------------------------------------------------!
            !     Decide how to output variable depending on the variable type.            !
            !------------------------------------------------------------------------------!
            select case (ndim(iv))
            case (2)
               !---------------------------------------------------------------------------!
               !    Two-dimensional array, no vertical information.                        !
               !---------------------------------------------------------------------------!


               !----- Set the number of levels. -------------------------------------------!
               nzvp(iv) = 1
               cvar     = trim(vp(iv))
               !---------------------------------------------------------------------------!


               !----- Adjust projection to GrADS. -----------------------------------------!
               call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng),nzvp(iv)      &
                                      ,nxgrads(ng),nygrads(ng),routgrads(ng)%rmi           &
                                      ,routgrads(ng)%iinf,routgrads(ng)%jinf,rout(ng)%r2   &
                                      ,routgrads(ng)%r2,proj)
               !---------------------------------------------------------------------------!


               !----- Dump array to output file. ------------------------------------------!
               call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng),nya(ng)  &
                             ,nyb(ng),1,1,routgrads(ng)%r2,nrec,cvar)
               !---------------------------------------------------------------------------!
            case (3)
               !---------------------------------------------------------------------------!
               !    Three-dimensional array.                                               !
               !---------------------------------------------------------------------------!


               !----- Alias for variable names. -------------------------------------------!
               cvar     = trim(vp(iv))
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Decide which vertical levels to use.                                  !
               !---------------------------------------------------------------------------!
               select case (ipresslev)
               case (0)
                  !------------------------------------------------------------------------!
                  !      Native coordinates.                                               !
                  !------------------------------------------------------------------------!


                  !----- Set the number of levels. ----------------------------------------!
                  nzvp(iv) = iep_nz(ng)
                  !------------------------------------------------------------------------!


                  !----- Adjust projection to GrADS. --------------------------------------!
                  call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng),nzvp(iv)   &
                                         ,nxgrads(ng),nygrads(ng),routgrads(ng)%rmi        &
                                         ,routgrads(ng)%iinf,routgrads(ng)%jinf            &
                                         ,rout(ng)%r3,routgrads(ng)%r3,proj)
                  !------------------------------------------------------------------------!

                  !----- Dump array to output file. ---------------------------------------!
                  call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)       &
                                ,nya(ng),nyb(ng),2,zlevmax(ng)+1,routgrads(ng)%r3,nrec,cvar)
                  !------------------------------------------------------------------------!
               case (1)
                  !------------------------------------------------------------------------!
                  !      Pressure levels.                                                  !
                  !------------------------------------------------------------------------!


                  !----- Set the number of levels. ----------------------------------------!
                  nzvp(iv) = inplevs
                  !------------------------------------------------------------------------!


                  !----- Interpolate to pressure levels. ----------------------------------!
                  write(unit=*,fmt='(a)') "        - PTRANSVAR"
                  call  ptransvar(rout(ng)%r3,iep_nx(ng),iep_ny(ng),iep_nz(ng),nzvp(iv)    &
                                 ,iplevs,rout(ng)%exner,dep_zlev(:,ng),rout(ng)%zplev      &
                                 ,rout(ng)%topo)
                  !------------------------------------------------------------------------!


                  !----- Adjust projection to GrADS. --------------------------------------!
                  write(unit=*,fmt='(a)') "        - PROJ_RAMS_TO_GRADS"
                  call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng),nzvp(iv)   &
                                         ,nxgrads(ng),nygrads(ng),routgrads(ng)%rmi        &
                                         ,routgrads(ng)%iinf,routgrads(ng)%jinf            &
                                         ,rout(ng)%r3,routgrads(ng)%r3,proj)
                  !------------------------------------------------------------------------!

                  !----- Dump array to output file. ---------------------------------------!
                  write(unit=*,fmt='(a)') "        - ED_PUTVAR"
                  call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)       &
                                ,nya(ng),nyb(ng),1,nzvp(iv),routgrads(ng)%r3,nrec,cvar)
                  !------------------------------------------------------------------------!
               case (2)
                  !------------------------------------------------------------------------!
                  !      Height levels.                                                    !
                  !------------------------------------------------------------------------!


                  !----- Set the number of levels. ----------------------------------------!
                  nzvp(iv) = inplevs
                  !------------------------------------------------------------------------!


                  !----- Interpolate to height levels. ------------------------------------!
                  call  ctransvar(iep_nx(ng),iep_ny(ng),iep_nz(ng),rout(ng)%r3             &
                                 ,rout(ng)%topo,nzvp(iv),iplevs,myztn(:,ng)                &
                                 ,myzmn(mynnzp(1)-1,1))
                  !------------------------------------------------------------------------!


                  !----- Adjust projection to GrADS. --------------------------------------!
                  call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng),nzvp(iv)   &
                                         ,nxgrads(ng),nygrads(ng),routgrads(ng)%rmi        &
                                         ,routgrads(ng)%iinf,routgrads(ng)%jinf            &
                                         ,rout(ng)%r3,routgrads(ng)%r3,proj)
                  !------------------------------------------------------------------------!


                  !----- Dump array to output file. ---------------------------------------!
                  call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)       &
                                ,nya(ng),nyb(ng),1,nzvp(iv),routgrads(ng)%r3,nrec,cvar)
                  !------------------------------------------------------------------------!
               case (3)
                  !------------------------------------------------------------------------!
                  !      Selected sigma-z levels.                                          !
                  !------------------------------------------------------------------------!


                  !----- Set the number of levels. ----------------------------------------!
                  nzvp(iv) = inplevs
                  !------------------------------------------------------------------------!


                  !----- Pick only the levels we are interested in. -----------------------!
                  call  select_sigmaz(iep_nx(ng),iep_ny(ng),iep_nz(ng),rout(ng)%r3         &
                                     ,inplevs,iplevs)
                  !------------------------------------------------------------------------!


                  !----- Adjust projection to GrADS. --------------------------------------!
                  call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng),nzvp(iv)   &
                                         ,nxgrads(ng),nygrads(ng),routgrads(ng)%rmi        &
                                         ,routgrads(ng)%iinf,routgrads(ng)%jinf            &
                                         ,rout(ng)%r3,routgrads(ng)%r3,proj)
                  !------------------------------------------------------------------------!


                  !----- Dump array to output file. ---------------------------------------!
                  call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)       &
                                ,nya(ng),nyb(ng),1,nzvp(iv),routgrads(ng)%r3,nrec,cvar)
                  !------------------------------------------------------------------------!

               end select

            case (6)
               !---------------------------------------------------------------------------!
               !    Four-dimensional array, fourth dimension is the cloud-level.  We will  !
               ! save these as independent variables.  Like in the 3-D case, we must also  !
               ! decide which vertical levels to plot.                                     !
               !---------------------------------------------------------------------------!

               !----- Update the number of variables. -------------------------------------!
               if (nfn == 1) nnvp = nnvp + iep_nc - 1
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Decide which vertical levels to use.                                  !
               !---------------------------------------------------------------------------!
               select case (ipresslev)
               case (0)
                  !------------------------------------------------------------------------!
                  !      Native coordinates.                                               !
                  !------------------------------------------------------------------------!


                  !----- Set the number of levels. ----------------------------------------!
                  nzvp(iv) = iep_nz(ng)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    Loop over clouds.                                                   !
                  !------------------------------------------------------------------------!
                  do ic = 1,iep_nc
                     write(cvar,fmt='(a,i2.2)') trim(vp(iv)),ic

                     !----- Convert the 4D array into a 3D. -------------------------------!
                     call s4d_to_3d(iep_nx(ng),iep_ny(ng),iep_nz(ng),iep_nc,ic             &
                                   ,rout(ng)%r6,rout(ng)%r3)
                     !---------------------------------------------------------------------!


                     !----- Adjust projection to GrADS. -----------------------------------!
                     call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng)         &
                                            ,nzvp(iv),nxgrads(ng),nygrads(ng)              &
                                            ,routgrads(ng)%rmi,routgrads(ng)%iinf          &
                                            ,routgrads(ng)%jinf,rout(ng)%r3                &
                                            ,routgrads(ng)%r3,proj)
                     !---------------------------------------------------------------------!

                     !----- Dump array to output file. ------------------------------------!
                     call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)    &
                                   ,nya(ng),nyb(ng),2,zlevmax(ng)+1,routgrads(ng)%r3,nrec  &
                                   ,cvar)
                     !---------------------------------------------------------------------!
                  end do
                  !------------------------------------------------------------------------!

               case (1)
                  !------------------------------------------------------------------------!
                  !      Pressure levels.                                                  !
                  !------------------------------------------------------------------------!


                  !----- Set the number of levels. ----------------------------------------!
                  nzvp(iv) = inplevs
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    Loop over clouds.                                                   !
                  !------------------------------------------------------------------------!
                  do ic = 1,iep_nc
                     write(cvar,fmt='(a,i2.2)') trim(vp(iv)),ic

                     !----- Convert the 4D array into a 3D. -------------------------------!
                     call s4d_to_3d(iep_nx(ng),iep_ny(ng),iep_nz(ng),iep_nc,ic             &
                                   ,rout(ng)%r6,rout(ng)%r3)
                     !---------------------------------------------------------------------!

                     !----- Interpolate to pressure levels. -------------------------------!
                     call  ptransvar(rout(ng)%r3,iep_nx(ng),iep_ny(ng),iep_nz(ng),inplevs  &
                                    ,iplevs,rout(ng)%exner,dep_zlev(:,ng),rout(ng)%zplev   &
                                    ,rout(ng)%topo)
                     !---------------------------------------------------------------------!


                     !----- Adjust projection to GrADS. -----------------------------------!
                     call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng)         &
                                            ,nzvp(iv),nxgrads(ng),nygrads(ng)              &
                                            ,routgrads(ng)%rmi,routgrads(ng)%iinf          &
                                            ,routgrads(ng)%jinf,rout(ng)%r3                &
                                            ,routgrads(ng)%r3,proj)
                     !---------------------------------------------------------------------!

                     !----- Dump array to output file. ------------------------------------!
                     call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)    &
                                   ,nya(ng),nyb(ng),1,nzvp(iv),routgrads(ng)%r3,nrec,cvar)
                     !---------------------------------------------------------------------!
                  end do
                  !------------------------------------------------------------------------!

               case (2)
                  !------------------------------------------------------------------------!
                  !      Height levels.                                                    !
                  !------------------------------------------------------------------------!


                  !----- Set the number of levels. ----------------------------------------!
                  nzvp(iv) = inplevs
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    Loop over clouds.                                                   !
                  !------------------------------------------------------------------------!
                  do ic = 1,iep_nc
                     write(cvar,fmt='(a,i2.2)') trim(vp(iv)),ic

                     !----- Convert the 4D array into a 3D. -------------------------------!
                     call s4d_to_3d(iep_nx(ng),iep_ny(ng),iep_nz(ng),iep_nc,ic             &
                                   ,rout(ng)%r6,rout(ng)%r3)
                     !---------------------------------------------------------------------!

                     !----- Interpolate to height levels. ---------------------------------!
                     call  ctransvar(iep_nx(ng),iep_ny(ng),iep_nz(ng),rout(ng)%r3          &
                                    ,rout(ng)%topo,inplevs,iplevs,myztn(:,ng)              &
                                    ,myzmn(mynnzp(1)-1,1))
                     !---------------------------------------------------------------------!


                     !----- Adjust projection to GrADS. -----------------------------------!
                     call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng)         &
                                            ,nzvp(iv),nxgrads(ng),nygrads(ng)              &
                                            ,routgrads(ng)%rmi,routgrads(ng)%iinf          &
                                            ,routgrads(ng)%jinf,rout(ng)%r3                &
                                            ,routgrads(ng)%r3,proj)
                     !---------------------------------------------------------------------!


                     !----- Dump array to output file. ------------------------------------!
                     call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)    &
                                   ,nya(ng),nyb(ng),1,nzvp(iv),routgrads(ng)%r3,nrec,cvar)
                     !---------------------------------------------------------------------!
                  end do
                  !------------------------------------------------------------------------!
               case (3)
                  !------------------------------------------------------------------------!
                  !      Selected sigma-z levels.                                          !
                  !------------------------------------------------------------------------!


                  !----- Set the number of levels. ----------------------------------------!
                  nzvp(iv) = inplevs
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    Loop over clouds.                                                   !
                  !------------------------------------------------------------------------!
                  do ic = 1,iep_nc
                     write(cvar,fmt='(a,i2.2)') trim(vp(iv)),ic

                     !----- Convert the 4D array into a 3D. -------------------------------!
                     call s4d_to_3d(iep_nx(ng),iep_ny(ng),iep_nz(ng),iep_nc,ic             &
                                   ,rout(ng)%r6,rout(ng)%r3)
                     !---------------------------------------------------------------------!

                     !----- Pick only the levels we are interested in. --------------------!
                     call  select_sigmaz(iep_nx(ng),iep_ny(ng),iep_nz(ng),rout(ng)%r3      &
                                        ,inplevs,iplevs)
                     !---------------------------------------------------------------------!


                     !----- Adjust projection to GrADS. -----------------------------------!
                     call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng)         &
                                            ,nzvp(iv),nxgrads(ng),nygrads(ng)              &
                                            ,routgrads(ng)%rmi,routgrads(ng)%iinf          &
                                            ,routgrads(ng)%jinf,rout(ng)%r3                &
                                            ,routgrads(ng)%r3,proj)
                     !---------------------------------------------------------------------!


                     !----- Dump array to output file. ------------------------------------!
                     call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)    &
                                   ,nya(ng),nyb(ng),1,nzvp(iv),routgrads(ng)%r3,nrec,cvar)
                     !---------------------------------------------------------------------!
                  end do
                  !------------------------------------------------------------------------!

               end select
            case (7)
               !---------------------------------------------------------------------------!
               !    Three-dimensional array, third dimension is the patch-level.  We will  !
               ! save these as independent variables.                                      !
               !---------------------------------------------------------------------------!

               !----- Update the number of variables. -------------------------------------!
               if (nfn == 1) nnvp = nnvp + iep_np - 1
               !---------------------------------------------------------------------------!

               !----- Set the number of levels. -------------------------------------------!
               nzvp(iv) = iep_np
               !---------------------------------------------------------------------------!


               !----- Adjust projection to GrADS. -----------------------------------------!
               call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng),nzvp(iv)      &
                                      ,nxgrads(ng),nygrads(ng),routgrads(ng)%rmi           &
                                      ,routgrads(ng)%iinf,routgrads(ng)%jinf,rout(ng)%r7   &
                                      ,routgrads(ng)%r7,proj)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Loop over all patches, and write the 2-D arrays.                       !
               !---------------------------------------------------------------------------!
               do ip=1,iep_np
                  write(cvar,fmt='(a,i2.2)') trim(vp(iv)),ip

                  call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)       &
                                ,nya(ng),nyb(ng),ip,ip,routgrads(ng)%r7,nrec,cvar)
               end do
               !---------------------------------------------------------------------------!

            case (8)
               !---------------------------------------------------------------------------!
               !    Soil variable that has layers and patches.                             !
               !---------------------------------------------------------------------------!


               !----- Update the number of variables. -------------------------------------!
               if (nfn == 1) nnvp = nnvp + iep_np - 1
               !---------------------------------------------------------------------------!


               !----- Set the number of levels. -------------------------------------------!
               nzvp(iv) = iep_ng
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Loop over patches.                                                     !
               !---------------------------------------------------------------------------!
               do ip = 1,iep_np
                  write(cvar,fmt='(a,i2.2)') trim(vp(iv)),ip

                  !----- Convert the 4D array into a 3D. ----------------------------------!
                  call s4d_to_3d(iep_nx(ng),iep_ny(ng),nzvp(iv),iep_np,ip                  &
                                ,rout(ng)%r8,rout(ng)%r10)
                  !------------------------------------------------------------------------!


                  !----- Adjust projection to GrADS. --------------------------------------!
                  call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng),nzvp(iv)   &
                                         ,nxgrads(ng),nygrads(ng),routgrads(ng)%rmi        &
                                         ,routgrads(ng)%iinf,routgrads(ng)%jinf            &
                                         ,rout(ng)%r10,routgrads(ng)%r10,proj)
                  !------------------------------------------------------------------------!


                  !----- Dump array to output file. ---------------------------------------!
                  call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)       &
                                ,nya(ng),nyb(ng),1,iep_ng,routgrads(ng)%r10,nrec,cvar)
                  !------------------------------------------------------------------------!
               end do
               !---------------------------------------------------------------------------!

            case (9)
               !---------------------------------------------------------------------------!
               !    Three-dimensional array, third dimension is the cloud-level.  We will  !
               ! save these as independent variables.                                      !
               !---------------------------------------------------------------------------!


               !----- Update the number of variables. -------------------------------------!
               if (nfn == 1) nnvp = nnvp + iep_nc - 1
               !---------------------------------------------------------------------------!


               !----- Set the number of levels. -------------------------------------------!
               nzvp(iv) = iep_nc
               !---------------------------------------------------------------------------!


               !----- Adjust projection to GrADS. -----------------------------------------!
               call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng),nzvp(iv)      &
                                      ,nxgrads(ng),nygrads(ng),routgrads(ng)%rmi           &
                                      ,routgrads(ng)%iinf,routgrads(ng)%jinf,rout(ng)%r9   &
                                      ,routgrads(ng)%r9,proj)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Loop over all clouds, and write the 2-D arrays.                        !
               !---------------------------------------------------------------------------!
               do ic=1,iep_nc
                  write(cvar,fmt='(a,i2.2)') trim(vp(iv)),ic

                  call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)       &
                                ,nya(ng),nyb(ng),ic,ic,routgrads(ng)%r9,nrec,cvar)
               end do
               !---------------------------------------------------------------------------!

            case (10)
               !---------------------------------------------------------------------------!
               !    Soil variable that has layers but no patches.                          !
               !---------------------------------------------------------------------------!
               cvar=trim(vp(iv))


               !----- Set the number of levels. -------------------------------------------!
               nzvp(iv) = iep_ng
               !---------------------------------------------------------------------------!


               !----- Adjust projection to GrADS. -----------------------------------------!
               call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng),nzvp(iv)      &
                                      ,nxgrads(ng),nygrads(ng),routgrads(ng)%rmi           &
                                      ,routgrads(ng)%iinf,routgrads(ng)%jinf,rout(ng)%r10  &
                                      ,routgrads(ng)%r10,proj)
               !---------------------------------------------------------------------------!


               !----- Dump array to output file. ------------------------------------------!
               call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)          &
                             ,nya(ng),nyb(ng),1,nzvp(iv),routgrads(ng)%r10,nrec,cvar)
               !---------------------------------------------------------------------------!

            case default
               !------ Invalid variable, remove one from the total count. -----------------!
               if (nfn == 1) nnvp=nnvp-1
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!
         end do varloop
         !---------------------------------------------------------------------------------!
      end do fileloop
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Close the binary file now.                                                     !
      !------------------------------------------------------------------------------------!
      close (unit=19,status='keep')
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Define the number of output grid points in X and Y.                            !
      !------------------------------------------------------------------------------------!
      nxpg = nxb(ng) - nxa(ng) + 1
      nypg= nyb(ng) - nya(ng) + 1
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Open the CTL file and start writing the header.                                !
      !------------------------------------------------------------------------------------!
      open (unit=20,file=trim(gprefix)//'_g'//cgrid//'.ctl',status='replace'               &
           ,action='write')
      write(unit=20,fmt='(a)') 'dset ^'//trim(gprefix)//'_g'//cgrid//'.gra'
      write(unit=20,fmt='(a,1x,es10.3)') 'undef',undefflg
      write(unit=20,fmt='(a)') 'title BRAMS-4.0.6 output'
      write(unit=20,fmt='(a,1x,i5,1x,a,2(1x,f14.5))')  'xdef',nxpg,'linear'                &
                                                      ,dep_glon(1,ng),dep_glon(2,ng)
      write(unit=20,fmt='(a,1x,i5,1x,a,2(1x,f14.5))')  'ydef',nypg,'linear'                &
                                                      ,dep_glat(1,ng),dep_glat(2,ng)

      !------------------------------------------------------------------------------------!
      !     Write the vertical coordinates.                                                !
      !------------------------------------------------------------------------------------!
      select case (ipresslev)
      case (0)
         !----- Native coordinates.  Break it in case there are more than 15 lines. -------!
         if (zlevmax(ng) > 15) then
            write (unit=20,fmt='(a,1x,i5,1x,a,15(1x,f14.5))') 'zdef',zlevmax(ng),'levels'  &
                                                             ,(dep_zlev(n,ng),n=2,15)
            write (unit=20,fmt='(200(1x,f14.5))') (dep_zlev(n,ng),n=16,zlevmax(ng)+1)
         else
            write (unit=20,fmt='(a,1x,i5,1x,a,200(1x,f14.5))') 'zdef',zlevmax(ng),'levels' &
                                                        ,(dep_zlev(n,ng),n=2,zlevmax(ng)+1)
         end if
         !---------------------------------------------------------------------------------!
      case (1,2)
         !----- Pressure or height coordinates.  ------------------------------------------!
         write (unit=20,fmt='(a,1x,i5,1x,a,200(1x,f14.5))') 'zdef',inplevs,'levels'        &
                                                           ,(iplevs(n)*1.0,n=1,inplevs)
         !---------------------------------------------------------------------------------!
      case (3)
         !----- Selected sigma-z coordinates.  --------------------------------------------!
         write (unit=20,fmt='(a,1x,i5,1x,a,200(1x,f14.5))') 'zdef',inplevs,'levels'        &
                                                      ,(dep_zlev(iplevs(n),ng),n=1,inplevs)
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Write time.                                                                    !
      !------------------------------------------------------------------------------------!
      write(unit=20,fmt='(a,1x,i5,3(1x,a))') 'tdef',nfiles,'linear',chdate,chstep
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Loop over variables.                                                           !
      !------------------------------------------------------------------------------------!
      write(unit=20,fmt='(a,1x,i5)') 'vars',nnvp
      varoutloop: do iv=1,nvp
         select case (ndim(iv))
         case (2)
            !------------------------------------------------------------------------------!
            !      2 D variable, no height reference.                                      !
            !------------------------------------------------------------------------------!
            write(unit=20,fmt='(a,2(1x,i5),4(1x,a))')  vp(iv),0,99,vpln(iv)                &
                                                      ,'[',vpun(iv),']'
         case (3)
            !------------------------------------------------------------------------------!
            !      3 D variable.  Check which vertical coordinate to use.                  !
            !------------------------------------------------------------------------------!
            select case (ipresslev)
            case (0)
               !----- Native coordinates. -------------------------------------------------!
               write(unit=20,fmt='(a,2(1x,i5),4(1x,a))') vp(iv),zlevmax(ng),99,vpln(iv)    &
                                                        ,'[',vpun(iv),']'
            case default
               !----- Other coordinates. --------------------------------------------------!
               write(unit=20,fmt='(a,2(1x,i5),4(1x,a))') vp(iv),inplevs,99,vpln(iv)        &
                                                        ,'[',vpun(iv),']'
            end select
            !------------------------------------------------------------------------------!
         case (6)
            !------------------------------------------------------------------------------!
            !      4 D variable.  Make one entry per cloud, and check which vertical       !
            ! coordinate to use.                                                           !
            !------------------------------------------------------------------------------!
            do ic = 1, iep_nc
               write(cldnumber,fmt='(i2.2)') ic
               tmpvar  = trim(vp(iv))//cldnumber
               tmpdesc = trim(vpln(iv))//': Cloud # '//cldnumber
               !---------------------------------------------------------------------------!
               !      Check which vertical coordinate to use.                              !
               !---------------------------------------------------------------------------!
               select case (ipresslev)
               case (0)
                  !----- Native coordinates. ----------------------------------------------!
                  write(unit=20,fmt='(a,2(1x,i5),4(1x,a))') tmpvar,zlevmax(ng),99,tmpdesc  &
                                                           ,'[',vpun(iv),']'
               case default
                  !----- Other coordinates. -----------------------------------------------!
                  write(unit=20,fmt='(a,2(1x,i5),4(1x,a))') tmpvar,inplevs,99,tmpdesc      &
                                                           ,'[',vpun(iv),']'
               end select
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!

         case (7)
            !------------------------------------------------------------------------------!
            !      3-D variable.  Make one entry per patch.                                !
            !------------------------------------------------------------------------------!
            do ip = 1, iep_np
               write(patchnumber,fmt='(i2.2)') ip
               tmpvar  = trim(vp(iv))//patchnumber
               tmpdesc = trim(vpln(iv))//': Patch # '//patchnumber
               write(unit=20,fmt='(a,2(1x,i5),4(1x,a))') tmpvar,0,99,tmpdesc               &
                                                        ,'[',vpun(iv),']'
            end do
            !------------------------------------------------------------------------------!

         case (8)
            !------------------------------------------------------------------------------!
            !      4-D variable.  Make one entry per patch, with vertical being number of  !
            ! soil layers.                                                                 !
            !------------------------------------------------------------------------------!
            do ip = 1, iep_np
               write(patchnumber,fmt='(i2.2)') ip
               tmpvar  = trim(vp(iv))//patchnumber
               tmpdesc = trim(vpln(iv))//': Patch # '//patchnumber
               write(unit=20,fmt='(a,2(1x,i5),4(1x,a))') tmpvar,nzvp(iv),99,tmpdesc        &
                                                        ,'[',vpun(iv),']'
            end do
            !------------------------------------------------------------------------------!

         case (9)
            !------------------------------------------------------------------------------!
            !      3-D variable.  Make one entry per cloud.                                !
            !------------------------------------------------------------------------------!
            do ic = 1, iep_nc
               write(cldnumber,fmt='(i2.2)') ic
               tmpvar  = trim(vp(iv))//cldnumber
               tmpdesc = trim(vpln(iv))//': Cloud # '//cldnumber
               write(unit=20,fmt='(a,2(1x,i5),4(1x,a))') tmpvar,0,99,tmpdesc               &
                                                        ,'[',vpun(iv),']'
            end do
            !------------------------------------------------------------------------------!

         case (10)
            !------------------------------------------------------------------------------!
            !      3-D variable, soil layers with no patch.                                !
            !------------------------------------------------------------------------------!
            write(unit=20,fmt='(a,2(1x,i5),4(1x,a))')  vp(iv),nzvp(iv),99,vpln(iv)         &
                                                      ,'[',vpun(iv),']'
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      end do varoutloop
      !------------------------------------------------------------------------------------!
      write(unit=20,fmt='(a)') 'endvars'
      close(unit=20,status='keep')
      write (unit=*,fmt='(92a)'    )    ('=',n=1,92)
      write (unit=*,fmt='(a)'      )    ' '
      !------------------------------------------------------------------------------------!
   end do gridloop
   !---------------------------------------------------------------------------------------!

   write(*,'(a)') ' ------ Ramspost execution ends ------'
   stop
end program ramspost
!==========================================================================================!
!==========================================================================================!
