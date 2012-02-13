!==========================================================================================!
!==========================================================================================!
!                                                                                          !
!                         RAMSPOST - RAMS Post Processor for GrADS.                        !
!                                                                                          !
!------------------------------------------------------------------------------------------!
program ramspost
   use rpost_coms
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
   integer                                                           :: nvp
   integer                                                           :: nfiles
   integer               , dimension(maxvars)                        :: nzvp
   integer                                                           :: nrec
   integer                                                           :: ipresslev
   integer               , dimension(nplmax)                         :: iplevs
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
   integer                                                           :: inplevsef
   integer                                                           :: nxpg
   integer                                                           :: nypg
   integer                                                           :: nzpg
   integer               , dimension(maxgrds)                        :: nxgrads
   integer               , dimension(maxgrds)                        :: nygrads
   integer               , dimension(maxgx,maxgy)                    :: iinf
   integer               , dimension(maxgx,maxgy)                    :: jinf
   integer               , dimension(maxgrds)                        :: nxa
   integer               , dimension(maxgrds)                        :: nxb
   integer               , dimension(maxgrds)                        :: nya
   integer               , dimension(maxgrds)                        :: nyb
   integer               , dimension(maxgrds)                        :: iep_nx
   integer               , dimension(maxgrds)                        :: iep_ny
   integer               , dimension(maxgrds)                        :: iep_nz
   real                  , dimension(maxgx,maxgy,4)                  :: rmi
   real                  , dimension(maxgx,maxgy,nzpmax)             :: routgrads
   real                  , dimension(nxpmax,nypmax,nzepmax)          :: a
   real                  , dimension(nxpmax,nypmax,nzepmax)          :: b
   real                  , dimension(nxpmax,nypmax,nzepmax)          :: rout
   real                  , dimension(nxpmax,nypmax,nplmax)           :: zplev
   real                  , dimension(nxpmax,nypmax)                  :: mytopo
   real                  , dimension(nxpmax,nypmax,nzpmax)           :: mypi
   real                  , dimension(nxpmax,nypmax)                  :: rlat
   real                  , dimension(nxpmax,nypmax)                  :: rlon
   real                  , dimension(maxgrds)                        :: lati
   real                  , dimension(maxgrds)                        :: latf
   real                  , dimension(maxgrds)                        :: loni
   real                  , dimension(maxgrds)                        :: lonf
   real                  , dimension(nxpmax,nypmax,nzgmax,maxpatch)  :: a2
   real                  , dimension(nxpmax,nypmax,nzgmax,maxpatch)  :: rout2
   real                  , dimension(nxpmax,nypmax,nzpmax,maxclouds) :: a6
   real                  , dimension(nxpmax,nypmax,nzpmax,maxclouds) :: rout6
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
   !      Grid loop.                                                                       !
   !---------------------------------------------------------------------------------------!
   gridloop: do ng=1,iep_ngrids
      write (unit=*,fmt='(92a)'    )    ('=',n=1,92)
      write (unit=*,fmt='(a)'      )    ' '
      write (unit=*,fmt='(a,1x,i5)')    ' + Writing Grid ',ng


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

      call ep_getvar('lat',rlat,a,b,iep_nx(ng),iep_ny(ng),1,ng,cfln,vpln(iv)               &
                    ,vpun(iv),n,iep_np,iep_nc,iep_ng,a2,rout2,a6,rout6)
      write(unit=*,fmt='(a,1x,i5)') '     # Output variable type:  ',n

      write(unit=*,fmt='(a)') ' '
      write(unit=*,fmt='(2(a,1x))') '   - Variable:  ','lon'

      call ep_getvar('lon',rlon,a,b,iep_nx(ng),iep_ny(ng),1,ng,cfln ,vpln(iv)              &
                    ,vpun(iv),n,iep_np,iep_nc,iep_ng,a2,rout2,a6,rout6)
      write(unit=*,fmt='(a,1x,i5)') '     # Output variable type:  ',n
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Find the dimensions of x and y domains depending on the sought projection, and !
      ! find the mapping to interpolate values to the regular lon-lat grid if needed.      !
      !------------------------------------------------------------------------------------!
      call geo_grid(iep_nx(ng),iep_ny(ng),rlat,rlon,dep_glon(1,ng),dep_glon(2,ng)          &
                   ,dep_glat(1,ng),dep_glat(2,ng),rlatmin,rlatmax,rlonmin,rlonmax          &
                   ,nxgrads(ng),nygrads(ng),proj)
      call array_interpol(ng,nxgrads(ng),nygrads(ng),iep_nx(ng),iep_ny(ng),dep_glat(1,ng)  &
                         ,dep_glat(2,ng),dep_glon(1,ng),dep_glon(2,ng),iinf,jinf,rmi,proj)

      call define_lim(ng,nxgrads(ng),nygrads(ng),dep_glat(1,ng),dep_glat(2,ng)             &
                     ,dep_glon(1,ng),dep_glon(2,ng),lati(ng),latf(ng),loni(ng),lonf(ng)    &
                     ,nxa(ng),nxb(ng),nya(ng),nyb(ng),proj,iep_nx(ng),iep_ny(ng),rlat,rlon)
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
            call ep_getvar('topo',mytopo,a,b,iep_nx(ng),iep_ny(ng),1,ng,cfln               &
                          ,vpln(iv),vpun(iv),n,iep_np,iep_nc,iep_ng,a2,rout2,a6,rout6)
            write(unit=*,fmt='(a,1x,i5)') '       # Output variable type:  ',n
            !------------------------------------------------------------------------------!

            !----- Load Exner function. ---------------------------------------------------!
            write(unit=*,fmt='(a)') ' '
            write(unit=*,fmt='(2(a,1x))') '     * Variable:  ','pi'
            call ep_getvar('pi',mypi,a,b,iep_nx(ng),iep_ny(ng),iep_nz(ng),ng,cfln,vpln(iv) &
                          ,vpun(iv),n,iep_np,iep_nc,iep_ng,a2,rout2,a6,rout6)
            write(unit=*,fmt='(a,1x,i5)') '       # Output variable type:  ',n
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Loop over all other variables.                                              !
         !---------------------------------------------------------------------------------!
         varloop: do iv=1,nvp
            write(unit=*,fmt='(a)') ' '
            write(unit=*,fmt='(2(a,1x))') '     * Variable:  ',trim(vp(iv))
            call ep_getvar(trim(vp(iv)),rout,a,b,iep_nx(ng),iep_ny(ng),iep_nz(ng),ng,cfln  &
                          ,vpln(iv),vpun(iv),ndim(iv),iep_np,iep_nc,iep_ng,a2,rout2        &
                          ,a6,rout6)
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
               !---------------------------------------------------------------------------!


               !----- Adjust projection to GrADS. -----------------------------------------!
               call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng),nzvp(iv)      &
                                      ,nxgrads(ng),nygrads(ng),rmi,iinf,jinf,rout          &
                                      ,routgrads,rlat,rlon,proj)
               !---------------------------------------------------------------------------!


               !----- Dump array to output file. ------------------------------------------!
               call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng),nya(ng)  &
                             ,nyb(ng),1,1,routgrads,nrec)
               !---------------------------------------------------------------------------!
            case (3)
               !---------------------------------------------------------------------------!
               !    Three-dimensional array.                                               !
               !---------------------------------------------------------------------------!


               !----- Set the number of levels. -------------------------------------------!
               nzvp(iv) = iep_nz(ng)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Decide which vertical levels to use.                                  !
               !---------------------------------------------------------------------------!
               select case (ipresslev)
               case (0)
                  !------------------------------------------------------------------------!
                  !      Native coordinates.                                               !
                  !------------------------------------------------------------------------!

                  !----- Adjust projection to GrADS. --------------------------------------!
                  call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng),nzvp(iv)   &
                                         ,nxgrads(ng),nygrads(ng),rmi,iinf,jinf,rout       &
                                         ,routgrads,rlat,rlon,proj)
                  !------------------------------------------------------------------------!

                  !----- Dump array to output file. ---------------------------------------!
                  call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)       &
                                ,nya(ng),nyb(ng),2,zlevmax(ng)+1,routgrads,nrec)
                  !------------------------------------------------------------------------!
               case (1)
                  !------------------------------------------------------------------------!
                  !      Pressure levels.                                                  !
                  !------------------------------------------------------------------------!

                  !----- Set the number of levels. ----------------------------------------!
                  inplevsef = inplevs
                  !------------------------------------------------------------------------!

                  !----- Interpolate to pressure levels. ----------------------------------!
                  call  ptransvar(rout,iep_nx(ng),iep_ny(ng),nzvp(iv),inplevs,iplevs,mypi  &
                                 ,dep_zlev(1,ng),zplev,mytopo)
                  !------------------------------------------------------------------------!


                  !----- Adjust projection to GrADS. --------------------------------------!
                  call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng),nzvp(iv)   &
                                         ,nxgrads(ng),nygrads(ng),rmi,iinf,jinf,rout       &
                                         ,routgrads,rlat,rlon,proj)
                  !------------------------------------------------------------------------!

                  !----- Dump array to output file. ---------------------------------------!
                  call ep_putvar(19,nxgrads(ng),nygrads(ng),inplevsef,nxa(ng),nxb(ng)      &
                                ,nya(ng),nyb(ng),1,inplevsef,routgrads,nrec)
                  !------------------------------------------------------------------------!
               case (2)
                  !------------------------------------------------------------------------!
                  !      Height levels.                                                    !
                  !------------------------------------------------------------------------!

                  !----- Set the number of levels. ----------------------------------------!
                  inplevsef = inplevs
                  !------------------------------------------------------------------------!

                  !----- Interpolate to height levels. ------------------------------------!
                  call  ctransvar(iep_nx(ng),iep_ny(ng),iep_nz(ng),rout,mytopo,inplevs     &
                                 ,iplevs,myztn(1,ng),myzmn(mynnzp(1)-1,1))
                  !------------------------------------------------------------------------!


                  !----- Adjust projection to GrADS. --------------------------------------!
                  call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng),nzvp(iv)   &
                                         ,nxgrads(ng),nygrads(ng),rmi,iinf,jinf,rout       &
                                         ,routgrads,rlat,rlon,proj)
                  !------------------------------------------------------------------------!


                  !----- Dump array to output file. ---------------------------------------!
                  call ep_putvar(19,nxgrads(ng),nygrads(ng),inplevsef,nxa(ng),nxb(ng)      &
                                ,nya(ng),nyb(ng),1,inplevsef,routgrads,nrec)
                  !------------------------------------------------------------------------!
               case (3)
                  !------------------------------------------------------------------------!
                  !      Selected sigma-z levels.                                          !
                  !------------------------------------------------------------------------!

                  !----- Set the number of levels. ----------------------------------------!
                  inplevsef = inplevs
                  !------------------------------------------------------------------------!

                  !----- Pick only the levels we are interested in. -----------------------!
                  call  select_sigmaz(iep_nx(ng),iep_ny(ng),iep_nz(ng),rout,inplevs,iplevs)
                  !------------------------------------------------------------------------!


                  !----- Adjust projection to GrADS. --------------------------------------!
                  call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng),nzvp(iv)   &
                                         ,nxgrads(ng),nygrads(ng),rmi,iinf,jinf,rout       &
                                         ,routgrads,rlat,rlon,proj)
                  !------------------------------------------------------------------------!


                  !----- Dump array to output file. ---------------------------------------!
                  call ep_putvar(19,nxgrads(ng),nygrads(ng),inplevsef,nxa(ng),nxb(ng)      &
                                ,nya(ng),nyb(ng),1,inplevsef,routgrads,nrec)
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


               !----- Set the number of levels. -------------------------------------------!
               nzvp(iv) = iep_nz(ng)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Decide which vertical levels to use.                                  !
               !---------------------------------------------------------------------------!
               select case (ipresslev)
               case (0)
                  !------------------------------------------------------------------------!
                  !      Native coordinates.                                               !
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    Loop over clouds.                                                   !
                  !------------------------------------------------------------------------!
                  do ic = 1,iep_nc
                     !----- Convert the 4D array into a 3D. -------------------------------!
                     call S4d_to_3d(iep_nx(ng),iep_ny(ng),nzvp(iv),iep_nc,ic,rout6,rout)
                     !---------------------------------------------------------------------!


                     !----- Adjust projection to GrADS. -----------------------------------!
                     call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng)         &
                                            ,nzvp(iv),nxgrads(ng),nygrads(ng),rmi,iinf     &
                                            ,jinf,rout,routgrads,rlat,rlon,proj)
                     !---------------------------------------------------------------------!

                     !----- Dump array to output file. ------------------------------------!
                     call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)    &
                                ,nya(ng),nyb(ng),2,zlevmax(ng)+1,routgrads,nrec)
                     !---------------------------------------------------------------------!
                  end do
                  !------------------------------------------------------------------------!

               case (1)
                  !------------------------------------------------------------------------!
                  !      Pressure levels.                                                  !
                  !------------------------------------------------------------------------!

                  !----- Set the number of levels. ----------------------------------------!
                  inplevsef = inplevs
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    Loop over clouds.                                                   !
                  !------------------------------------------------------------------------!
                  do ic = 1,iep_nc
                     !----- Convert the 4D array into a 3D. -------------------------------!
                     call S4d_to_3d(iep_nx(ng),iep_ny(ng),nzvp(iv),iep_nc,ic,rout6,rout)
                     !---------------------------------------------------------------------!

                     !----- Interpolate to pressure levels. -------------------------------!
                     call  ptransvar(rout,iep_nx(ng),iep_ny(ng),nzvp(iv),inplevs,iplevs    &
                                    ,mypi,dep_zlev(1,ng),zplev,mytopo)
                     !---------------------------------------------------------------------!


                     !----- Adjust projection to GrADS. -----------------------------------!
                     call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng)         &
                                            ,nzvp(iv),nxgrads(ng),nygrads(ng),rmi,iinf     &
                                            ,jinf,rout,routgrads,rlat,rlon,proj)
                     !---------------------------------------------------------------------!

                     !----- Dump array to output file. ------------------------------------!
                     call ep_putvar(19,nxgrads(ng),nygrads(ng),inplevsef,nxa(ng),nxb(ng)   &
                                   ,nya(ng),nyb(ng),1,inplevsef,routgrads,nrec)
                     !---------------------------------------------------------------------!
                  end do
                  !------------------------------------------------------------------------!

               case (2)
                  !------------------------------------------------------------------------!
                  !      Height levels.                                                    !
                  !------------------------------------------------------------------------!

                  !----- Set the number of levels. ----------------------------------------!
                  inplevsef = inplevs
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    Loop over clouds.                                                   !
                  !------------------------------------------------------------------------!
                  do ic = 1,iep_nc
                     !----- Convert the 4D array into a 3D. -------------------------------!
                     call S4d_to_3d(iep_nx(ng),iep_ny(ng),nzvp(iv),iep_nc,ic,rout6,rout)
                     !---------------------------------------------------------------------!

                     !----- Interpolate to height levels. ---------------------------------!
                     call  ctransvar(iep_nx(ng),iep_ny(ng),iep_nz(ng),rout,mytopo,inplevs  &
                                    ,iplevs,myztn(1,ng),myzmn(mynnzp(1)-1,1))
                     !---------------------------------------------------------------------!


                     !----- Adjust projection to GrADS. -----------------------------------!
                     call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng)         &
                                            ,nzvp(iv),nxgrads(ng),nygrads(ng),rmi,iinf     &
                                            ,jinf,rout,routgrads,rlat,rlon,proj)
                     !---------------------------------------------------------------------!


                     !----- Dump array to output file. ------------------------------------!
                     call ep_putvar(19,nxgrads(ng),nygrads(ng),inplevsef,nxa(ng),nxb(ng)   &
                                   ,nya(ng),nyb(ng),1,inplevsef,routgrads,nrec)
                     !---------------------------------------------------------------------!
                  end do
                  !------------------------------------------------------------------------!
               case (3)
                  !------------------------------------------------------------------------!
                  !      Selected sigma-z levels.                                          !
                  !------------------------------------------------------------------------!

                  !----- Set the number of levels. ----------------------------------------!
                  inplevsef = inplevs
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    Loop over clouds.                                                   !
                  !------------------------------------------------------------------------!
                  do ic = 1,iep_nc
                     !----- Convert the 4D array into a 3D. -------------------------------!
                     call S4d_to_3d(iep_nx(ng),iep_ny(ng),nzvp(iv),iep_nc,ic,rout6,rout)
                     !---------------------------------------------------------------------!

                     !----- Pick only the levels we are interested in. --------------------!
                     call  select_sigmaz(iep_nx(ng),iep_ny(ng),iep_nz(ng),rout,inplevs     &
                                        ,iplevs)
                     !---------------------------------------------------------------------!


                     !----- Adjust projection to GrADS. -----------------------------------!
                     call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng)         &
                                            ,nzvp(iv),nxgrads(ng),nygrads(ng),rmi,iinf     &
                                            ,jinf,rout,routgrads,rlat,rlon,proj)
                     !---------------------------------------------------------------------!


                     !----- Dump array to output file. ------------------------------------!
                     call ep_putvar(19,nxgrads(ng),nygrads(ng),inplevsef,nxa(ng),nxb(ng)   &
                                   ,nya(ng),nyb(ng),1,inplevsef,routgrads,nrec)
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
                                      ,nxgrads(ng),nygrads(ng),rmi,iinf,jinf,rout          &
                                      ,routgrads,rlat,rlon,proj)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Loop over all patches, and write the 2-D arrays.                       !
               !---------------------------------------------------------------------------!
               do ip=1,iep_np
                  call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)       &
                                   ,nya(ng),nyb(ng),ip,ip,routgrads,nrec)
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
                  !----- Convert the 4D array into a 3D. ----------------------------------!
                  call S4d_to_3d(iep_nx(ng),iep_ny(ng),nzvp(iv),iep_np,ip,rout2,rout)
                  !------------------------------------------------------------------------!


                  !----- Adjust projection to GrADS. --------------------------------------!
                  call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng),nzvp(iv)   &
                                         ,nxgrads(ng),nygrads(ng),rmi,iinf,jinf,rout       &
                                         ,routgrads,rlat,rlon,proj)
                  !------------------------------------------------------------------------!


                  !----- Dump array to output file. ---------------------------------------!
                  call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)       &
                                   ,nya(ng),nyb(ng),1,iep_ng,routgrads,nrec)
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
                                      ,nxgrads(ng),nygrads(ng),rmi,iinf,jinf,rout          &
                                      ,routgrads,rlat,rlon,proj)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Loop over all clouds, and write the 2-D arrays.                        !
               !---------------------------------------------------------------------------!
               do ic=1,iep_nc
                  call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)       &
                                   ,nya(ng),nyb(ng),ic,ic,routgrads,nrec)
               end do
               !---------------------------------------------------------------------------!

            case (10)
               !---------------------------------------------------------------------------!
               !    Soil variable that has layers but no patches (obsolete).               !
               !---------------------------------------------------------------------------!
               !----- Set the number of levels. -------------------------------------------!
               nzvp(iv) = iep_ng
               !---------------------------------------------------------------------------!

               !----- Adjust projection to GrADS. -----------------------------------------!
               call proj_rams_to_grads(vp(iv),ndim(iv),iep_nx(ng),iep_ny(ng),nzvp(iv)      &
                                      ,nxgrads(ng),nygrads(ng),rmi,iinf,jinf,rout          &
                                      ,routgrads,rlat,rlon,proj)
               !---------------------------------------------------------------------------!


               !----- Dump array to output file. ------------------------------------------!
               call ep_putvar(19,nxgrads(ng),nygrads(ng),nzvp(iv),nxa(ng),nxb(ng)          &
                             ,nya(ng),nyb(ng),1,nzvp(iv),routgrads,nrec)
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
      write(unit=20,fmt='(a)') 'undef -9.99e33'
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
            !      3-D variable, soil layers with no patch (obsolete).                     !
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






!==========================================================================================!
!==========================================================================================!
subroutine ep_getvar(cvar,rout,a,b,nx,ny,nz,ng,fn,cdname,cdunits,itype,npatch,nclouds,nzg  &
                    ,a2,rout2,a6,rout6)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   character(len=*)            , intent(in)    :: cvar
   character(len=*)            , intent(in)    :: fn
   character(len=*)            , intent(out)   :: cdname
   character(len=*)            , intent(out)   :: cdunits
   integer                     , intent(in)    :: nx
   integer                     , intent(in)    :: ny
   integer                     , intent(in)    :: nz
   integer                     , intent(in)    :: ng
   integer                     , intent(out)   :: itype
   integer                     , intent(in)    :: npatch
   integer                     , intent(in)    :: nzg
   integer                     , intent(in)    :: nclouds
   real   , dimension(*)       , intent(inout) :: a
   real   , dimension(*)       , intent(inout) :: b
   real   , dimension(*)       , intent(inout) :: a2
   real   , dimension(*)       , intent(inout) :: a6
   real   , dimension(*)       , intent(inout) :: rout
   real   , dimension(*)       , intent(inout) :: rout2
   real   , dimension(*)       , intent(inout) :: rout6
   !---------------------------------------------------------------------------------------!

   !----- Load the variable. --------------------------------------------------------------!
   call RAMS_varlib(cvar,nx,ny,nz,nzg,npatch,nclouds,ng,fn,cdname,cdunits,itype,a,b,a2,a6)   

   !----- Copy to the appropriate scratch. ------------------------------------------------!
   select case (itype)
   case (2)
      call atob(nx*ny,a,rout)
   case (3)
      call atob(nx*ny*nz,a,rout)
   case (6)
      call atob(nx*ny*nz*nclouds,a6,rout6)
   case (7)
      call atob(nx*ny*npatch,a,rout)
   case (8)
      call atob(nx*ny*nzg*npatch,a2,rout2)
   case (9)
      call atob(nx*ny*nclouds,a,rout)
   case (10)
      call atob(nx*ny*nzg,a,rout)
   end select
   return
end subroutine ep_getvar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ep_setdate(iyear1,imonth1,idate1,strtim,itrec)
   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   integer              , intent(in)  :: iyear1
   integer              , intent(in)  :: imonth1
   integer              , intent(in)  :: idate1
   real                 , intent(in)  :: strtim
   integer, dimension(6), intent(out) :: itrec
   !---------------------------------------------------------------------------------------!

   itrec(1) = iyear1
   itrec(2) = imonth1
   itrec(3) = idate1
   itrec(4) = int(mod(strtim,24.))
   itrec(5) = int(mod(strtim,1.)*60)
   itrec(6) = int(mod( (strtim) *3600.,60.))

   return
end subroutine ep_setdate
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine dumps the array to the output binary file (gra file).               !
!------------------------------------------------------------------------------------------!
subroutine ep_putvar(iunit,nxp,nyp,nzp,xa,xz,ya,yz,za,zz,array3d,irec)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: iunit
   integer                     , intent(in)    :: nxp
   integer                     , intent(in)    :: nyp
   integer                     , intent(in)    :: nzp
   integer                     , intent(in)    :: xa
   integer                     , intent(in)    :: xz
   integer                     , intent(in)    :: ya
   integer                     , intent(in)    :: yz
   integer                     , intent(in)    :: za
   integer                     , intent(in)    :: zz
   real, dimension(nxp,nyp,nzp), intent(in)    :: array3d
   integer                     , intent(inout) :: irec
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: x
   integer                                     :: y
   integer                                     :: z
   real, dimension(nxp,nyp)                    :: mat
   !---------------------------------------------------------------------------------------!
   do z=za,zz
      do y=1,nyp
         do x=1,nxp
            mat(x,y) = array3d(x,y,z)
         end do
      end do

      irec=irec+1
      write (unit=iunit,rec=irec) ((mat(x,y),x=xa,xz),y=ya,yz)
   end do

   return
end subroutine ep_putvar
!==========================================================================================!
!==========================================================================================!

!-------------------------------------------------------------------
!
Subroutine array_interpol(ng,nxg,nyg,nxr,nyr,rlat1,dlat, &
     rlon1,dlon,iinf,jinf,rmi,proj)
  use rpost_coms
  use brams_data
  use misc_coms, only : glong, glatg

  Dimension rmi(nxg,nyg,4),iinf(nxg,nyg),jinf(nxg,nyg)
  character(len=*) :: proj
  !
  if(trim(proj) == 'no') RETURN

  !       Construcao da matriz de interpolacao.
  !       Flag para pontos do grads fora do dominio do modelo
  undef=-9.99e33
  do i=1,nxg
     do j=1,nyg
        iinf(i,j)=1
        jinf(i,j)=1
        do l=1,4
           rmi(i,j,l)=undef
        enddo
     enddo
  enddo
  do i=1,nxg
     do j=1,nyg
        !       Encontra posicao do ponto de grade do GRADS na grade do RAMS
        !
        !        glatg(i)=-37.113
        !        glong(j)=-79.128
        !        xlat=glatg(i)
        !        xlon= glong(j)
        !        Call ll_xy(glatg(i),glong(j),polelat,polelon,x,y)
        !       call getops(pla,plo,xlat,xlon,polelat,polelon)
        !       call pstoxy(x,y,pla,plo,6376000.)
        !
        call ge_to_xy(polelat,polelon,glong(i),glatg(j),x,y)
        !
        !        print*,x,y,glong(i),glatg(j),polelat,polelon

        !       Elimina pontos fora:
        if(x.lt.xtn(1,ng).or.x.gt.xtn(nxr,ng)) go to 777
        if(y.lt.ytn(1,ng).or.y.gt.ytn(nyr,ng)) go to 777
        !        
        do ix=1,nxr
           if(x.le.xtn(ix,ng)) go to 555
        enddo
555     continue
        i1=ix-1
        i2=ix
        iinf(i,j)=i1         

        do iy=1,nyr
           if(y.le.ytn(iy,ng)) go to 666
        enddo
666     continue
        j1=iy-1
        j2=iy                
        jinf(i,j)=j1
        !        
        rmi(i,j,1)=(x-xtn(i1,ng))/deltaxn(ng)
        rmi(i,j,2)=1.-rmi(i,j,1)
        !         
        rmi(i,j,3)=(y-ytn(j1,ng))/deltayn(ng)
        rmi(i,j,4)=1.-rmi(i,j,3)
777     continue
        !
        !         print*,i,j,rmi(i,j,1),rmi(i,j,2),rmi(i,j,3),rmi(i,j,4)
        !
        !        
     enddo
  enddo
  return
end Subroutine array_interpol

!*************************************************************************

Subroutine proj_rams_to_grads(vp,n,nxr,nyr,nzz,nxg,nyg,     &
     rmi,iinf,jinf,                &
     rout,routgrads,rlat,rlon,proj)

  character*(*) proj
  character*10 vp
  Dimension rlat(nxr,nyr),rlon(nxr,nyr)
  Dimension rout(nxr,nyr,nzz),routgrads(nxg,nyg,nzz)
  Dimension rmi(nxg,nyg,4),iinf(nxg,nyg),jinf(nxg,nyg)


  if(trim(proj) == 'no') then
     if (nxg /= nxr .or. nyg /= nyr) then
        call abort_run  ('Projection with problems...','proj_rams_to_grads','rpost_main.f90')
     endif
     call rout_to_routgrads(nxr*nyr*nzz,rout,routgrads)
     return
  endif

  do i=1,nxg
     do j=1,nyg

        !
        !         print*,i,j,rmi(i,j,1),rmi(i,j,2),rmi(i,j,3),rmi(i,j,4)
        !
        !        
        r1= rmi(i,j,1)
        r2= rmi(i,j,2)
        r3= rmi(i,j,3)
        r4= rmi(i,j,4)
        i1= iinf(i,j)
        i2= i1+1
        j1= jinf(i,j)
        j2= j1+1


        do k=1,nzz
           rr1=   rout(i1,j1,k)*(1.-r1)+rout(i2,j1,k)*(1.-r2) 
           rr2=   rout(i1,j2,k)*(1.-r1)+rout(i2,j2,k)*(1.-r2) 
	   routgrads(i,j,k)=rr1*(1.-r3)+          rr2*(1.-r4)

           !	   print*,rr1,rr2,rout(i1,j1,k),routgrads(i,j,k)

           if(abs(routgrads(i,j,k)).gt.1.E+06)  &
                routgrads(i,j,k)=-9.99E+33
           !
           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           !   write(2,0998)
           !   write(2,0999) i,j,i1,j1,glatg(j),glong(i)
           !   write(2,1000) rlat(i1,j1),rlat(i2,j1),rlat(i1,j2),rlat(i2,j2)
           !   write(2,1001) rlon(i1,j1),rlon(i2,j1),rlon(i1,j2),rlon(i2,j2)
           !   write(2,1002) rout(i1,j1,k),rout(i2,j1,k),rout(i1,j2,k),&
           !		  rout(i2,j2,k), routgrads(i,j,k)
           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        enddo
     enddo
  enddo

0998 format(1x,'---------------------------------------------')
0999 format(1x,4i3,2f10.2)
1000 format(1x,4f10.2)
1001 format(1x,4f10.2)
1002 format(1x,4f10.2,f16.3)

  !xxxxxxxxxxxxxxxxxxxxxxxxx
  !      k=1
  !      do jj=1,nyr
  !      do ii=1,nxr
  !         write(10,'(2i3,3f8.1)')ii,jj,rlat(ii,jj),rlon(ii,jj)
  !     +             ,rout(ii,jj,k)
  !      enddo
  !      enddo
  !      do jj=1,nyg
  !      do ii=1,nxg
  !         write(11,'(2i3,3f8.1)')ii,jj,glatg(jj),glong(ii)
  !     +             ,routgrads(ii,jj,k)
  !      enddo
  !      enddo
  !xxxxxxxxxxxxxxxxxxxxxxxxx
  return
end Subroutine proj_rams_to_grads

!---------------------------------------------------------------------
!---------------------------------------------------------------------                  
subroutine ge_to_xy(polelat,polelon,xlon,xlat,x,y)

  parameter(rt=6367000.00)
  p=3.14159265360/180.00

  !       transformacao horizontal:
  b = 1.0+sin(p*xlat)*sin(p*polelat)+                 &
       cos(p*xlat)*cos(p*polelat)*cos(p*(xlon-polelon))

  f = 2.00*rt/b


  y = f*(cos(p*polelat)*sin(p*xlat) -                 &
       sin(p*polelat)*cos(p*xlat)*cos(p*(xlon-polelon)))

  x = f*(cos(p*xlat)*sin(p*(xlon - polelon)))

  return
end subroutine ge_to_xy

!---------------------------------------------------------------------

Subroutine geo_grid(nx,ny,rlat,rlon,dep_glon1,dep_glon2,  &
     dep_glat1,dep_glat2,		         &
     rlatmin,rlatmax,rlonmin,rlonmax,  	 &
     nxg,nyg,proj)
  use rpost_dims
  use misc_coms, only : glong, glatg
  Dimension rlat(nx,ny),rlon(nx,ny)
  character*(*) proj

  dep_glon1=rlon(1,1)
  dep_glon2=rlon(nx,1)
  do n=1,ny
     if(rlon(1,n).gt.dep_glon1)  dep_glon1=rlon(1,n)
     if(rlon(nx,n).lt.dep_glon2) dep_glon2=rlon(nx,n)
  enddo
  dep_glon2= (dep_glon2-dep_glon1)/(nx-1)


  dep_glat1=rlat(1,1)
  dep_glat2=rlat(1,ny)
  do n=1,nx
     if(rlat(n,1).gt.dep_glat1)  dep_glat1=rlat(n,1)
     if(rlat(n,ny).lt.dep_glat2) dep_glat2=rlat(n,ny)
  enddo
  dep_glat2= (dep_glat2-dep_glat1)/(ny-1)

  !10/08/98
  x=0
  xx=0
  do n=1,ny
     x=x+rlon(1,n)
     xx=xx+ (rlon(nx,n)-rlon(1,n))/(nx-1)
  enddo
  dep_glon1= x/ny
  dep_glon2=xx/ny


  x=0
  xx=0
  do n=1,nx
     x=x+rlat(n,1)
     xx=xx+ (rlat(n,ny)-rlat(n,1))/(ny-1)
  enddo
  dep_glat1= x/nx
  dep_glat2=xx/nx

  if(proj == 'no') then
     nxg=nx
     nyg=ny

  elseif (proj == 'yes') then

     !...... Grade para o GRADS:

     rlatmin=rlat(1,1)
     rlatmax=rlat(1,1)
     rlonmin=rlon(1,1)
     rlonmax=rlon(1,1)
     do i=1,nx
        do j=1,ny
           rlatmin=min(rlatmin,rlat(i,j))
           rlatmax=max(rlatmax,rlat(i,j))
           rlonmin=min(rlonmin,rlon(i,j))
           rlonmax=max(rlonmax,rlon(i,j))
        enddo
     enddo

     !...... Definicao da grade do GRADS 
     !
     ! Para testar dependencia com a resolucao da grade
     !
     !            dep_glon2=0.5*dep_glon2
     !            dep_glat2=0.5*dep_glat2

     nxg=int((rlonmax-rlonmin)/dep_glon2+0.5)+4
     nyg=int((rlatmax-rlatmin)/dep_glat2+0.5)+4
     rlon1=rlonmin-(nxg-nx-1)*dep_glon2
     rlat1=rlatmin-(nyg-ny-1)*dep_glat2

     !            rlon1=rlonmin-2*dep_glon2
     !            rlat1=rlatmin-2*dep_glat2

     rlon1=rlonmin-dep_glon2
     rlat1=rlatmin-dep_glat2
     dep_glat1= rlat1
     dep_glon1= rlon1

     !            print*,rlonmin,rlonmax,rlatmin,rlatmax
     !            print*,nxg,nyg,dep_glon2,dep_glat2,dep_glat1,dep_glon1
     !            stop
  else 
     call abort_run   ('Invalid value for iproj: '//trim(proj)//'...' &
                      ,'geo_grid','rpost_main.f90')
  endif




  !	Define	grade do GRADS

  do i=1,nxg
     glong(i)=dep_glon1+float(i-1)*dep_glon2
     !        print*,' i lon=',i,glong(i)
  enddo

  do j=1,nyg
     glatg(j)=dep_glat1+float(j-1)*dep_glat2
     !        print*,' j lat=',j,glatg(j)
  enddo

  return
end Subroutine geo_grid

!---------------------------------------------------------------------                  
subroutine rout_to_routgrads(nxyz,rinp,rout)
  dimension rinp(nxyz),rout(nxyz)
  do i=1,nxyz
     rout(i)=rinp(i)
  enddo
  return
end subroutine rout_to_routgrads

!---------------------------------------------------------------------               

