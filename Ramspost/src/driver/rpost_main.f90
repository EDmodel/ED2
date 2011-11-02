!==========================================================================================!
!==========================================================================================!
! -                                                           -
! - RAMSPOST - RAMS Post Processor for GrADS.                 -
! -                                                           -
! -------------------------------------------------------------
! - Adapted for RAMS 4.3 by Saulo R Freitas (SP/14/05/1998)   -
! - Adapted for RAMS 6.0 by Saulo R Freitas (SP/06/03/2005)   -
! -------------------------------------------------------------
program ramspost

  ! -----------------------
  ! -   BASIC PARAMETERS  -
  ! -----------------------

  use rpost_coms
  use rpost_dims, only  : str_len            ! ! intent(in)
  use brams_data
  use leaf_coms , only  : sfclyr_init_params ! ! sub-routine
  character(len=str_len), dimension(maxfiles) :: fln
  character(len=str_len)                      :: inp
  character(len=str_len)                      :: fprefix
  character(len=str_len)                      :: gprefix
  character(len=str_len)                      :: cfln
  character(len=40), dimension(200)           :: vpln
  character(len=40)                           :: cdum1
  character(len=20), dimension(200)           :: vp
  character(len=20), dimension(200)           :: vpun
  character(len=20)                           :: cdum2
  character(len=20)                           :: proj
  character(len=20)                           :: anl2gra
  character(len=1)                            :: cgrid
  character(len=2)                            :: ccgrid
  character(len=2)                            :: patchnumber
  character(len=2)                            :: cldnumber
  character(len=3), dimension(12)             :: cmo
  character*15 chdate,chstep,xchstep
  integer nvp,nfiles,nzvp(200),nrec,ipresslev,iplevs(nplmax)
  integer inplevs,zlevmax(maxgrds),ndim(200),iproj,ianl2gra,icld
  integer fim_inp
  integer hunit

  real a(nxpmax,nypmax,nzpmax),b(nxpmax,nypmax,nzpmax),        &
       rout(nxpmax,nypmax,nzpmax),                             &
       zplev(nxpmax,nypmax,nplmax),                            &
       mytopo(nxpmax,nypmax),mypi(nxpmax,nypmax,nzpmax),       &
       rlat(nxpmax,nypmax),rlon(nxpmax,nypmax),	               &
       lati(maxgrds),latf(maxgrds),loni(maxgrds),lonf(maxgrds),&
       a2(nxpmax,nypmax,nzgmax,maxpatch),                      &
       rout2(nxpmax,nypmax,nzgmax,maxpatch),                   &
       a6(nxpmax,nypmax,nzpmax,maxclouds),                     &
       rout6(nxpmax,nypmax,nzpmax,maxclouds)

  integer nxgrads(maxgrds),nygrads(maxgrds),            &
       iinf (maxgx,maxgy), jinf (maxgx,maxgy),       &
       nxa(maxgrds),nxb(maxgrds),nya(maxgrds),nyb(maxgrds) 
  real rmi(maxgx,maxgy,4),routgrads(maxgx,maxgy,nzpmax)
  !  ---

  namelist/rp_input/ fprefix,nvp,vp,gprefix,anl2gra,proj,lati,latf,  &
       loni,lonf,zlevmax,ipresslev,inplevs,iplevs

  data cmo/'jan','feb','mar','apr','may','jun','jul','aug','sep', &
       'oct','nov','dec'/	

  dimension dep_zlev(nzpmax,maxgrds),iep_nx(maxgrds),       &
       iep_ny(maxgrds),iep_nz(maxgrds),dep_glat(2,maxgrds), dep_glon(2,maxgrds)

  character(len=str_len), dimension(maxfiles) ::  wfln

  ! -----------------------------
  ! -   INITIALIZING ROUTINES   -
  ! -----------------------------
  print*,'############################################'
  print*,' RamsPost - GrADS Visualization for RAMS    '
  print*,'############################################'

  call  getarg(1,inp)
  fim_inp=index(inp,' ')
  if (fim_inp<3) inp='ramspost.inp'

  print *, ' '
  print *, 'Opening '//trim(inp)//' file'

  open(5,file=trim(inp),status='old')
  read(5,rp_input)
  cgrid='0'
  iproj   =lastchar(proj)
  ianl2gra=lastchar(anl2gra)

  ! --- frequencia com as analises serao escrita     
  nstep = 1    

  nrec=0
  ic=lastchar(gprefix)
  call RAMS_anal_init(nfiles,fln,fprefix,        &
       dep_zlev,iep_nx,iep_ny,iep_nz,iep_ng,iep_np,iep_nc,   &
       iep_ngrids)
  chdate='00:00z00mmm1900'
  call RAMS_get_time_init(1,iyear,imonth,idate,ihour,imin)
  call RAMS_get_time_step(iistep,hunit,nfiles)

!  print*,iyear,imonth,idate,ihour,imin      
  write(chdate(1:2),'(i2.2)') ihour
  write(chdate(4:5),'(i2.2)') imin
  write(chdate(7:8),'(i2.2)') idate
  write(chdate(12:15),'(i4.2)') iyear
  chdate(9:11)=cmo(imonth)(1:3)

  if(hunit.eq.1) chstep='          sec'
  if(hunit.eq.2) chstep='          mn'
  if(hunit.eq.3) chstep='          hr'
  write(chstep(8:10),'(i3)') iistep

  ! -----------------
  ! -   GRID LOOP   -
  ! -----------------


  do ng=1,iep_ngrids
     write (unit=*,fmt='(a)')       ' '
     write (unit=*,fmt='(a)')       '========================================================='
     write (unit=*,fmt='(a,1x,i5)'),' + Writing Grid ',ng
     !.................
     nnvp=nvp
     iv=1
     nfn=1
     cfln=fln(nfn)
     ip=lastchar(cfln)-9
     !.................
     !   rlat and rlon = lat and lon of "thermodynamic points" of RAMS model.
     !
     write(unit=*,fmt='(a)') ' '
     write(unit=*,fmt='(2(a,1x))') '     * Variable:  ','lat'

     call ep_getvar('lat',                    &
          rlat,a,b,iep_nx(ng),iep_ny(ng),         &
          1,ng,cfln(1:ip),vpln(iv),		     &
          vpun(iv),n,iep_np,iep_nc,iep_ng,a2,rout2,a6,rout6)
     write(unit=*,fmt='(a,1x,i5)') '       # Output variable type:  ',n

     write(unit=*,fmt='(a)') ' '
     write(unit=*,fmt='(2(a,1x))') '     * Variable:  ','lon'

     call ep_getvar('lon',                    &
          rlon,a,b,iep_nx(ng),iep_ny(ng),         &
          1,ng,cfln(1:ip),vpln(iv),  	     &
          vpun(iv),n,iep_np,iep_nc,iep_ng,a2,rout2,a6,rout6)
     write(unit=*,fmt='(a,1x,i5)') '       # Output variable type:  ',n


     !.................
     call geo_grid(iep_nx(ng),iep_ny(ng),rlat,rlon,  &
          dep_glon(1,ng),dep_glon(2,ng),     &
          dep_glat(1,ng),dep_glat(2,ng),     &
          rlatmin,rlatmax,rlonmin,rlonmax,   &
          nxgrads(ng),nygrads(ng),           &
          proj(1:iproj))
     !            print*,nxgrads(ng),nygrads(ng),iep_nx(ng),iep_ny(ng)
     !            print*, proj(1:iproj),iproj
     !.................            
     !
     Call Matriz_interp(ng,nxgrads(ng),nygrads(ng),   &
          iep_nx(ng),iep_ny(ng),	     &
          dep_glat(1,ng),dep_glat(2,ng),     &
          dep_glon(1,ng),dep_glon(2,ng),     &
          iinf,jinf,rmi,                     &
          proj(1:iproj))

     !_................
     Call define_lim(ng,nxgrads(ng),nygrads(ng),          &
          dep_glat(1,ng),dep_glat(2,ng),	 &
          dep_glon(1,ng),dep_glon(2,ng),	 &
          lati(ng),latf(ng),loni(ng),lonf(ng),  &
          nxa(ng),nxb(ng),nya(ng),nyb(ng),proj(1:iproj), &
          iep_nx(ng),iep_ny(ng),rlat,rlon)
     !            print*,nxgrads(ng),nygrads(ng),iep_nx(ng),iep_ny(ng) &
     !	           ,iep_nz(ng),iv,vpln(iv), nxa(ng),nxb(ng),nya(ng),nyb(ng),&
     !     	            vpun(iv),n
     !     ----------------------
     !     - WRITE GRADS BINARY -
     !     ----------------------
     !     ----------------------
     write(cgrid,'(i1)')ng
     if(anl2gra(1:ianl2gra) .ne. 'ONE' .and. &
        anl2gra(1:ianl2gra) .ne. 'one' ) then
       open(19,file=gprefix(1:ic)//'_g'//cgrid//'.gra',         &
            form='unformatted',access='direct',status='unknown',  &
            recl=4*(nxb(ng)-nxa(ng)+1)*(nyb(ng)-nya(ng)+1))	  
       nrec=0
     endif

     do nfn=1,nfiles,nstep
        write(unit=*,fmt='(a)')        ' '
        write(unit=*,fmt='(a,1x,i5)')  '   - Timestep: ',nfn

        cfln=fln(nfn)
        ip=lastchar(cfln)-9
        !.................
!.................
! open arquivos individuais
     !     ----------------------
     !     - WRITE GRADS BINARY -
     !     ----------------------
       if(anl2gra(1:ianl2gra) .eq. 'ONE' .or. &
          anl2gra(1:ianl2gra) .eq. 'one' ) then
           !print*,ftimes(nfn),ifdates(nfn),iftimes(nfn),startutc,httop 
           call date1(ifdates(nfn),iyear,imon,idate)
           call makefnam(wfln(nfn),gprefix(1:ic)//' ',0.,iyear,imon,idate,  &
                        iftimes(nfn),'A','g'//cgrid,'gra')
           !print*,iyear,imon,idate,iftimes(nfn),wfln(nfn)
           iunit=19
           open(iunit,file=wfln(nfn),form='unformatted',access='direct'         &
               ,status='unknown',recl=4*(nxb(ng)-nxa(ng)+1)*(nyb(ng)-nya(ng)+1))	  
           nrec=0
	endif
!.................

        if(ipresslev.gt.0) then
           write(unit=*,fmt='(a)') ' '
           write(unit=*,fmt='(2(a,1x))') '     * Variable:  ','topo'

           call ep_getvar('topo',mytopo,a,b,iep_nx(ng),iep_ny(ng), &
                1,ng,cfln(1:ip),vpln(iv),		    &
                vpun(iv),n,iep_np,iep_nc,iep_ng,a2,rout2,a6,rout6)
           write(unit=*,fmt='(a,1x,i5)') '       # Output variable type:  ',n

           write(unit=*,fmt='(a)') ' '
           write(unit=*,fmt='(2(a,1x))') '     * Variable:  ','pi'

           call ep_getvar('pi',mypi,a,b,iep_nx(ng),iep_ny(ng),     &
                iep_nz(ng),ng,cfln(1:ip),vpln(iv),       &
                vpun(iv),n,iep_np,iep_nc,iep_ng,a2,rout2,a6,rout6)
           write(unit=*,fmt='(a,1x,i5)') '       # Output variable type:  ',n

        endif

      DO iv=1,nvp
           write(unit=*,fmt='(a)') ' '
           write(unit=*,fmt='(2(a,1x))') '     * Variable:  ',vp(iv)
           !.................
           call ep_getvar(vp(iv),                               &
                rout,a,b,iep_nx(ng),iep_ny(ng),        &
                iep_nz(ng),ng,cfln(1:ip),vpln(iv),	 &
                vpun(iv),n,iep_np,iep_nc,iep_ng,a2,rout2,a6,rout6)
           write(unit=*,fmt='(a,1x,i5)') '       # Output variable type:  ',n

           ndim(iv)=n
           !................
           !.....            
           IF(ndim(iv).eq.8) then
              nzvp(iv)=iep_ng
              do ipatch=1,iep_np
                 Call S4d_to_3d(iep_nx(ng),iep_ny(ng),nzvp(iv),iep_np, &
                      ipatch,rout,rout2)
                 Call proj_rams_to_grads(vp(iv),ndim(iv),        &
                      iep_nx(ng),iep_ny(ng),nzvp(iv),  &
                      nxgrads(ng),nygrads(ng),	      &
                      rmi,iinf,jinf,		      &
                      rout,routgrads,rlat,rlon,proj(1:iproj))
                 Call ep_putvar(routgrads,a,nxgrads(ng),nygrads(ng),   &
                      nxa(ng),nxb(ng),nya(ng),nyb(ng),	    &
                      nzvp(iv),nrec,1,iep_ng)
                 if(nfn.eq.1) nnvp=nnvp+1
              enddo
              if(nfn.eq.1) nnvp=nnvp-1
              !.....
              !
 	    ELSEIF(ndim(iv).eq.10) then
!!..
              nzvp(iv)=iep_ng
	      Call proj_rams_to_grads(vp(iv),ndim(iv),       &
			     iep_nx(ng),iep_ny(ng),nzvp(iv),  &
			     nxgrads(ng),nygrads(ng),	      &
			     rmi,iinf,jinf,		      &
			     rout,routgrads,rlat,rlon,proj(1:iproj))

	      Call ep_putvar(routgrads,a,nxgrads(ng),nygrads(ng),  &
			     nxa(ng),nxb(ng),nya(ng),nyb(ng),	    &
			     nzvp(iv),nrec,1,nzvp(iv))
         ELSEIF(ndim(iv).eq.6) then
              nzvp(iv)=iep_nz(ng)
              do icld=1,iep_nc
                 Call S4d_to_3d(iep_nx(ng),iep_ny(ng),nzvp(iv),iep_nc, &
                      icld,rout,rout6)
                 Call proj_rams_to_grads(vp(iv),ndim(iv),        &
                      iep_nx(ng),iep_ny(ng),nzvp(iv),  &
                      nxgrads(ng),nygrads(ng),	      &
                      rmi,iinf,jinf,		      &
                      rout,routgrads,rlat,rlon,proj(1:iproj))
                 Call ep_putvar(routgrads,a,nxgrads(ng),nygrads(ng),   &
                      nxa(ng),nxb(ng),nya(ng),nyb(ng),	    &
                      nzvp(iv),nrec,1,nzvp(iv))
                 if(nfn.eq.1) nnvp=nnvp+1
              enddo
              if(nfn.eq.1) nnvp=nnvp-1
              !.....
              !

          ELSEIF(ndim(iv).eq.2)  THEN
              nzvp(iv)=1
              Call proj_rams_to_grads(vp(iv),ndim(iv),     &
                   iep_nx(ng),iep_ny(ng),nzvp(iv),&
                   nxgrads(ng),nygrads(ng),	    &
                   rmi,iinf,jinf,		    &
                   rout,routgrads,rlat,rlon,proj(1:iproj))
              Call ep_putvar(routgrads,a,nxgrads(ng),nygrads(ng),  &
                   nxa(ng),nxb(ng),nya(ng),nyb(ng),	    &
                   nzvp(iv),nrec,1,1)
              !!.....
              !
           ELSEIF(ndim(iv).eq.7)  THEN
              nzvp(iv)=iep_np
	      Call proj_rams_to_grads(vp(iv),ndim(iv),      &
                   iep_nx(ng),iep_ny(ng),nzvp(iv),  &
                   nxgrads(ng),nygrads(ng),	    &
                   rmi,iinf,jinf,		    &
                   rout,routgrads,rlat,rlon,proj(1:iproj))
	      do ipatch=1,iep_np
                 Call ep_putvar(routgrads,a,nxgrads(ng),nygrads(ng), &
                      nxa(ng),nxb(ng),nya(ng),nyb(ng),	   &
                      nzvp(iv),nrec,ipatch,ipatch)
                 if(nfn.eq.1) nnvp=nnvp+1
              enddo
              if(nfn.eq.1) nnvp=nnvp-1
           !MLO - For cloud type (x,y,cloud) dependent variable
	   ELSEIF(ndim(iv).eq.9)  THEN
	     nzvp(iv)=iep_nc
	     Call proj_rams_to_grads(vp(iv),ndim(iv),      &
			   iep_nx(ng),iep_ny(ng),nzvp(iv),  &
			   nxgrads(ng),nygrads(ng),	    &
			   rmi,iinf,jinf,		    &
			   rout,routgrads,rlat,rlon,proj(1:iproj))
	     do icld=1,iep_nc
	       Call ep_putvar(routgrads,a,nxgrads(ng),nygrads(ng), &
			       nxa(ng),nxb(ng),nya(ng),nyb(ng),	   &
		 	       nzvp(iv),nrec,icld,icld)
	       if(nfn.eq.1) nnvp=nnvp+1
	     enddo
	     if(nfn.eq.1) nnvp=nnvp-1
!
   !!.....
              !
           ELSEIF(ndim(iv).eq.5)  THEN
              !	     nzvp(iv)=iep_ng
              !	       Call proj_rams_to_grads(vp(iv),ndim(iv),       &
              !			     iep_nx(ng),iep_ny(ng),nzvp(iv),  &
              !			     nxgrads(ng),nygrads(ng),	      &
              !			     rmi,iinf,jinf,		      &
              !			     rout,routgrads,rlat,rlon,proj(1:iproj))
              !
              !		Call ep_putvar(routgrads,a,nxgrads(ng),nygrads(ng),  &
              !			     nxa(ng),nxb(ng),nya(ng),nyb(ng),	    &
              !			      nzvp(iv),nrec,1,iep_ng)
              !!.....
              !
	   ELSEIF(ndim(iv).eq.3)  THEN
              nzvp(iv)=iep_nz(ng)
              !..
              if(ipresslev.eq.1) then

                 inplevsef=inplevs
                 Call  ptransvar(rout,iep_nx(ng),iep_ny(ng),nzvp(iv),  &
                      inplevs,iplevs,mypi,dep_zlev(1,ng),zplev,mytopo)

                 Call proj_rams_to_grads(vp(iv),ndim(iv),	      &
                      iep_nx(ng),iep_ny(ng),nzvp(iv),   &
                      nxgrads(ng),nygrads(ng),	     &
                      rmi,iinf,jinf,		     &
                      rout,routgrads,rlat,rlon,proj(1:iproj))

                 Call ep_putvar(routgrads,a,nxgrads(ng),nygrads(ng), &
                      nxa(ng),nxb(ng),nya(ng),nyb(ng),	 &
                      inplevsef,nrec,1,inplevsef)

              elseif(ipresslev.eq.2) then

                 inplevsef=inplevs
                 Call  ctransvar(iep_nx(ng),iep_ny(ng),iep_nz(ng),rout &
                      ,mytopo,inplevs,iplevs,ztn(1,ng),zmn(nnzp(1)-1,1))

                 Call proj_rams_to_grads(vp(iv),ndim(iv),     &
                      iep_nx(ng),iep_ny(ng),nzvp(iv), &
                      nxgrads(ng),nygrads(ng),	   &
                      rmi,iinf,jinf,		   &
                      rout,routgrads,rlat,rlon,proj(1:iproj))

                 Call ep_putvar(routgrads,a,nxgrads(ng),nygrads(ng),  &
                      nxa(ng),nxb(ng),nya(ng),nyb(ng),	 &
                      inplevsef,nrec,1,inplevsef)

              elseif(ipresslev.eq.3) then

                 inplevsef=inplevs
                 Call  select_sigmaz(iep_nx(ng),iep_ny(ng),iep_nz(ng),rout &
                      ,inplevs,iplevs)

                 Call proj_rams_to_grads(vp(iv),ndim(iv),     &
                      iep_nx(ng),iep_ny(ng),nzvp(iv), &
                      nxgrads(ng),nygrads(ng),	   &
                      rmi,iinf,jinf,		   &
                      rout,routgrads,rlat,rlon,proj(1:iproj))

                 Call ep_putvar(routgrads,a,nxgrads(ng),nygrads(ng),  &
                      nxa(ng),nxb(ng),nya(ng),nyb(ng),	 &
                      inplevsef,nrec,1,inplevsef)
                 !!..
              else
                 !!..
                 Call proj_rams_to_grads(vp(iv),ndim(iv),       &
                      iep_nx(ng),iep_ny(ng),nzvp(iv),  &
                      nxgrads(ng),nygrads(ng),	      &
                      rmi,iinf,jinf,		      &
                      rout,routgrads,rlat,rlon,proj(1:iproj))
                 if(zlevmax(ng).ge.iep_nz(ng)) zlevmax(ng)=iep_nz(ng)-1

                 Call ep_putvar(routgrads,a,nxgrads(ng),nygrads(ng),  &
                      nxa(ng),nxb(ng),nya(ng),nyb(ng),	    &
                      nzvp(iv),nrec,2,zlevmax(ng)+1)
                 !!..
              endif
           ELSE
             ! Removing one variable from the counting...
             nnvp=nnvp-1
              !!
	   Endif
           !!................
           !
           !
           !

        ENDDO
        !.................

        if(anl2gra(1:ianl2gra).eq.'ONE'.or.anl2gra(1:ianl2gra).eq.'one') close(19)
       enddo ! enddo do NFILES

      if(anl2gra(1:ianl2gra).ne.'ONE'.and.anl2gra(1:ianl2gra).ne.'one')  close(19)


     !.................
     !     -----------------------
     !     - WRITE GRADS CONTROL -
     !     -----------------------

     write(cgrid,'(i1)')ng
     iunit = 20
     
     do nfn=1,nfiles,nstep

     iuniti=iunit
     iunitf=iunit

! case 1 : all analysis at only one grads file
!----
      if(anl2gra(1:ianl2gra).ne.'ONE'.and.anl2gra(1:ianl2gra).ne.'one' .and. &
         nfn == 1 ) then 
           open(iunit,file=gprefix(1:ic)//'_g'//cgrid//'.ctl', &
                status='unknown')
           write(iunit,2001) '^'//gprefix(1:ic)//'_g'//cgrid//'.gra'
! case 2 : one analysis at one grads file
      elseif(anl2gra(1:ianl2gra).eq.'ONE'.or.anl2gra(1:ianl2gra).eq.'one') then
	   call date1(ifdates(nfn),iyear,imon,idate)	   
           call makefnam(wfln(nfn),gprefix(1:ic)//' ',0.,iyear,imon,idate,  &
                        iftimes(nfn),'A','g'//cgrid,'ctl')      
           open(iunit,file=wfln(nfn),status='unknown')
           write(iunit,2001) '^'//wfln(nfn)(1:lastchar(wfln(nfn))-3)//'gra'
! case 2 with template
           if(nfn==1) then
	     call date1(ifdates(nfn),iyear,imon,idate)	   
             call makefnam(wfln(nfn),gprefix(1:ic)//'-template'//' ',0.,iyear,imon,idate,  &
                        iftimes(nfn),'A','g'//cgrid,'ctl')      
             open(iunit+1,file=wfln(nfn),status='unknown')
!valido somente para hora cheia  --------------------------------------vvvv-
             write(iunit+1,2001) '^'//gprefix(1:ic)//'-A-'//'%y4-%m2-%d2-%h20000-'//'g'//cgrid//'.gra'
             write(iunit+1,2002) 'options template'
	     iunitf=iunit+1
	   endif   
	   
      endif
!----

      do iunit=iuniti,iunitf
           
       write(iunit,2002) 'undef -9.99e33'
       write(iunit,2002) 'title RAMS 4.2 Output'
       write(iunit,2003) nxb(ng)-nxa(ng)+1,(dep_glon(i,ng),i=1,2)
       write(iunit,2004) nyb(ng)-nya(ng)+1,(dep_glat(i,ng),i=1,2)
       if(ipresslev.gt.0.and.ipresslev.le.2) then
         write(iunit,2005) inplevs,(iplevs(i)*1.0,i=1,inplevs)
       elseif(ipresslev.eq.3) then
         write(iunit,2005) inplevs,(dep_zlev(iplevs(i),ng),i=1,inplevs)
       else
        if(zlevmax(ng)+1.lt.15) then
           write(iunit,2005) zlevmax(ng), &
                (dep_zlev(i,ng),i=2,zlevmax(ng)+1)
        else
           write(iunit,2005) zlevmax(ng),(dep_zlev(i,ng),i=2,15)
           write(iunit,2055) (dep_zlev(i,ng),i=16,zlevmax(ng)+1)
        endif
       endif
! case 1
       if(anl2gra(1:ianl2gra).ne.'ONE'.and.anl2gra(1:ianl2gra).ne.'one' .and. &
         nfn == 1 ) then 
         write(iunit,2006) nfiles,chdate,chstep
! case 2
       elseif(anl2gra(1:ianl2gra).eq.'ONE'.or.anl2gra(1:ianl2gra).eq.'one') then

        call RAMS_get_time_init(nfn,iyear,imonth,idate,ihour,imin)
        write(chdate(1:2),'(i2.2)') ihour
    	write(chdate(4:5),'(i2.2)') imin
    	write(chdate(7:8),'(i2.2)') idate
    	write(chdate(12:15),'(i4.2)') iyear
    	chdate(9:11)=cmo(imonth)(1:3)

       if(iunitf== iuniti)   write(iunit,2006) 1,chdate,chstep
       if(iunitf== iuniti+1) write(iunit,2006) nfiles,chdate,chstep ! para template
       
       endif
!----
       write(iunit,2007) nnvp
       do i=1,nvp

        if(ipresslev.gt.0.and.nzvp(i).eq.iep_nz(ng)) then
           write(iunit,2008) vp(i),inplevs,vpln(i),vpun(i)
        else
           if(ndim(i).eq.3) &
                write(iunit,2008) vp(i),zlevmax(ng),vpln(i),vpun(i)

           if(ndim(i).eq.2) &
                write(iunit,2008) vp(i),0,vpln(i),vpun(i)

           if(ndim(i).eq.6) then
              il =lastchar(vp(i))
              il2=lastchar(vpln(i))
              do icld=1,iep_nc
                 write(cldnumber,'(i2.2)')icld
                 cdum2=vp(i)(1:il)//cldnumber
                 cdum1=vpln(i)(1:il2)//': Cloud # '//cldnumber
                 write(iunit,2008) cdum2,nzvp(i),cdum1,vpun(i)
              enddo
           endif

           if(ndim(i).eq.7) then
              il =lastchar(vp(i))
              il2=lastchar(vpln(i))
              do ipatch=1,iep_np
                 write(patchnumber,'(i2.2)')ipatch
                 cdum2=vp(i)(1:il)//patchnumber
                 cdum1=vpln(i)(1:il2)//': patch # '//patchnumber
                 write(iunit,2008) cdum2,0,cdum1,vpun(i)
              enddo
           endif
           if(ndim(i).eq.8) then
              il =lastchar(vp(i))
              il2=lastchar(vpln(i))
              do ipatch=1,iep_np
                 write(patchnumber,'(i2.2)')ipatch
                 cdum2=vp(i)(1:il)//patchnumber
                 cdum1=vpln(i)(1:il2)//': patch # '//patchnumber
                 write(iunit,2008) cdum2,nzvp(i),cdum1,vpun(i)
              enddo
           endif
              if(ndim(i).eq.9) then
                il =lastchar(vp(i))
                il2=lastchar(vpln(i))
                do icld=1,iep_nc
                 write(cldnumber,'(i2.2)')icld
                 cdum2=vp(i)(1:il)//cldnumber
                 cdum1=vpln(i)(1:il2)//': Cloud # '//cldnumber
                 write(20,2008) cdum2,0,cdum1,vpun(i)
                enddo
              endif

           if(ndim(i).eq.5) &
                write(iunit,2008) vp(i),nzvp(i),vpln(i),vpun(i)

              if(ndim(i).eq.10) &
                write(20,2008) vp(i),nzvp(i),vpln(i),vpun(i)

         endif
       enddo
       write(iunit,2002) 'endvars'
       close(iunit)

      enddo ! enddo nas unidades de escrita

!     endif
  
     if(anl2gra(1:ianl2gra).ne.'ONE'.and.anl2gra(1:ianl2gra).ne.'one' .and. &
        nfn == 1 ) exit
   enddo ! enddo do NFILES
   write (unit=*,fmt='(a)')       '========================================================='

  enddo  !enddo do NGRIDS

2001 format('dset ',a)
2002 format(a)
2003 format('xdef ',i4,' linear ',2f15.3)
2004 format('ydef ',i4,' linear ',2f15.3)
2005 format('zdef ',i4,' levels ',60f10.1)
2006 format('tdef ',i4,' linear ',2a15)
2007 format('vars ',i4)
2008 format(a,1x,i4,1x,' 99    - RAMS : ',1x,a,'[',1x,a8,1x,']')
2055 format(60f7.0)
close (iunit+1)

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
      call atob(nx*ny*nzg*nclouds,a6,rout6)
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
  real time

  integer itrec(6)
  itrec(1)=iyear1
  itrec(2)=imonth1
  itrec(3)=idate1
  itrec(4)=int(mod(strtim,24.))
  itrec(5)=int(mod(strtim,1.)*60)
  itrec(6)=int(mod( (strtim) *3600.,60.))

  !      print*,'---------------------------------------'
  !      print*,itrec(1),itrec(2),itrec(3),itrec(4),itrec(5),
  !     +       itrec(6)
  !      print*,'---------------------------------------'

  return
end subroutine ep_setdate

!*****************************************************************************

! --------------------------------------------------------
! -   SUBROUTINE EP_PUTVAR : WRITE ARRAY TO GRADS FILE   -
! --------------------------------------------------------

subroutine ep_putvar(rout,a,nx,ny,nxa,nxb,nya,nyb,   &
                     nz,nrec,istartz,iendz)
  dimension a(nx,ny),rout(nx,ny,nz)
  integer istartz,iendz
  ! 
  do k=istartz,iendz
     do j=1,ny
        do i=1,nx
           a(i,j)=rout(i,j,k)
           !cc
           !            print*,'PUT VAR=',i,j,k,a(i,j)
           !cc
        enddo
     enddo
     nrec=nrec+1
     write (19,rec=nrec) ((a(i,j),i=nxa,nxb),j=nya,nyb)
     !       write(19,rec=nrec) a
  enddo

  !      k=1
  !        write(18,'(59f10.3)')((rout(ii,jj,k),ii=1,nx),jj=1,ny)
  !        
  return
end subroutine ep_putvar

!-------------------------------------------------------------------
!
Subroutine Matriz_interp(ng,nxg,nyg,nxr,nyr,rlat1,dlat, &
     rlon1,dlon,iinf,jinf,rmi,proj)
  use rpost_coms
  use brams_data
  use misc_coms, only : glong, glatg

  Dimension rmi(nxg,nyg,4),iinf(nxg,nyg),jinf(nxg,nyg)
  character(len=*) :: proj
  !
  if(proj.ne.'YES'.AND.proj.ne.'yes') RETURN

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
end Subroutine Matriz_interp

!*************************************************************************

Subroutine proj_rams_to_grads(vp,n,nxr,nyr,nzz,nxg,nyg,     &
     rmi,iinf,jinf,                &
     rout,routgrads,rlat,rlon,proj)

  character*(*) proj
  character*10 vp
  Dimension rlat(nxr,nyr),rlon(nxr,nyr)
  Dimension rout(nxr,nyr,nzz),routgrads(nxg,nyg,nzz)
  Dimension rmi(nxg,nyg,4),iinf(nxg,nyg),jinf(nxg,nyg)


  if(proj.ne.'YES'.AND.proj.ne.'yes') then
     if(nxg.ne.nxr.AND.nyg.ne.nyr) then
        print*,'Projection with problems nxr nxg ...'
        stop
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

  if(proj.ne.'YES'.and.proj.ne.'yes') then
     nxg=nx
     nyg=ny

  else

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

