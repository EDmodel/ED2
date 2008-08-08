subroutine micro_opt()

  use mem_basic, only : &
       basic_g            ! INTENT(INOUT)

  use mem_micro, only:  &
       micro_g            ! INTENT(INOUT)

  ! ALF - For optimization
  use mem_micro_opt, only:   &
       indexing,             & ! Subroutine
       ij_last,              & ! INTENT(IN)
       indextile,            & ! INTENT(OUT)
       step_limit,           & ! INTENT(IN)
       ij_total,             & ! INTENT(OUT)
       ij_final,             & ! INTENT(OUT)
       alloc_micro_tile,     & ! Subroutine
       dealloc_micro_tile,   & ! Subroutine
       init_copy_mic_block,  & ! Subroutine
       final_copy_mic_block, & ! Subroutine
       copyback_opt,         & ! Subroutine
       range_check_opt,      & ! Subroutine
       mcphys_opt              ! Subroutine 

  use mem_grid, only:   &
       ngrids,          & ! INTENT(IN)
       ngrid,           & ! INTENT(IN)
       zm,              & ! INTENT(IN)
       dzt,             & ! INTENT(IN)
       dtlt,            & ! INTENT(IN)
       maxnzp,          & ! INTENT(IN)
       time,            & ! INTENT(IN)
       zt,              & ! INTENT(IN)
       if_adap,         & ! INTENT(IN)
       grid_g             ! INTENT(IN)

  use mem_radiate, only:&
       radiate_g          ! INTENT(INOUT)

  use node_mod, only :  &
       mmzp,            & ! INTENT(IN)
       mzp,             & ! INTENT(IN)
       mxp,             & ! INTENT(IN)
       myp,             & ! INTENT(IN)
       ja,              & ! INTENT(IN)
       jz,              & ! INTENT(IN)
       ia,              & ! INTENT(IN)
       iz,              & ! INTENT(IN)
       mynum              ! INTENT(IN)

  use micphys, only:    &
       maxgrds,         & ! INTENT(IN)
       level,           & ! INTENT(IN)
       nhcat,           & ! INTENT(IN)
       ch3,             & ! INTENT(OUT)
       pwvt,            & ! INTENT(IN)
       pwmasi,          & ! INTENT(IN)
       cdp1,            & ! INTENT(OUT)
       pwvtmasi,        & ! INTENT(OUT)
       ch2,             & ! INTENT(OUT)
       dispemb1,        & ! INTENT(IN)
       dispemb0           ! INTENT(IN)

  implicit none

  ! Local Variables:

  integer :: nembfall,maxkfall,ngr,lhcat,i,j,k
  integer, dimension(10)  :: k1,k2,k3
  integer, save :: ncall = 0
  integer, save, dimension(15)  :: ncall2g = 0
  real :: dtlti
  type pcp_tab_type
     real, pointer, dimension(:,:,:,:) :: pcpfillc,pcpfillr
     real, pointer, dimension(:,:,:)   :: sfcpcp
  end type pcp_tab_type
  type (pcp_tab_type), save :: pcp_tab(maxgrds)

  integer :: ij, step, ijm

  if (level /= 3) return

  nembfall = 20
  maxkfall = 4

  if(ncall == 0) then
     ncall = 1

     do ngr = 1,ngrids
        allocate (pcp_tab(ngr)%pcpfillc(mmzp(ngr),maxkfall,nembfall,nhcat))
        allocate (pcp_tab(ngr)%pcpfillr(mmzp(ngr),maxkfall,nembfall,nhcat))
        allocate (pcp_tab(ngr)%sfcpcp(maxkfall,nembfall,nhcat))
     enddo

     call micinit()
     call make_autotab()
     call haznuc()
     call tabmelt()
     call tabhab()

     do lhcat = 1,nhcat
        ch3(lhcat) = pwvt(lhcat) * pwmasi(lhcat)
        cdp1(lhcat) = pwmasi(lhcat) * (1.5 + .5 * pwvt(lhcat))
        pwvtmasi(lhcat) = pwvt(lhcat) * pwmasi(lhcat)
     enddo
  endif

  if (ncall2g(ngrid) .ne. 5) then
     ncall2g(ngrid) = 5

     call mksedim_tab(mzp,mxp,myp,ngrid,nembfall,maxkfall,zm,dzt  &
          ,pcp_tab(ngrid)%pcpfillc(1,1,1,1),pcp_tab(ngrid)%pcpfillr(1,1,1,1)  &
	  ,pcp_tab(ngrid)%sfcpcp(1,1,1))

     do lhcat = 1,nhcat
        ch2(lhcat,ngrid) = float(nembfall-1) &
             / log10(dispemb1(lhcat,ngrid) / dispemb0(lhcat,ngrid))
     enddo

     call homfrzcl(dtlt,ngrid)
  endif


  call each_call(mzp,dtlt)
  dtlti = 1. / dtlt
  ngr = ngrid


  ! Calling the routine to claculate indexes in the current grid
  call indexing(ngr, ia, iz, ja, jz)


  ! Optimized loop
  step     = min(step_limit, (mxp*myp))

  indextile = 1

  do ijm=1,ij_total,step


     ! Defining ij loop limit (ij_final)
     ij_final = ij_last(indextile,ngr)

     ! Aloca tile
     call alloc_micro_tile(mzp, mxp, myp)

     ! Copiando dados anteriores para os arrays de otimizacao e inicializando
     call init_copy_mic_block(ngr, mzp, mxp, ia, iz, myp, ja, jz, ijm,    &
          micro_g(ngr), basic_g(ngr), grid_g(ngr)%rtgt, grid_g(ngr)%flpw)


     ! ALF - Calling a optimized subroutine for vector machines
     call range_check_opt(mzp)


     call mcphys_opt(mzp,ijm,step,ngrid,maxnzp,          &
          nembfall,maxkfall,mynum,dtlt,dtlti,time,zm,zt, & 
          radiate_g(ngr),pcp_tab(ngr)%pcpfillc,          &
          pcp_tab(ngr)%pcpfillr,pcp_tab(ngr)%sfcpcp,     &
          grid_g(ngr)%glat,grid_g(ngr)%topt,if_adap)


     call copyback_opt()


     call final_copy_mic_block(ngr, mzp, micro_g(ngr), &
          basic_g(ngr), grid_g(ngr)%rtgt, grid_g(ngr)%flpw)


     ! dealoca tile
     call dealloc_micro_tile()

     indextile = indextile + 1


  enddo

  return
end subroutine micro_opt


