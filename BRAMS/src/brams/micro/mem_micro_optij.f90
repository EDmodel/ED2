! BRAMS 3.3 - CATT-BRAMS
! ALF

Module mem_micro_opt

  ! Variables to be dimensioned by (nxp*nyp,7) - Original dimension(7)

  real,    allocatable :: accpx_mod(:,:)
  real,    allocatable :: pcprx_mod(:,:)

  ! Integer Variables to be dimensioned by (nxp*nyp,10)
  ! Original dimension(10)

  integer, allocatable :: k1_mod(:,:)
  integer, allocatable :: k2_mod(:,:)
  integer, allocatable :: k3_mod(:,:)

  ! Variables to be dimensioned by (nxp*nyp,nzp)
  ! Original dimension(nzpmax)

  real,    allocatable :: cccnx_mod(:,:)
  real,    allocatable :: cifnx_mod(:,:)
  real,    allocatable :: scrmic1_mod(:,:)
  real,    allocatable :: scrmic2_mod(:,:)
  real,    allocatable :: scrmic3_mod(:,:)

  ! Variables to be dimensioned by (nxp*nyp,nhcat)
  ! Original dimension(nhcat)

  real,    allocatable :: ch1_mod(:,:)
  real,    allocatable :: cfvt_mod(:,:)

  ! Variables to be dimensioned by (nxp*nyp,nzp,ncat)
  ! Original dimension(nzpmax,ncat)

  real,    allocatable :: rx_mod(:,:,:)
  real,    allocatable :: cx_mod(:,:,:)
  real,    allocatable :: qr_mod(:,:,:)
  real,    allocatable :: qx_mod(:,:,:)
  real,    allocatable :: vap_mod(:,:,:)
  real,    allocatable :: tx_mod(:,:,:)
  real,    allocatable :: emb_mod(:,:,:)
  real,    allocatable :: vterm_mod(:,:,:)
  real,    allocatable :: sd_mod(:,:,:)
  real,    allocatable :: se_mod(:,:,:)
  real,    allocatable :: sf_mod(:,:,:)
  real,    allocatable :: sg_mod(:,:,:)
  real,    allocatable :: sm_mod(:,:,:)
  real,    allocatable :: ss_mod(:,:,:)
  real,    allocatable :: su_mod(:,:,:)
  real,    allocatable :: sw_mod(:,:,:)
  real,    allocatable :: sy_mod(:,:,:)
  real,    allocatable :: sz_mod(:,:,:)
  real,    allocatable :: wct1_mod(:,:,:)
  real,    allocatable :: wct2_mod(:,:,:)

  integer, allocatable :: ict1_mod(:,:,:)
  integer, allocatable :: ict2_mod(:,:,:)

  ! Variables to be dimensioned by (nxp*nyp,nzp,10)
  ! Original dimension(nzpmax,10)

  real, allocatable :: eff_mod(:,:,:)

  ! Variables to be dimensioned by (nxp*nyp,nzp,ncat,ncat)
  ! Original dimension(nzpmax,ncat,ncat)

  real,    allocatable :: rxfer_mod(:,:,:,:)
  real,    allocatable :: qrxfer_mod(:,:,:,:)
  real,    allocatable :: enxfer_mod(:,:,:,:)

  ! Variable to be dimensioned by (nxp*nyp,nzp) - Original (nzpmax)

  real,    allocatable :: pitot_mod(:,:)
  real,    allocatable :: press_mod(:,:)
  real,    allocatable :: tair_mod(:,:)
  real,    allocatable :: til_mod(:,:)
  real,    allocatable :: rliq_mod(:,:)
  real,    allocatable :: rice_mod(:,:)
  real,    allocatable :: qhydm_mod(:,:)
  real,    allocatable :: rvstr_mod(:,:)
  real,    allocatable :: tairstrc_mod(:,:)
  real,    allocatable :: rvlsair_mod(:,:)
  real,    allocatable :: rvisair_mod(:,:)
  real,    allocatable :: dn0i_mod(:,:)
  real,    allocatable :: tairc_mod(:,:)
  real,    allocatable :: thrmcon_mod(:,:)
  real,    allocatable :: dynvisc_mod(:,:)
  real,    allocatable :: vapdif_mod(:,:)
  real,    allocatable :: rdynvsci_mod(:,:)
  real,    allocatable :: denfac_mod(:,:)
  real,    allocatable :: colfacr_mod(:,:)
  real,    allocatable :: colfacr2_mod(:,:)
  real,    allocatable :: colfacc_mod(:,:)
  real,    allocatable :: colfacc2_mod(:,:)
  real,    allocatable :: sumuy_mod(:,:)
  real,    allocatable :: sumuz_mod(:,:)
  real,    allocatable :: sumvr_mod(:,:)
  real,    allocatable :: rvs0_mod(:,:)

  ! Variable to be dimensioned by (nxp*nyp,nzp,ncat)
  ! Original (nzpmax,ncat)

  real,    allocatable :: sh_mod(:,:,:)
  integer, allocatable :: jhcat_mod(:,:,:)

  ! Variable to be dimensioned by (nxp*nyp,nzp,9) - Original (nzpmax,9)

  real,    allocatable :: sa_mod(:,:,:)

  ! Variable to be dimensioned by (nxp*nyp,nzp,2) - Original (nzpmax,2)

  real,    allocatable :: tref_mod(:,:,:)
  real,    allocatable :: rvsref_mod(:,:,:)
  real,    allocatable :: rvsrefp_mod(:,:,:)



  ! tile arrays

  ! Variables to be dimensioned by (nxp*nyp,nzp)(nzp,nxp,nyp)

  real, allocatable :: rcp_tile(:,:)
  real, allocatable :: rrp_tile(:,:)
  real, allocatable :: rpp_tile(:,:)
  real, allocatable :: rsp_tile(:,:)
  real, allocatable :: rap_tile(:,:)
  real, allocatable :: rgp_tile(:,:)
  real, allocatable :: rhp_tile(:,:)
  real, allocatable :: ccp_tile(:,:)
  real, allocatable :: crp_tile(:,:)
  real, allocatable :: cpp_tile(:,:)
  real, allocatable :: csp_tile(:,:)
  real, allocatable :: cap_tile(:,:)
  real, allocatable :: cgp_tile(:,:)
  real, allocatable :: chp_tile(:,:)
  real, allocatable :: cccnp_tile(:,:)
  real, allocatable :: cifnp_tile(:,:)
  real, allocatable :: q2_tile(:,:)
  real, allocatable :: q6_tile(:,:)
  real, allocatable :: q7_tile(:,:)


  ! Variables to be dimensioned by (nxp*nyp,nzp) - Original (n1,n2,n3)
  ! Copy from mem_basic
  
  real, allocatable :: pi0_opt_tile(:,:)
  real, allocatable :: pp_opt_tile(:,:)
  real, allocatable :: theta_opt_tile(:,:)
  real, allocatable :: thp_opt_tile(:,:)
  real, allocatable :: rv_opt_tile(:,:)
  real, allocatable :: rtp_opt_tile(:,:)
  real, allocatable :: wp_opt_tile(:,:)
  real, allocatable :: dn0_opt_tile(:,:)


  ! Variables to be dimensioned by (nxp*nyp) - Original (n2,n3)
  
  real, allocatable :: rtgt_opt_tile(:)
  
  ! Variables to be dimensioned by (nnxp*nyp) - Original (nnxp*nyp)
  
  real, allocatable :: accpr_tile(:)
  real, allocatable :: accpp_tile(:)
  real, allocatable :: accps_tile(:)
  real, allocatable :: accpa_tile(:)
  real, allocatable :: accpg_tile(:)
  real, allocatable :: accph_tile(:)
  real, allocatable :: pcprr_tile(:)
  real, allocatable :: pcprp_tile(:)
  real, allocatable :: pcprs_tile(:)
  real, allocatable :: pcpra_tile(:)
  real, allocatable :: pcprg_tile(:)
  real, allocatable :: pcprh_tile(:)
  real, allocatable :: pcpg_opt_tile(:)
  real, allocatable :: qpcpg_opt_tile(:)
  real, allocatable :: dpcpg_opt_tile(:)



  ! Variable to be dimensioned by *(nxp*nyp) - Original (nxp,nyp)
  integer, allocatable :: lpw_opt(:)

  integer, parameter :: step_limit=256 ! Limit to Vector Pipelines in SX-6
  integer :: i_last=1, j_last=1
  integer :: ij_final, ij_total, ij_dimension, dimtile
  integer, allocatable :: indicei(:,:,:), indicej(:,:,:) ! (ij,tile,grid)
  integer, allocatable :: ijcall_ngrid(:) !(ngrid)
  integer, allocatable :: ij_last(:,:) ! (tile,grid)
  integer              :: indextile !, maxtiles
  integer, allocatable :: maxtiles(:) !, ij_final_lasttile(:)
  integer :: ncall7 = 0 ! To use with subroutine "effxy_opt"
  real    :: old_dtlt = 0 ! New variable introduced in "auto_accret_opt"

contains

  ! ***************************************************************************

  subroutine alloc_micro_opt(n1,n2,n3)

    use micphys, only: ncat, nhcat, level,  &
         irain, ipris, isnow, iaggr, igraup, ihail, icloud  ! INTENT(IN)
    use mem_grid, only: ngrids             ! INTENT(IN)

    implicit none

    ! Arguments:
    integer, dimension(:), intent(in) :: n1,n2,n3

    ! Local Variables:
    integer :: ng, maxn1, maxn2, maxn3 !, maxtiles

    maxn1 = maxval(n1(1:ngrids))
    maxn2 = maxval(n2(1:ngrids))
    maxn3 = maxval(n3(1:ngrids))

    ij_dimension = min(step_limit, maxn2*maxn3)

    ! Allocation for the adjacent matrix
    ! (indexes for maping i,j in ij betwin tiles and grids)
    allocate (maxtiles(ngrids))
    allocate (ijcall_ngrid(ngrids))

    ijcall_ngrid(:) = 0 ! Initializing

    do ng = 1, ngrids
       maxtiles(ng) = ceiling(real( (n2(ng)-2)*(n3(ng)-2) ) / &
            real(min(step_limit, (n2(ng)-2)*(n3(ng)-2) ) ) )
    enddo

    allocate (ij_last(maxval(maxtiles), ngrids))
    allocate (indicei(ij_dimension, maxval(maxtiles), ngrids))
    allocate (indicej(ij_dimension, maxval(maxtiles), ngrids))

    ! Allocating module global variables
    allocate (accpx_mod(ij_dimension,7))
    allocate (pcprx_mod(ij_dimension,7))
    allocate (k1_mod(ij_dimension, 10))
    allocate (k2_mod(ij_dimension, 10))
    allocate (k3_mod(ij_dimension, 10))
    allocate (cccnx_mod(ij_dimension,maxn1))
    allocate (cifnx_mod(ij_dimension,maxn1))
    allocate (scrmic1_mod(ij_dimension,maxn1))
    allocate (scrmic2_mod(ij_dimension,maxn1))
    allocate (scrmic3_mod(ij_dimension,maxn1))
    allocate (ch1_mod(ij_dimension,nhcat))
    allocate (cfvt_mod(ij_dimension,nhcat))
    allocate (rx_mod(ij_dimension,maxn1,ncat))
    allocate (cx_mod(ij_dimension,maxn1,ncat))
    allocate (qr_mod(ij_dimension,maxn1,ncat))
    allocate (qx_mod(ij_dimension,maxn1,ncat))
    allocate (vap_mod(ij_dimension,maxn1,ncat))
    allocate (tx_mod(ij_dimension,maxn1,ncat))
    allocate (emb_mod(ij_dimension,maxn1,ncat))
    allocate (vterm_mod(ij_dimension,maxn1,ncat))
    allocate (sd_mod(ij_dimension,maxn1,ncat))
    allocate (se_mod(ij_dimension,maxn1,ncat))
    allocate (sf_mod(ij_dimension,maxn1,ncat))
    allocate (sg_mod(ij_dimension,maxn1,ncat))
    allocate (sm_mod(ij_dimension,maxn1,ncat))
    allocate (ss_mod(ij_dimension,maxn1,ncat))
    allocate (su_mod(ij_dimension,maxn1,ncat))
    allocate (sw_mod(ij_dimension,maxn1,ncat))
    allocate (sy_mod(ij_dimension,maxn1,ncat))
    allocate (sz_mod(ij_dimension,maxn1,ncat))
    allocate (wct1_mod(ij_dimension,maxn1,ncat))
    allocate (wct2_mod(ij_dimension,maxn1,ncat))
    allocate (ict1_mod(ij_dimension,maxn1,ncat))
    allocate (ict2_mod(ij_dimension,maxn1,ncat))
    allocate (eff_mod(ij_dimension,maxn1,10))
    allocate (rxfer_mod(ij_dimension,maxn1,ncat,ncat))
    allocate (qrxfer_mod(ij_dimension,maxn1,ncat,ncat))
    allocate (enxfer_mod(ij_dimension,maxn1,ncat,ncat))
    allocate (pitot_mod(ij_dimension,maxn1))
    allocate (press_mod(ij_dimension,maxn1))
    allocate (tair_mod(ij_dimension,maxn1))
    allocate (til_mod(ij_dimension,maxn1))
    allocate (rliq_mod(ij_dimension,maxn1))
    allocate (rice_mod(ij_dimension,maxn1))
    allocate (qhydm_mod(ij_dimension,maxn1))
    allocate (rvstr_mod(ij_dimension,maxn1))
    allocate (tairstrc_mod(ij_dimension,maxn1))
    allocate (rvlsair_mod(ij_dimension,maxn1))
    allocate (rvisair_mod(ij_dimension,maxn1))
    allocate (dn0i_mod(ij_dimension,maxn1))
    allocate (tairc_mod(ij_dimension,maxn1))
    allocate (thrmcon_mod(ij_dimension,maxn1))
    allocate (dynvisc_mod(ij_dimension,maxn1))
    allocate (vapdif_mod(ij_dimension,maxn1))
    allocate (rdynvsci_mod(ij_dimension,maxn1))
    allocate (denfac_mod(ij_dimension,maxn1))
    allocate (colfacr_mod(ij_dimension,maxn1))
    allocate (colfacr2_mod(ij_dimension,maxn1))
    allocate (colfacc_mod(ij_dimension,maxn1))
    allocate (colfacc2_mod(ij_dimension,maxn1))
    allocate (sumuy_mod(ij_dimension,maxn1))
    allocate (sumuz_mod(ij_dimension,maxn1))
    allocate (sumvr_mod(ij_dimension,maxn1))
    allocate (rvs0_mod(ij_dimension,maxn1))
    allocate (sh_mod(ij_dimension,maxn1,ncat))
    allocate (jhcat_mod(ij_dimension,maxn1,ncat))
    allocate (sa_mod(ij_dimension,maxn1,9))
    allocate (tref_mod(ij_dimension,maxn1,2))
    allocate (rvsref_mod(ij_dimension,maxn1,2))
    allocate (rvsrefp_mod(ij_dimension,maxn1,2))

    ! Allocation for lpw_opt
    allocate (lpw_opt(ij_dimension))


  end subroutine alloc_micro_opt

  ! ***************************************************************************

  subroutine dealloc_micro_opt()

    use mem_grid, only: ngrids             ! INTENT(IN)

    implicit none

    ! Local Variables:
    integer :: ng

    ! Adjacent matrix - indexes for maping i,j in ij
    deallocate (indicei)
    deallocate (indicej)

    ! LPW
    deallocate (lpw_opt)

    ! Deallocating module global variables
    if (allocated(accpx_mod)) deallocate(accpx_mod)
    if (allocated(pcprx_mod)) deallocate(pcprx_mod)
    if (allocated(k1_mod)) deallocate(k1_mod)
    if (allocated(k2_mod)) deallocate(k2_mod)
    if (allocated(k3_mod)) deallocate(k3_mod)
    if (allocated(cccnx_mod)) deallocate(cccnx_mod)
    if (allocated(cifnx_mod)) deallocate(cifnx_mod)
    if (allocated(scrmic1_mod)) deallocate(scrmic1_mod)
    if (allocated(scrmic2_mod)) deallocate(scrmic2_mod)
    if (allocated(scrmic3_mod)) deallocate(scrmic3_mod)
    if (allocated(ch1_mod)) deallocate(ch1_mod)
    if (allocated(cfvt_mod)) deallocate(cfvt_mod)
    if (allocated(rx_mod)) deallocate(rx_mod)
    if (allocated(cx_mod)) deallocate(cx_mod)
    if (allocated(qr_mod)) deallocate(qr_mod)
    if (allocated(qx_mod)) deallocate(qx_mod)
    if (allocated(vap_mod)) deallocate(vap_mod)
    if (allocated(tx_mod)) deallocate(tx_mod)
    if (allocated(emb_mod)) deallocate(emb_mod)
    if (allocated(vterm_mod)) deallocate(vterm_mod)
    if (allocated(sd_mod)) deallocate(sd_mod)
    if (allocated(se_mod)) deallocate(se_mod)
    if (allocated(sf_mod)) deallocate(sf_mod)
    if (allocated(sg_mod)) deallocate(sg_mod)
    if (allocated(sm_mod)) deallocate(sm_mod)
    if (allocated(ss_mod)) deallocate(ss_mod)
    if (allocated(su_mod)) deallocate(su_mod)
    if (allocated(sw_mod)) deallocate(sw_mod)
    if (allocated(sy_mod)) deallocate(sy_mod)
    if (allocated(sz_mod)) deallocate(sz_mod)
    if (allocated(wct1_mod)) deallocate(wct1_mod)
    if (allocated(wct2_mod)) deallocate(wct2_mod)
    if (allocated(ict1_mod)) deallocate(ict1_mod)
    if (allocated(ict2_mod)) deallocate(ict2_mod)
    if (allocated(eff_mod)) deallocate(eff_mod)
    if (allocated(rxfer_mod)) deallocate(rxfer_mod)
    if (allocated(qrxfer_mod)) deallocate(qrxfer_mod)
    if (allocated(enxfer_mod)) deallocate(enxfer_mod)
    if (allocated(pitot_mod)) deallocate(pitot_mod)
    if (allocated(press_mod)) deallocate(press_mod)
    if (allocated(tair_mod)) deallocate(tair_mod)
    if (allocated(til_mod)) deallocate(til_mod)
    if (allocated(rliq_mod)) deallocate(rliq_mod)
    if (allocated(rice_mod)) deallocate(rice_mod)
    if (allocated(qhydm_mod)) deallocate(qhydm_mod)
    if (allocated(rvstr_mod)) deallocate(rvstr_mod)
    if (allocated(tairstrc_mod)) deallocate(tairstrc_mod)
    if (allocated(rvlsair_mod)) deallocate(rvlsair_mod)
    if (allocated(rvisair_mod)) deallocate(rvisair_mod)
    if (allocated(dn0i_mod)) deallocate(dn0i_mod)
    if (allocated(tairc_mod)) deallocate(tairc_mod)
    if (allocated(thrmcon_mod)) deallocate(thrmcon_mod)
    if (allocated(dynvisc_mod)) deallocate(dynvisc_mod)
    if (allocated(vapdif_mod)) deallocate(vapdif_mod)
    if (allocated(rdynvsci_mod)) deallocate(rdynvsci_mod)
    if (allocated(denfac_mod)) deallocate(denfac_mod)
    if (allocated(colfacr_mod)) deallocate(colfacr_mod)
    if (allocated(colfacr2_mod)) deallocate(colfacr2_mod)
    if (allocated(colfacc_mod)) deallocate(colfacc_mod)
    if (allocated(colfacc2_mod)) deallocate(colfacc2_mod)
    if (allocated(sumuy_mod)) deallocate(sumuy_mod)
    if (allocated(sumuz_mod)) deallocate(sumuz_mod)
    if (allocated(sumvr_mod)) deallocate(sumvr_mod)
    if (allocated(rvs0_mod)) deallocate(rvs0_mod)
    if (allocated(sh_mod)) deallocate(sh_mod)
    if (allocated(jhcat_mod)) deallocate(jhcat_mod)
    if (allocated(sa_mod)) deallocate(sa_mod)
    if (allocated(tref_mod)) deallocate(tref_mod)
    if (allocated(rvsref_mod)) deallocate(rvsref_mod)
    if (allocated(rvsrefp_mod)) deallocate(rvsrefp_mod)
    
  end subroutine dealloc_micro_opt

  !****************************************************************************

  subroutine alloc_micro_tile(n1, n2, n3)
    use micphys, only: &
         level,        &
         irain,        &
         ipris,        &
         isnow,        &
         iaggr,        &
         igraup,       &
         ihail,        &
         icloud 

    implicit none
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: n3

    ij_dimension = min(step_limit, (n2*n3))
    
    allocate (pi0_opt_tile(ij_dimension,n1))
    allocate (pp_opt_tile(ij_dimension,n1))
    allocate (theta_opt_tile(ij_dimension,n1))
    allocate (thp_opt_tile(ij_dimension,n1))
    allocate (rv_opt_tile(ij_dimension,n1))
    allocate (rtp_opt_tile(ij_dimension,n1))
    allocate (wp_opt_tile(ij_dimension,n1))
    allocate (dn0_opt_tile(ij_dimension,n1))
    allocate (rtgt_opt_tile(ij_dimension))
    
    ! Allocate arrays based on options (if necessary)
    if (level >= 2 ) then
       allocate (rcp_tile(ij_dimension,n1))
    endif
    if (level >= 3) then
       if(irain >= 1)  then
          allocate (rrp_tile(ij_dimension,n1))
          allocate (accpr_tile(ij_dimension))
          allocate (pcprr_tile(ij_dimension))
          allocate (q2_tile(ij_dimension,n1))
       endif
       if(ipris >= 1)  then
          allocate (rpp_tile(ij_dimension,n1))
          allocate (accpp_tile(ij_dimension))
          allocate (pcprp_tile(ij_dimension))
       endif
       if(isnow >= 1)  then
          allocate (rsp_tile(ij_dimension,n1))
          allocate (accps_tile(ij_dimension))
          allocate (pcprs_tile(ij_dimension))
       endif
       if(iaggr >= 1)  then
          allocate (rap_tile(ij_dimension,n1))
          allocate (accpa_tile(ij_dimension))
          allocate (pcpra_tile(ij_dimension))
       endif
       if(igraup >= 1) then
          allocate (rgp_tile(ij_dimension,n1))
          allocate (accpg_tile(ij_dimension))
          allocate (pcprg_tile(ij_dimension))
          allocate (q6_tile(ij_dimension,n1))
       endif
       if(ihail >= 1)  then
          allocate (rhp_tile(ij_dimension,n1))
          allocate (accph_tile(ij_dimension))
          allocate (pcprh_tile(ij_dimension))
          allocate (q7_tile(ij_dimension,n1))
       endif
       if(icloud == 5) allocate (ccp_tile(ij_dimension,n1))
       if(irain == 5)  allocate (crp_tile(ij_dimension,n1))
       if(ipris == 5)  allocate (cpp_tile(ij_dimension,n1))
       if(isnow == 5)  allocate (csp_tile(ij_dimension,n1))
       if(iaggr == 5)  allocate (cap_tile(ij_dimension,n1))
       if(igraup == 5) allocate (cgp_tile(ij_dimension,n1))
       if(ihail == 5)  allocate (chp_tile(ij_dimension,n1))
       allocate (cccnp_tile(ij_dimension,n1))
       allocate (cifnp_tile(ij_dimension,n1))
       allocate (pcpg_opt_tile(ij_dimension))
       allocate (qpcpg_opt_tile(ij_dimension))
       allocate (dpcpg_opt_tile(ij_dimension))
       
    endif

  end subroutine alloc_micro_tile

  !****************************************************************************

  subroutine dealloc_micro_tile
    implicit none
    if (allocated(pi0_opt_tile))     deallocate (pi0_opt_tile)
    if (allocated(pp_opt_tile))      deallocate (pp_opt_tile)
    if (allocated(theta_opt_tile))   deallocate (theta_opt_tile)
    if (allocated(thp_opt_tile))     deallocate (thp_opt_tile)
    if (allocated(rv_opt_tile))      deallocate (rv_opt_tile)
    if (allocated(rtp_opt_tile))     deallocate (rtp_opt_tile)
    if (allocated(wp_opt_tile))      deallocate (wp_opt_tile)
    if (allocated(dn0_opt_tile))     deallocate (dn0_opt_tile)
    if (allocated(rtgt_opt_tile))    deallocate (rtgt_opt_tile)
    if (allocated(rcp_tile))         deallocate (rcp_tile)
    if (allocated(rrp_tile))         deallocate (rrp_tile)
    if (allocated(rpp_tile))         deallocate (rpp_tile)
    if (allocated(rsp_tile))         deallocate (rsp_tile)
    if (allocated(rap_tile))         deallocate (rap_tile)
    if (allocated(rgp_tile))         deallocate (rgp_tile)
    if (allocated(rhp_tile))         deallocate (rhp_tile)
    if (allocated(ccp_tile))         deallocate (ccp_tile)
    if (allocated(crp_tile))         deallocate (crp_tile)
    if (allocated(cpp_tile))         deallocate (cpp_tile)
    if (allocated(csp_tile))         deallocate (csp_tile)
    if (allocated(cap_tile))         deallocate (cap_tile)
    if (allocated(cgp_tile))         deallocate (cgp_tile)
    if (allocated(chp_tile))         deallocate (chp_tile)
    if (allocated(cccnp_tile))       deallocate (cccnp_tile)
    if (allocated(cifnp_tile))       deallocate (cifnp_tile)
    if (allocated(q2_tile))          deallocate (q2_tile)
    if (allocated(q6_tile))          deallocate (q6_tile)
    if (allocated(q7_tile))          deallocate (q7_tile)
    if (allocated(accpr_tile))       deallocate (accpr_tile)
    if (allocated(accpp_tile))       deallocate (accpp_tile)
    if (allocated(accps_tile))       deallocate (accps_tile)
    if (allocated(accpa_tile))       deallocate (accpa_tile)
    if (allocated(accpg_tile))       deallocate (accpg_tile)
    if (allocated(accph_tile))       deallocate (accph_tile)
    if (allocated(pcprr_tile))       deallocate (pcprr_tile)
    if (allocated(pcprp_tile))       deallocate (pcprp_tile)
    if (allocated(pcprs_tile))       deallocate (pcprs_tile)
    if (allocated(pcpra_tile))       deallocate (pcpra_tile)
    if (allocated(pcprg_tile))       deallocate (pcprg_tile)
    if (allocated(pcprh_tile))       deallocate (pcprh_tile)
    if (allocated(pcpg_opt_tile))    deallocate (pcpg_opt_tile)
    if (allocated(qpcpg_opt_tile))   deallocate (qpcpg_opt_tile)
    if (allocated(dpcpg_opt_tile))   deallocate (dpcpg_opt_tile)
  end subroutine dealloc_micro_tile

  ! ***************************************************************************


  subroutine range_check_opt(m1)
    ! Using: k1_mod, k2_mod, k3_mod 
    ! INTENT(OUT)
    ! Using lpw_opt substituing lpw

    use micphys, only: &
         ncat,         & ! INTENT(IN)
         jnmb            ! INTENT(IN)

    implicit none

    ! Arguments
    integer, intent(in)               :: m1

    ! Local Variables
    integer                           :: k, lcat, l, jcat, ij


    ! zero out microphysics scratch arrays for the present i,j column

    do lcat = 1,ncat
       do k = 2,m1-1
          do ij=1,ij_final
             rx_mod(ij,k,lcat)  = 0.
             cx_mod(ij,k,lcat)  = 0.
             qr_mod(ij,k,lcat)  = 0.
             qx_mod(ij,k,lcat)  = 0.
             vap_mod(ij,k,lcat) = 0.
             tx_mod(ij,k,lcat)  = 0.
          enddo
       enddo

       if (jnmb(lcat) >= 3) then
          do k = 2,m1-1
             do ij=1,ij_final
                emb_mod(ij,k,lcat) = 0.
             enddo
          enddo
       endif

       do jcat = 1,ncat
          do k = 2,m1-1
             do ij=1,ij_final
                rxfer_mod(ij,k,lcat,jcat)  = 0.
                qrxfer_mod(ij,k,lcat,jcat) = 0.
                enxfer_mod(ij,k,lcat,jcat) = 0.
             enddo
          enddo
       enddo
    enddo

    do l = 1,7
       do ij=1,ij_final
          k1_mod(ij,l) = lpw_opt(ij)
          k2_mod(ij,l) = 1
       enddo
    enddo

    ! fill scratch arrays for cloud water

    if (jnmb(1) >= 1) then
       do k = max(minval(lpw_opt),1),m1-1
          do ij=1,ij_final
             if (k >= lpw_opt(ij)) then
                if (rcp_tile(ij,k) >= 1.e-9) then
                   k2_mod(ij,1) = k
                   rx_mod(ij,k,1) = rcp_tile(ij,k)
                   if (jnmb(1) >= 5)  &
                        cx_mod(ij,k,1)  = ccp_tile(ij,k)
                   if (jnmb(1) == 7)  &
                        cccnx_mod(ij,k) = cccnp_tile(ij,k)
                else
                   if (k2_mod(ij,1) == 1) &
                        k1_mod(ij,1) = k+1
                endif
             endif
          enddo
       enddo
    endif

    ! fill scratch arrays for rain

    if (jnmb(2) >= 1) then
       do k = minval(lpw_opt),m1-1
          do ij=1,ij_final
             if (k >= lpw_opt(ij)) then
                if (rrp_tile(ij,k) >= 1.e-9) then
                   k2_mod(ij,2) = k
                   rx_mod(ij,k,2) = rrp_tile(ij,k)
                   qx_mod(ij,k,2) = q2_tile(ij,k)
                   qr_mod(ij,k,2) =                   &
                        qx_mod(ij,k,2) *              &
                        rx_mod(ij,k,2)
                   if (jnmb(2) >= 5) cx_mod(ij,k,2) = &
                        crp_tile(ij,k)
                else
                   if (k2_mod(ij,2) == 1) &
                        k1_mod(ij,2) = k+1
                endif
             endif
          enddo
       enddo
    endif

    ! fill scratch arrays for pristine ice

    if (jnmb(3) >= 1) then
       do k = minval(lpw_opt),m1-1
          do ij=1,ij_final
             if (k >= lpw_opt(ij)) then
                if (rpp_tile(ij,k) >= 1.e-9) then
                   k2_mod(ij,3) = k
                   rx_mod(ij,k,3) = rpp_tile(ij,k)
                   cx_mod(ij,k,3) = cpp_tile(ij,k)
                   if (jnmb(3) == 7)  &
                        cifnx_mod(ij,k) = cifnp_tile(ij,k)
                else
                   if (k2_mod(ij,3) == 1) &
                        k1_mod(ij,3) = k+1
                endif
             endif
          enddo
       enddo
    endif

    ! fill scratch arrays for snow

    if (jnmb(4) >= 1) then
       do k = minval(lpw_opt),m1-1
          do ij=1,ij_final
             if (k >= lpw_opt(ij)) then
                if (rsp_tile(ij,k) >= 1.e-9) then
                   k2_mod(ij,4) = k
                   rx_mod(ij,k,4) = rsp_tile(ij,k)
                   if (jnmb(4) >= 5) cx_mod(ij,k,4) = &
                        csp_tile(ij,k)
                else
                   if (k2_mod(ij,4) == 1) &
                        k1_mod(ij,4) = k+1
                endif
             endif
          enddo
       enddo
    endif

    ! fill scratch arrays for aggregates

    if (jnmb(5) >= 1) then
       do k = minval(lpw_opt),m1-1
          do ij=1,ij_final
             if (k >= lpw_opt(ij)) then
                if (rap_tile(ij,k) >= 1.e-9) then
                   k2_mod(ij,5) = k
                   rx_mod(ij,k,5) = rap_tile(ij,k)
                   if (jnmb(5) >= 5) cx_mod(ij,k,5) = &
                        cap_tile(ij,k)
                else
                   if (k2_mod(ij,5) == 1) &
                        k1_mod(ij,5) = k+1
                endif
             endif
          enddo
       enddo
    endif

    ! fill scratch arrays for graupel

    if (jnmb(6) >= 1) then
       do k = minval(lpw_opt),m1-1
          do ij=1,ij_final
             if (k >= lpw_opt(ij)) then
                if (rgp_tile(ij,k) >= 1.e-9) then
                   k2_mod(ij,6) = k
                   rx_mod(ij,k,6) = rgp_tile(ij,k)
                   qx_mod(ij,k,6) = q6_tile(ij,k)
                   qr_mod(ij,k,6) =                   &
                        qx_mod(ij,k,6) *              &
                        rx_mod(ij,k,6)
                   if (jnmb(6) >= 5) cx_mod(ij,k,6) = &
                        cgp_tile(ij,k)
                else
                   if (k2_mod(ij,6) == 1) &
                        k1_mod(ij,6) = k+1
                endif
             endif
          enddo
       enddo
    endif

    ! fill scratch arrays for hail

    if (jnmb(7) >= 1) then
       do k = minval(lpw_opt),m1-1
          do ij=1,ij_final
             if (k >= lpw_opt(ij)) then
                if (rhp_tile(ij,k) >= 1.e-9) then
                   k2_mod(ij,7) = k
                   rx_mod(ij,k,7) = rhp_tile(ij,k)
                   qx_mod(ij,k,7) = q7_tile(ij,k)
                   qr_mod(ij,k,7) =                   &
                        qx_mod(ij,k,7) *              &
                        rx_mod(ij,k,7)
                   if (jnmb(7) >= 5) cx_mod(ij,k,7) = &
                        chp_tile(ij,k)
                else
                   if (k2_mod(ij,7) == 1) &
                        k1_mod(ij,7) = k+1
                endif
             endif
          enddo
       enddo
    endif


    do ij=1,ij_final

       k3_mod(ij,1)  = k2_mod(ij,1)
       k3_mod(ij,3)  = k2_mod(ij,3)
       k1_mod(ij,8)  = min(k1_mod(ij,1), &
            k1_mod(ij,2))
       k2_mod(ij,8)  = max(k2_mod(ij,1), &
            k2_mod(ij,2))
       k1_mod(ij,9)  = min(k1_mod(ij,3), &
            k1_mod(ij,4), k1_mod(ij,5),  &
            k1_mod(ij,6), k1_mod(ij,7))
       k2_mod(ij,9)  = max(k2_mod(ij,3), &
            k2_mod(ij,4), k2_mod(ij,5),  &
            k2_mod(ij,6), k2_mod(ij,7))
       k1_mod(ij,10) = min(k1_mod(ij,8), &
            k1_mod(ij,9))
       k2_mod(ij,10) = max(k2_mod(ij,8), &
            k2_mod(ij,9))

    enddo

    return

  end subroutine range_check_opt

  !****************************************************************************

  subroutine copyback_opt()

    use mem_micro, only: &
         micro_vars        ! INTENT(IN) ! Only a type structure

    use micphys, only:   &
         jnmb              ! INTENT(IN)

    implicit none

    ! Arguments:

    ! Local Variables:
    integer :: ij, lpw_min, k2_max, k3_max


    lpw_min = minval(lpw_opt)
    k2_max  = maxval(k2_mod(:,10))
    k3_max  = maxval(k3_mod(:,1))

    if (jnmb(1) >= 1) then
       call copmic_k3(lpw_min, k3_max, 1, rcp_tile, rx_mod, 1)

       if (jnmb(1) >= 5)  &
            call copmic_k3(lpw_min, k3_max, 1, ccp_tile, cx_mod,1)
    endif

    if (jnmb(2) >= 1) then
       call copmic_k2(lpw_min, k2_max, 10, q2_tile, qx_mod, 2)
       call copmic_k2(lpw_min, k2_max, 10, rrp_tile, rx_mod, 2)

       do ij=1,ij_final
          accpr_tile(ij) = accpr_tile(ij) + accpx_mod(ij,2)
          pcprr_tile(ij) = pcprx_mod(ij,2)
       enddo

       if (jnmb(2) >= 5)  &
            call copmic_k3(lpw_min, k3_max, 1, crp_tile, cx_mod,2)
    endif

    if (jnmb(3) >= 1) then
       k3_max  = maxval(k3_mod(:,3))

       call copmic_k3(lpw_min, k3_max, 3, rpp_tile, rx_mod, 3)

       do ij=1, ij_final
          accpp_tile(ij) = accpp_tile(ij) + accpx_mod(ij,3)
          pcprp_tile(ij) = pcprx_mod(ij,3)
       enddo

       if (jnmb(3) >= 5)  &
            call copmic_k3(lpw_min, k3_max, 3, cpp_tile, cx_mod,3)

    endif

    if (jnmb(4) >= 1) then
       call copmic_k2(lpw_min, k2_max, 10, rsp_tile, rx_mod, 4)

       do ij=1,ij_final
          accps_tile(ij) = accps_tile(ij) + accpx_mod(ij,4)
          pcprs_tile(ij) = pcprx_mod(ij,4)
       enddo

       if (jnmb(4) >= 5)  &
            call copmic_k2(lpw_min, k2_max,10, csp_tile, cx_mod,4)
    endif

    if (jnmb(5) >= 1) then
       call copmic_k2(lpw_min, k2_max, 10, rap_tile, rx_mod, 5)

       do ij=1,ij_final
          accpa_tile(ij) = accpa_tile(ij) + accpx_mod(ij,5)
          pcpra_tile(ij) = pcprx_mod(ij,5)
       enddo

       if (jnmb(5) >= 5)  &
            call copmic_k2(lpw_min, k2_max, 10, cap_tile,cx_mod,5)
    endif

    if (jnmb(6) >= 1) then
       call copmic_k2(lpw_min, k2_max, 10, rgp_tile, rx_mod, 6)

       call copmic_k2(lpw_min, k2_max, 10, q6_tile, qx_mod, 6)

       do ij=1,ij_final
          accpg_tile(ij) = accpg_tile(ij) + accpx_mod(ij,6)
          pcprg_tile(ij) = pcprx_mod(ij,6)
       enddo

       if (jnmb(6) >= 5)   &
            call copmic_k2(lpw_min, k2_max, 10, cgp_tile,cx_mod,6)
    endif

    if (jnmb(7) >= 1) then
       call copmic_k2(lpw_min, k2_max, 10, rhp_tile, rx_mod, 7)

       call copmic_k2(lpw_min, k2_max, 10, q7_tile, qx_mod, 7)

       do ij=1, ij_final
          accph_tile(ij) = accph_tile(ij) + accpx_mod(ij,7)
          pcprh_tile(ij) = pcprx_mod(ij,7)
       enddo

       if (jnmb(7) >= 5)   &
            call copmic_k2(lpw_min, k2_max,10, chp_tile, cx_mod,7)
    endif

    return
  end subroutine copyback_opt

  !****************************************************************************

  subroutine copmic_k2(k_min, k_max, iden, cr3, cr1, dim3)

    implicit none

    ! Arguments
    integer, intent(in)                  :: &
         k_min, k_max, iden, dim3
    real, dimension(:, :), intent(out)   :: cr3
    real, dimension(:, :, :), intent(in) :: cr1

    ! Local Variables
    integer :: k, ij

    do k = k_min, k_max
       do ij=1, ij_final
          if (k>=lpw_opt(ij) .and. k<=k2_mod(ij,iden))  then
             cr3(ij, k) = cr1(ij, k, dim3)
          endif
       enddo
    enddo

    return
  end subroutine copmic_k2

  !****************************************************************************

  subroutine copmic_k3(k_min, k_max, iden, cr3, cr1, dim3)

    implicit none

    ! Arguments
    integer, intent(in)                  :: &
         k_min, k_max, iden, dim3
    real, dimension(:, :), intent(out)   :: cr3
    real, dimension(:, :, :), intent(in) :: cr1

    ! Local Variables
    integer :: k, ij

    do k = k_min, k_max
       do ij=1, ij_final
          if (k>=lpw_opt(ij) .and. k<=k3_mod(ij,iden))  then
             cr3(ij, k) = cr1(ij, k, dim3)
          endif
       enddo
    enddo

    return
  end subroutine copmic_k3

  !***************************************************************************

  subroutine indexing(ngr, ia, iz, ja, jz)

    ! Arguments
    integer :: ngr, ia, iz, ja, jz

    ! Local Varables
    integer :: indextile, step, ijm, ij, i, j


    ! Calculating the max number of elements
    ij_total  = ((iz-ia+1)*(jz-ja+1)) !(mxp*myp)


    if (ijcall_ngrid(ngr) == 0) then !Checking if current grid is indexed

       ijcall_ngrid(ngr) = 1

       ! Calculating the apropriated step
       step      = min(step_limit, ij_total)  !min(step_limit, (mxp*myp))

       ! Calculating max number of tiles
       dimtile   = ceiling(real(ij_total)/real(step))

       ! Finding limits for each tile
       indextile = 1
       limits: do ijm=1,ij_total,step
          ! Defining ij loop limit (ij_final)
          ij_final = min(step,(ij_total-ijm+1)) 
          !                   ,((iz-i_last+1)*(jz-j_last+1)))
          ij_last(indextile,ngr) = ij_final
          indextile = indextile + 1
       enddo limits

       ! Indexing Tiles for current grid
       indextile = 1
       ij = 0

       column: do j=ja, jz

          line: do i=ia, iz

             ij = ij + 1
             if (ij > ij_last(indextile,ngr)) then
                ij = 1
                indextile = indextile + 1
                if (indextile > dimtile) then
                   stop "Erro!!!!! Indexing!!!"
                end if
             end if
             ! 1 <= ij <= ij_final
             indicei(ij,indextile,ngr) = i
             indicej(ij,indextile,ngr) = j
             !ijlast(indextile,ngr) = ij

          enddo line

       enddo column

    endif

  end subroutine indexing

  !***************************************************************************

  subroutine init_copy_mic_block(ngr, m1, m2, ia, iz, m3, ja, jz, ijm, &
       micro, basic, rtgt, flpw)

    use mem_micro, only:  &
         micro_vars            ! INTENT(IN)

    use mem_basic, only:  &
         basic_vars            ! INTENT(IN)

    use micphys, only:  &
         level,icloud,irain,ipris,isnow,iaggr,igraup,ihail,ncat, & !INTENT(IN)
         scrmic1,         & !INTENT(IN)
         scrmic2,         & !INTENT(IN)
         emb,             & !INTENT(IN)
         pitot,           & !INTENT(IN)
         press,           & !INTENT(IN)
         tair,            & !INTENT(IN)
         til,             & !INTENT(IN)
         rliq,            & !INTENT(IN)
         rice,            & !INTENT(IN)
         cx,              & !INTENT(IN)
         rx,              & !INTENT(IN)
         qx,              & !INTENT(IN)
         sm,              & !INTENT(IN)
         vap,             & !INTENT(IN)
         qhydm,           & !INTENT(IN)
         rvstr,           & !INTENT(IN)
         sa,              & !INTENT(IN)
         tairstrc,        & !INTENT(IN)
         tref,            & !INTENT(IN)
         vapdif,          & !INTENT(IN)
         eff,             & !INTENT(IN)
         sh,              & !INTENT(IN)
         dn0i,            & !INTENT(IN)
         jhcat              !INTENT(IN)

    use mem_grid, only: ngrids !INTENT(IN)

    implicit none

    ! Arguments
    type (micro_vars), intent(in)      :: micro 
    type (basic_vars), intent(in)      :: basic
    integer, intent(in)                :: m1,m2,ia,iz,m3,ja,jz,ijm, ngr
    real, intent(in)                   :: flpw(:,:)
    real, intent(in)                   :: rtgt(:,:)

    ! Local Variables
    integer :: i, j, k, ij, cat

    ! *** 2D Arrays:

    do ij=1,ij_final

       lpw_opt(ij) = nint(flpw(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)))
       rtgt_opt_tile(ij) = rtgt(indicei(ij,indextile,ngr),                 &
            indicej(ij,indextile,ngr))

       if (level >= 3) then
          if(irain >= 1)  then
             accpr_tile(ij) = micro%accpr(indicei(ij,indextile,ngr),       &
                  indicej(ij,indextile,ngr))
             pcprr_tile(ij) = micro%pcprr(indicei(ij,indextile,ngr),       &
                  indicej(ij,indextile,ngr))
          endif
          if(ipris >= 1)  then
             accpp_tile(ij) = micro%accpp(indicei(ij,indextile,ngr),       &
                  indicej(ij,indextile,ngr))
             pcprp_tile(ij) = micro%pcprp(indicei(ij,indextile,ngr),       &
                  indicej(ij,indextile,ngr))
          endif
          if(isnow >= 1)  then
             accps_tile(ij) = micro%accps(indicei(ij,indextile,ngr),       &
                  indicej(ij,indextile,ngr))
             pcprs_tile(ij) = micro%pcprs(indicei(ij,indextile,ngr),       &
                  indicej(ij,indextile,ngr))
          endif
          if(iaggr >= 1)  then
             accpa_tile(ij) = micro%accpa(indicei(ij,indextile,ngr),       &
                  indicej(ij,indextile,ngr))
             pcpra_tile(ij) = micro%pcpra(indicei(ij,indextile,ngr),       &
                  indicej(ij,indextile,ngr))
          endif
          if(igraup >= 1) then
             accpg_tile(ij) = micro%accpg(indicei(ij,indextile,ngr),       &
                  indicej(ij,indextile,ngr))
             pcprg_tile(ij) = micro%pcprg(indicei(ij,indextile,ngr),       &
                  indicej(ij,indextile,ngr))
          endif
          if(ihail >= 1)  then
             accph_tile(ij) = micro%accph(indicei(ij,indextile,ngr),       &
                  indicej(ij,indextile,ngr))
             pcprh_tile(ij) = micro%pcprh(indicei(ij,indextile,ngr),       &
                  indicej(ij,indextile,ngr))
          endif
          pcpg_opt_tile(ij)    = micro%pcpg(indicei(ij,indextile,ngr),     &
               indicej(ij,indextile,ngr))
          qpcpg_opt_tile(ij)   = micro%qpcpg(indicei(ij,indextile,ngr),    &
               indicej(ij,indextile,ngr))
          dpcpg_opt_tile(ij)   = micro%dpcpg(indicei(ij,indextile,ngr),    &
               indicej(ij,indextile,ngr))
       endif

    enddo

    ! *** 3D Arrays:

    cccnx_mod(:,:) = 0. ! cldnuc
    cifnx_mod(:,:) = 0. ! icenuc

    do k = 1,m1

       do ij=1,ij_final

          ! Input data used in thrmstr_opt
          pp_opt_tile(ij,k)    = basic%pp(k,indicei(ij,indextile,ngr),     &
               indicej(ij,indextile,ngr))
          thp_opt_tile(ij,k)   = basic%thp(k,indicei(ij,indextile,ngr),    &
               indicej(ij,indextile,ngr))
          pi0_opt_tile(ij,k)   = basic%pi0(k,indicei(ij,indextile,ngr),    &
               indicej(ij,indextile,ngr))
          rtp_opt_tile(ij,k)   = basic%rtp(k,indicei(ij,indextile,ngr),    &
               indicej(ij,indextile,ngr))
          theta_opt_tile(ij,k) = basic%theta(k,indicei(ij,indextile,ngr),  &
               indicej(ij,indextile,ngr))
          rv_opt_tile(ij,k)    = basic%rv(k,indicei(ij,indextile,ngr),     &
               indicej(ij,indextile,ngr))
          ! Checando em each_column
          dn0_opt_tile(ij,k)   = basic%dn0(k,indicei(ij,indextile,ngr),    &
               indicej(ij,indextile,ngr))
          ! cldnuc
          wp_opt_tile(ij,k)    = basic%wp(k,indicei(ij,indextile,ngr),     &
               indicej(ij,indextile,ngr))

          if (level >= 2 ) then
             rcp_tile(ij,k)   = micro%rcp(k,indicei(ij,indextile,ngr),     &
                  indicej(ij,indextile,ngr))
          endif
          if (level >= 3) then
             if(irain >= 1) then
                rrp_tile(ij,k)   = micro%rrp(k,indicei(ij,indextile,ngr),  &
                     indicej(ij,indextile,ngr))
                q2_tile(ij,k)    = micro%q2(k,indicei(ij,indextile,ngr),   &
                     indicej(ij,indextile,ngr))
             endif
             if(ipris >= 1)  then
                rpp_tile(ij,k)   = micro%rpp(k,indicei(ij,indextile,ngr),  &
                     indicej(ij,indextile,ngr))
             endif
             if(isnow >= 1)  then
                rsp_tile(ij,k)   = micro%rsp(k,indicei(ij,indextile,ngr),  &
                     indicej(ij,indextile,ngr))
             endif
             if(iaggr >= 1)  then
                rap_tile(ij,k)   = micro%rap(k,indicei(ij,indextile,ngr),  &
                     indicej(ij,indextile,ngr))
             endif
             if(igraup >= 1) then
                rgp_tile(ij,k)   = micro%rgp(k,indicei(ij,indextile,ngr),  &
                     indicej(ij,indextile,ngr))
                q6_tile(ij,k)    = micro%q6(k,indicei(ij,indextile,ngr),   &
                     indicej(ij,indextile,ngr))
                rhp_tile(ij,k)   = micro%rhp(k,indicei(ij,indextile,ngr),  &
                     indicej(ij,indextile,ngr))
                q7_tile(ij,k)    = micro%q7(k,indicei(ij,indextile,ngr),   &
                     indicej(ij,indextile,ngr))
             endif
             if(icloud == 5) then
                ccp_tile(ij,k)   = micro%ccp(k,indicei(ij,indextile,ngr),  &
                     indicej(ij,indextile,ngr))
             endif
             if(irain == 5) then
                crp_tile(ij,k)   = micro%crp(k,indicei(ij,indextile,ngr),  &
                     indicej(ij,indextile,ngr))
             endif
             if(ipris == 5) then
                cpp_tile(ij,k)   = micro%cpp(k,indicei(ij,indextile,ngr),  &
                     indicej(ij,indextile,ngr))
             endif
             if(isnow == 5) then
                csp_tile(ij,k)   = micro%csp(k,indicei(ij,indextile,ngr),  &
                     indicej(ij,indextile,ngr))
             endif
             if(iaggr == 5) then
                cap_tile(ij,k)   = micro%cap(k,indicei(ij,indextile,ngr),  &
                     indicej(ij,indextile,ngr))
             endif
             if(igraup == 5) then
                cgp_tile(ij,k)   = micro%cgp(k,indicei(ij,indextile,ngr),  &
                     indicej(ij,indextile,ngr))
             endif
             if(ihail == 5) then
                chp_tile(ij,k)   = micro%chp(k,indicei(ij,indextile,ngr),  &
                     indicej(ij,indextile,ngr))
             endif
             cccnp_tile(ij,k)   = micro%cccnp(k,indicei(ij,indextile,ngr),  &
                  indicej(ij,indextile,ngr))
             cifnp_tile(ij,k)   = micro%cifnp(k,indicei(ij,indextile,ngr),  &
                  indicej(ij,indextile,ngr))
          endif

       enddo

    enddo

    ! *** 3D Arrays (ij,10):

    do k = 1,10
       do ij=1,ij_final

          k1_mod(ij,k) = 0
          k2_mod(ij,k) = 0
          k3_mod(ij,k) = 0

       enddo
    enddo

    ! *** 3D Arrays:
    do k = 1,m1
       pitot_mod(:,k)    = pitot(k)    !(OUT)
       press_mod(:,k)    = press(k)    !(OUT)
       tair_mod(:,k)     = 0.
       til_mod(:,k)      = 0.
       rliq_mod(:,k)     = 0.
       rice_mod(:,k)     = 0.
       qhydm_mod(:,k)    = 0.
       rvstr_mod(:,k)    = 0.
       tairstrc_mod(:,k) = 0.
       vapdif_mod(:,k)   = 0.
       rdynvsci_mod(:,k) = 0.
       denfac_mod(:,k) = 0.
       colfacr_mod(:,k) = 0.
       colfacr2_mod(:,k) = 0.
       colfacc_mod(:,k) = 0.
       colfacc2_mod(:,k) = 0.
       sumuy_mod(:,k) = 0.
       sumuz_mod(:,k) = 0.
       sumvr_mod(:,k) = 0.
       rvs0_mod(:,k) = 0.
       dn0i_mod(:,k) = dn0i(k)
       scrmic1_mod(:,k) = 0.
       scrmic2_mod(:,k) = 0.
       scrmic3_mod(:,k) = 0.

    enddo

    ! *** 4D Arrays:
    do cat = 1, ncat
       do k = 1,m1

          ! Input data used in thrmstr_opt
          cx_mod(:,k,cat) = cx(k,cat) !(IN)

          rx_mod(:,k,cat) = rx(k,cat) !(IN)

          qx_mod(:,k,cat) = 0 !(INOUT)

          ! enemb, each_call
          emb_mod(:,k,cat)= emb(k,cat)

          vterm_mod(:,k,cat) = 0.;

          jhcat_mod(:,k,cat) = jhcat(k,cat)

          sh_mod(:,k,cat) = sh(k,cat)

          emb_mod(:,k,cat) = emb(k,cat)

          vap_mod(:,k,cat) = 0.

          sd_mod(:,k,cat) = 0.
          se_mod(:,k,cat) = 0.
          sf_mod(:,k,cat) = 0.
          sg_mod(:,k,cat) = 0.
          sm_mod(:,k,cat) = sm(k,cat)
          ss_mod(:,k,cat) = 0.
          su_mod(:,k,cat) = 0.
          sw_mod(:,k,cat) = 0.
          sy_mod(:,k,cat) = 0.
          sz_mod(:,k,cat) = 0.
          ict1_mod(:,k,cat) = 0
          ict2_mod(:,k,cat) = 0
          wct1_mod(:,k,cat) = 0.
          wct2_mod(:,k,cat) = 0.

       enddo
    enddo

    ! *** 5D Arrays:
    rxfer_mod(:,:,:,:)  = 0.
    qrxfer_mod(:,:,:,:) = 0.
    enxfer_mod(:,:,:,:) = 0

    do k = 1, m1
       eff_mod(:,k,1) = 1.0
    enddo

    do cat = 2, 10
       do k = 1, m1
          eff_mod(:,k,cat) = 0.
       enddo
    enddo

    do cat = 1, 9
       do k = 1,m1
          ! Input data used in thrmstr_opt
          sa_mod(:,k,cat) = 0.
       enddo
    enddo

    do cat = 1, 2
       do k = 1,m1
          tref_mod(:,k,cat) = 0.
          rvsref_mod(:,k,cat) = 0.
          rvsrefp_mod(:,k,cat) = 0.
       enddo
    enddo

    ! movi de mic_driv_new

    qx_mod(:,:,:) = 0.

  end subroutine init_copy_mic_block

  !****************************************************************************

  subroutine mcphys_opt(m1, ijm, step, ngr, maxnzp, nembfall, maxkfall, mynum,  &
       dtlt, dtlti, time, zm, zt, radiate, pcpfillc, pcpfillr, sfcpcp, glat,    &
       topt, if_adap)

    use mem_radiate, only: &
         radiate_vars,     & ! INTENT(IN)
         iswrtyp,          & ! INTENT(IN)
         ilwrtyp,          & ! INTENT(IN)
         radfrq              ! INTENT(IN)

    use micphys, only:     &
         nhcat,            & ! INTENT(IN)
         jnmb,             & ! INTENT(IN)
         ncat,             & ! INTENT(IN)
         scrmic1,          & ! INTENT(INOUT)
         scrmic2,          & ! INTENT(INOUT)
         tairc,            & ! INTENT(OUT)
         ch1,              & ! INTENT(OUT)
         cfvt,             & ! INTENT(IN)
         scrmic3             ! INTENT(OUT)

    use node_mod, only :  &
         mzp,             & ! INTENT(IN)
         mxp,             & ! INTENT(IN)
         myp                ! INTENT(IN)

    implicit none

    ! Arguments:
    integer, intent(in)                   :: m1, ijm, step, ngr, maxnzp, & 
         nembfall, maxkfall, mynum
    real, intent(in)                      :: dtlt, dtlti
    real(kind=8), intent(in)              :: time
    real, dimension(m1), intent(in)       :: zm
    type (radiate_vars), intent(inout)    :: radiate
    integer, intent(in)                   :: if_adap
    ! Variables needed for Harrington radiation scheme
    real, dimension(m1,maxkfall,nembfall,nhcat), intent(in) :: &
         pcpfillc, pcpfillr
    real, dimension(maxkfall,nembfall,nhcat), intent(in)    :: sfcpcp
    real, intent(in)                                        :: &
         glat(:,:), topt(:,:)
    real, dimension(m1), intent(in)                         :: zt

    ! Local Variables:
    integer :: ij
    integer :: k,jflag,lcat,icv,icx,mc1,mc2,mc3,mc4  &
         ,mcat,lhcat
    integer, dimension(ij_final) :: k1cnuc, k2cnuc
    integer, dimension(7), parameter :: mivap  = (/1,3,4,5,2,6, 7/)
    integer, dimension(7), parameter :: mcats  = (/0,3,0,0,6,7,10/)
    integer, dimension(4), parameter :: mcat33 = (/0,0,4,5/)

    integer, dimension(9,4), parameter :: mcat1 = reshape ( &
         (/3,3,3,4,4,4,5,5,6,  &
         5,6,7,5,6,7,6,7,7,    &
         5,6,7,5,6,7,6,7,7,    &
         4,7,8,5,7,8,7,8,8/),  &
         (/ 9, 4/) ) 
    integer, dimension(7,2), parameter :: mcat2 = reshape ( &
         (/0,0,0,6,6,7,7,  &
         0,0,0,2,2,9,9/),  &
         (/ 7, 2/) ) 

    integer, dimension(7), parameter   :: mix02 = (/3,1,4,5,6,7,2/)

    real, dimension(7), parameter      :: dpcp0 = &
         (/.001,.001,.010,.010,.010,.003,.001/)

    ! Local variables necessary in radcalc3
    integer :: local_lpw, i, j
    real    :: local_rtgt, local_rv(m1), local_dn0(m1)

    ! Local variables necessary in sedim_opt
    real    :: local_ch1(ij_final,nhcat)
    real    :: time_rfrq

    !RETIRAR DEPOIS
    sm_mod(:,:,3:7) = 0. !OK
    sh_mod(:,:,3:5) = 0. ! OK 
    eff_mod(:,:,2:) = 0.


    call thrmstr_opt(ngr, m1)

    call each_column_opt(ngr, m1)

    ! Diagnose hydrometeor mean mass emb, and if necessary, number concentration.

    jflag = 1

    do lcat = 1,7
       if (jnmb(lcat) .ge. 1) then
          call enemb_opt(ngr, lcat, jflag)
       endif
    enddo

    ! Evaluate radiative heating rates if using Harrington radiation scheme

!    if (iswrtyp==3 .or. ilwrtyp==3) then
!
!       time_rfrq = real(dmod(time,dble(radfrq)))
!
!       if (mod(time_rfrq + .001,radfrq)<dtlt .or. time<.001) then
!
!          ! Using Harrington radiation scheme
!          ! Loop through 'ij' vector
!          ! Converting data 'ij,k' to 'k,i,j'
!
!          do ij=1, ij_final
!
!             i = indicei(ij, indextile, ngr)
!             j = indicej(ij, indextile, ngr)
!
!             local_lpw = lpw_opt(ij)
!             local_rtgt = rtgt_opt_tile(ij)
!             call copy_ijk2k(1, m1, ij, rv_opt_tile, local_rv)
!             call copy_ijk2k(1, m1, ij, dn0_opt_tile, local_dn0)
!
!             call radcalc3(m1,maxnzp,ncat,iswrtyp,ilwrtyp,if_adap,local_lpw,    &
!                  glat, local_rtgt, topt,                                       &
!                  radiate%albedt(i,j), radiate%cosz(i,j),                       &
!                  radiate%rlongup(i,j), radiate%rshort(i,j),                    &
!                  radiate%rlong(i,j),                                           &
!                  zm, zt, local_rv(1), local_dn0(1), radiate%fthrd(1,i,j), i, j,&
!                  time, ngr)
!
!          enddo
!
!       endif
!    endif

    do lcat = 1,7
       if (jnmb(lcat) .ge. 1) then
          call diffprep_opt(ngr, lcat) 
       endif
    enddo

    call vapdiff_opt(ngr, 10)

    do icv = 1,7
       lcat = mivap(icv)

       if (jnmb(lcat) .ge. 1) then
          call vapflux_opt(ngr, lcat)
       endif
    enddo

    if (jnmb(4) .ge. 1) then
       call psxfer_opt()
    endif

    jflag = 2
    do lcat = 1,7
       if (jnmb(lcat) .ge. 1) then
          call enemb_opt(ngr, lcat, jflag)
          call getict_opt(lcat)
       endif
    enddo

    call newtemp_opt(ngr)

    if (jnmb(2) .ge. 1) then
       call auto_accret_opt(ngr, dtlt)
    endif

    call effxy_opt(m1)

    ! Self collection of rain, aggregates, graupel, hail:  number change only

    do lcat = 2,7

       if (lcat==3 .or. lcat==4) cycle
       mc1 = mcats(lcat)
       if (jnmb(lcat)>=5) then
          call cols_opt(lcat, mc1, lcat)
       endif

    enddo

    ! Self collection of pristine ice, snow

    do lcat = 3,4
       mc1 = mcat33(lcat)
       if (jnmb(lcat)>=1 .and. jnmb(5)>=1) then
          call col3344_opt(lcat, 5, mc1, lcat)
       endif
    enddo

    ! Collection between pristine ice and snow

    if (jnmb(5) .ge. 1) then
       call col3443_opt(3, 4, 5)
    endif

    ! Ice-ice collisions

    do icx = 1,9
       mc1 = mcat1(icx,1)
       mc2 = mcat1(icx,2)
       mc3 = mcat1(icx,3)
       mc4 = mcat1(icx,4)

       if (jnmb(mc1)>=1 .and. jnmb(mc3)>=1) then
          call col1_opt(mc1, mc2, mc3, mc4, mc1, mc2)
       endif
    enddo

    ! Ice-cloud collisions

    do lcat = 4,7
       mc1 = mcat2(lcat,1)
       mc2 = mcat2(lcat,2)

       if (jnmb(lcat)>=1 .and. jnmb(mc1)>=1) then
          call col2_opt    (ngr,1,lcat,mc1,mc2,dtlt,1,lcat)
       endif
    enddo

    ! Ice-rain collisions

    do lcat = 3,7
       if (jnmb(lcat)>=1 .and. jnmb(7)>=1) then
          call col3_opt(2, lcat, 7, 2, lcat)
       endif
    enddo

    call colxfers_opt()

    do mcat = 1,7
       lcat = mix02(mcat)
       if (jnmb(lcat)>=1) &
            call x02_opt(ngr, lcat)
    enddo

    !
    k1cnuc(:) = 0.
    k2cnuc(:) = 0.
    !

    if (jnmb(1)>=1) then
       call cldnuc_opt(ngr,m1,k1cnuc,k2cnuc)
    endif

    do ij = 1, ij_final
       k1_mod(ij,1) = min(k1_mod(ij,1), k1cnuc(ij))
       k2_mod(ij,1) = max(k2_mod(ij,1), k2cnuc(ij))
       k3_mod(ij,1) = max(k2_mod(ij,1), &
            k3_mod(ij,1))
    enddo

    if (jnmb(1) .ge. 1) then
       call c03_opt(ngr, 1)
    endif

    k1cnuc(:) = 0.
    k2cnuc(:) = 0.

    if (jnmb(3)>=1) then
       call icenuc_opt(ngr,m1,k1cnuc,k2cnuc,dtlt)
    endif

    do ij=1, ij_final
       k1_mod(ij,3) = min(k1_mod(ij,3), k1cnuc(ij))
       k2_mod(ij,3) = max(k2_mod(ij,3), k2cnuc(ij))
       k3_mod(ij,3) = max(k2_mod(ij,3), &
            k3_mod(ij,3))
    enddo

    do lcat = 3, 1, -2
       if (jnmb(lcat)>=1) then
          call pc03_opt(ngr, lcat)
       endif
    enddo

    !  Zero out precip arrays.
    pcpg_opt_tile(:)  = 0.
    qpcpg_opt_tile(:) = 0.
    dpcpg_opt_tile(:) = 0.

    ! tairc is used here to accumulate changes to thp from sedim

    do k = minval(lpw_opt(:)), m1
       do ij = 1, ij_final
          if (k>=lpw_opt(ij)) tairc_mod(ij,k) = 0.
       enddo
    enddo

    !loop inserted inside sedim_opt for lhcat= 2,7
    do lhcat = 2,nhcat
       ch1(lhcat) = dtlt * cfvt(lhcat) / rtgt_opt_tile(ij_final)
    enddo


    ! New stp for calculating array: local_ch1
    ! coping values from ch1
    do lhcat = 1, nhcat
       local_ch1(:, lhcat) = ch1(lhcat)
    enddo
    !local_ch1(:, 1) = ch1(1)
    ! calculating new values
    do lhcat = 2, nhcat
       do ij = 1, ij_final
          local_ch1(ij, lhcat) = dtlt*cfvt(lhcat)/rtgt_opt_tile(ij)
       enddo
    enddo

    do lcat = 2,7
       if (jnmb(lcat)>=1) then
          call sedim_opt(m1, lcat, ngr, nembfall, maxkfall, dpcp0(lcat),       &
               dtlti, pcpfillc,pcpfillr,sfcpcp, dtlt, local_ch1, ij_final,     &
               pcprx_mod,                                         &
               k1_mod, k2_mod,                       &
               scrmic1_mod, scrmic2_mod,             &
               scrmic3_mod,                                       &
               jhcat_mod,                                         &
               rx_mod, cx_mod, qx_mod,  &
               emb_mod, dn0i_mod,                    &
               accpx_mod, tairc_mod,                 &
               tair_mod, dn0_opt_tile,                &
               rtgt_opt_tile, pcpg_opt_tile,           &
               qpcpg_opt_tile, dpcpg_opt_tile,         &
               rtp_opt_tile, thp_opt_tile,             &
               theta_opt_tile)
       endif
    enddo

    do k = minval(lpw_opt(:)), m1
       do ij =1, ij_final
          if (k>=lpw_opt(ij)) then
             thp_opt_tile(ij,k) = thp_opt_tile(ij,k) + &
                  tairc_mod(ij,k)
          endif
       enddo
    enddo

    return
  end subroutine mcphys_opt

  ! **************************************************************************

  subroutine thrmstr_opt(ngr,m1)

    use rconstants, only: &
         p00,             & !INTENT(IN)
         cpi,             & !INTENT(IN)
         cpor,            & !INTENT(IN)
         alvl,            & !INTENT(IN)
         alvi,            & !INTENT(IN)
         cpi4,            & !INTENT(IN)
         cp253i             !INTENT(IN)

    implicit none

    ! Arguments:
    integer, intent(in)                :: ngr,m1

    ! Local Variables:
    integer :: k,lcat
    real    :: fracliq,tcoal,tairstr
    integer :: ij

    do k = minval(lpw_opt),m1
       do ij=1,ij_final
          if (k>=lpw_opt(ij)) then
             pitot_mod(ij,k) = pi0_opt_tile(ij,k) + pp_opt_tile(ij,k)
             press_mod(ij,k) = p00 * (pitot_mod(ij,k) * cpi) ** cpor
             tair_mod(ij,k)  = theta_opt_tile(ij,k)*pitot_mod(ij,k) * cpi
          endif
       enddo
    enddo


    do k = 1, maxval(k1_mod(:,10))-1  !1,k1(10)-1
       do ij=1,ij_final
          if (k <= (k1_mod(ij,10)-1)) then
             theta_opt_tile(ij,k) = thp_opt_tile(ij,k)
             rv_opt_tile(ij,k)    = rtp_opt_tile(ij,k)
          endif
       enddo
    enddo


    do k = minval(k2_mod(:,10))+1,m1 !k2(10)+1,m1
       do ij=1,ij_final
          if (k >= (k2_mod(ij,10)+1)) then
             theta_opt_tile(ij,k) = thp_opt_tile(ij,k)
             rv_opt_tile(ij,k)    = rtp_opt_tile(ij,k)
          endif
       enddo
    enddo


    do k = minval(k1_mod(:,10)),maxval(k2_mod(:,10))
       !k1(10),k2(10)
       do ij=1,ij_final
          if ((k>=k1_mod(ij,10)).and.(k<=k2_mod(ij,10))) then
             til_mod(ij,k) = thp_opt_tile(ij,k) * pitot_mod(ij,k) * cpi
             rliq_mod(ij,k) = 0.
             rice_mod(ij,k) = 0.
          endif
       enddo
    enddo


    do lcat = 1,2
       do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))
          !k1(lcat),k2(lcat)
          do ij=1,ij_final
             if ((k>=k1_mod(ij,lcat)).and.(k<=k2_mod(ij,lcat))) then
                rliq_mod(ij,k) = rliq_mod(ij,k) + rx_mod(ij,k,lcat)
             endif
          enddo
       enddo
    enddo


    do lcat = 3,5
       do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))
          !k1(lcat),k2(lcat)
          do ij=1,ij_final
             if ((k>=k1_mod(ij,lcat)).and.(k<=k2_mod(ij,lcat))) then
                rice_mod(ij,k) = rice_mod(ij,k) + rx_mod(ij,k,lcat)
             endif
          enddo
       enddo
    enddo


    do lcat = 6,7
       do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))
          !k1(lcat),k2(lcat)
          do ij=1,ij_final
             if ((k>=k1_mod(ij,lcat)).and.(k<=k2_mod(ij,lcat))) then
                call qtc(qx_mod(ij,k,lcat),tcoal,fracliq)
                rliq_mod(ij,k) = rliq_mod(ij,k) + rx_mod(ij,k,lcat)*fracliq
                rice_mod(ij,k) = rice_mod(ij,k) + rx_mod(ij,k,lcat)*(1.-fracliq)
             endif
          enddo
       enddo
    enddo


    do k = minval(k1_mod(:,10)),maxval(k2_mod(:,10))
       !k1(10),k2(10)
       do ij=1,ij_final
          if ((k>=k1_mod(ij,10)).and.(k<=k2_mod(ij,10))) then
             qhydm_mod(ij,k) = alvl * rliq_mod(ij,k) + alvi * rice_mod(ij,k)
             rvstr_mod(ij,k) = rtp_opt_tile(ij,k) - rliq_mod(ij,k) - &
                  rice_mod(ij,k)
             sa_mod(ij,k,1) = til_mod(ij,k) * qhydm_mod(ij,k) / &
                  (1.e-12 + rliq_mod(ij,k) + rice_mod(ij,k))
          endif
       enddo
    enddo


    do k = minval(k1_mod(:,10)),maxval(k2_mod(:,10))
       !k1(10),k2(10)
       do ij=1,ij_final
          if ((k>=k1_mod(ij,10)).and.(k<=k2_mod(ij,10))) then
             if (tair_mod(ij,k) > 253.) then
                tairstr = 0.5 * (til_mod(ij,k) + sqrt(til_mod(ij,k) * &
                     (til_mod(ij,k) + cpi4 * qhydm_mod(ij,k))))
                sa_mod(ij,k,1)  = sa_mod(ij,k,1) * cpi / &
                     (2. * tairstr - til_mod(ij,k))
             else
                tairstr = til_mod(ij,k) * (1. + qhydm_mod(ij,k) * cp253i)
                sa_mod(ij,k,1)  = sa_mod(ij,k,1) * cp253i
             endif
             tairstrc_mod(ij,k) = tairstr - 273.16
          endif
       enddo
    enddo

    return
  end subroutine thrmstr_opt

  !***************************************************************************

  subroutine each_column_opt(ngr, m1)

    use rconstants, only : &
         alvl,             & ! INTENT(IN)
         alvi                ! INTENT(IN)

    use micphys, only :    &
         jhabtab,          & ! INTENT(IN)
         colf                ! INTENT(IN)

    implicit none

    ! Arguments
    integer, intent(in)                     :: ngr, m1

    ! Local Variables
    integer         :: k, nt, ns
    real            :: elsref, elsrefp, dplinv, eisref, eisrefp, dpiinv, relhum
    real, parameter :: ck1=-4.818544e-3, ck2= 1.407892e-4, ck3=-1.249986e-7
    integer         :: ij

    ! External Functions
    real :: rslf,rsif,eslf,eslpf,esif,esipf

    do k = minval(lpw_opt),m1-1
       do ij=1,ij_final
          if (k >= lpw_opt(ij)) then
             rvlsair_mod(ij,k) = rslf(press_mod(ij,k), tair_mod(ij,k))

             rvisair_mod(ij,k) = rsif(press_mod(ij,k), tair_mod(ij,k))

             dn0i_mod(ij,k)    = 1. / dn0_opt_tile(ij,k)

             tairc_mod(ij,k)   = tair_mod(ij,k) - 273.16

             tx_mod(ij,k,1)    = tairc_mod(ij,k)

             thrmcon_mod(ij,k) = ck1 + (ck2 + ck3 * tair_mod(ij,k)) * &
                  tair_mod(ij,k)

             dynvisc_mod(ij,k) = .1718e-4 + .49e-7 * tairc_mod(ij,k)

             ! Diagnose habit of pristine ice and snow

             nt                = max(1, min(31, -nint(tairc_mod(ij,k))))
             relhum            = min(1., rv_opt_tile(ij,k) / rvlsair_mod(ij,k))

             ns                = max(1, nint(100. * relhum))

             jhcat_mod(ij,k,3) = jhabtab(nt,ns,1)
             jhcat_mod(ij,k,4) = jhabtab(nt,ns,2)
          endif
       enddo
    enddo

    do k = minval(k1_mod(:,10)),maxval(k2_mod(:,10))
       do ij=1,ij_final
          if ((k>=k1_mod(ij,10)).and.(k<=k2_mod(ij,10))) then
             vapdif_mod(ij,k)   = 2.14 * (tair_mod(ij,k) / 273.15) ** 1.94 / &
                  press_mod(ij,k)
             rdynvsci_mod(ij,k) = sqrt(1. / dynvisc_mod(ij,k))
             denfac_mod(ij,k)   = sqrt(dn0i_mod(ij,k))

             colfacr_mod(ij,k)  = colf * denfac_mod(ij,k) * dn0_opt_tile(ij,k)
             colfacr2_mod(ij,k) = 2. * colfacr_mod(ij,k)
             colfacc_mod(ij,k)  = colfacr_mod(ij,k) * dn0_opt_tile(ij,k)
             colfacc2_mod(ij,k) = 2. * colfacc_mod(ij,k)

             tref_mod(ij,k,1)   = tairc_mod(ij,k) - &
                  min(25., 700. * (rvlsair_mod(ij,k) - rv_opt_tile(ij,k)))
             sa_mod(ij,k,2) = thrmcon_mod(ij,k) * sa_mod(ij,k,1)
             sa_mod(ij,k,3) = thrmcon_mod(ij,k) * (tairstrc_mod(ij,k) + &
                  sa_mod(ij,k,1) * rvstr_mod(ij,k))

             sumuy_mod(ij,k) = 0.
             sumuz_mod(ij,k) = 0.
             sumvr_mod(ij,k) = 0.
          endif
       enddo
    enddo

    do k = minval(k1_mod(:,8)), maxval(k2_mod(:,8))
       !k = k1(8),k2(8)
       do ij=1,ij_final
          if ((k>=k1_mod(ij,8)).and.(k<=k2_mod(ij,8))) then
             elsref  = eslf(tref_mod(ij,k,1))
             elsrefp = eslpf(tref_mod(ij,k,1))
             dplinv  = 1. / (press_mod(ij,k) - elsref)
             rvsref_mod(ij,k,1)  = .622 * elsref * dplinv
             rvsrefp_mod(ij,k,1) = .622 * elsrefp * dplinv * &
                  (1. + elsref * dplinv)

             sa_mod(ij,k,4)      = rvsrefp_mod(ij,k,1) * tref_mod(ij,k,1) - &
                  rvsref_mod(ij,k,1)
             sa_mod(ij,k,6)      = alvl * rvsrefp_mod(ij,k,1)
             sa_mod(ij,k,8)      = alvl * sa_mod(ij,k,4)
          endif
       enddo
    enddo

    do k = minval(k1_mod(:,9)), maxval(k2_mod(:,9))
       !k = k1(9),k2(9)
       do ij=1,ij_final
          if ((k>=k1_mod(ij,9)).and.(k<=k2_mod(ij,9))) then
             tref_mod(ij,k,2) = min(0., tref_mod(ij,k,1))
             eisref  = esif (tref_mod(ij,k,2))
             eisrefp = esipf(tref_mod(ij,k,2))
             dpiinv  = 1. / (press_mod(ij,k) - eisref)
             rvsref_mod (ij,k,2) = .622 * eisref  * dpiinv
             rvsrefp_mod(ij,k,2) = .622 * eisrefp * dpiinv * &
                  (1. + eisref * dpiinv)
             rvs0_mod(ij,k)      = 379.4 / (press_mod(ij,k) - 610.)

             sa_mod(ij,k,5) = rvsrefp_mod(ij,k,2) * tref_mod(ij,k,2) - &
                  rvsref_mod(ij,k,2)
             sa_mod(ij,k,7) = alvi * rvsrefp_mod(ij,k,2)
             sa_mod(ij,k,9) = alvi * sa_mod(ij,k,5)
             sh_mod(ij,k,3) = 0.
             sh_mod(ij,k,4) = 0.
             sh_mod(ij,k,5) = 0.
          endif
       enddo
    enddo

    return
  end subroutine each_column_opt

  !***************************************************************************

  subroutine enemb_opt(ngr, lcat, jflag)

    use micphys, only : &
         jnmb,          & ! INTENT(IN)
         cfemb0,        & ! INTENT(IN)
         pwemb0,        & ! INTENT(IN)
         cfen0,         & ! INTENT(IN)
         pwen0,         & ! INTENT(IN)
         parm,          & ! INTENT(IN)
         emb0,          & ! INTENT(IN)
         emb1,          & ! INTENT(IN)
         enmlttab         ! INTENT(IN)

    implicit none

    ! Arguments
    integer, intent(in)             :: ngr, lcat, jflag

    ! Local Variables
    integer :: k,lhcat, ij
    real :: embi(ij_dimension),parmi,fracmass,cxloss


    if (jnmb(lcat) == 2) then
       do ij =1, ij_final
          embi(ij) = 1. / emb_mod(ij,2,lcat)
       enddo
       do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))        !k1,k2
          do ij =1, ij_final
             if (k>=k1_mod(ij,lcat).and.k<=k2_mod(ij,lcat)) &
                  cx_mod(ij,k,lcat) = rx_mod(ij,k,lcat) * embi(ij)
          enddo
       enddo
    elseif (jnmb(lcat) == 3) then
       do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))        !k1,k2
          do ij =1, ij_final
             if (k>=k1_mod(ij,lcat).and. k<=k2_mod(ij,lcat)) then
                lhcat = jhcat_mod(ij,k,lcat)
                emb_mod(ij,k,lcat) = cfemb0(lhcat) * (dn0_opt_tile(ij,k) * &
                     rx_mod(ij,k,lcat)) ** pwemb0(lhcat)
                cx_mod(ij,k,lcat) = cfen0(lhcat) * dn0i_mod(ij,k) * &
                     (dn0_opt_tile(ij,k) * rx_mod(ij,k,lcat)) ** pwen0(lhcat)
             endif
          enddo
       enddo
    elseif (jnmb(lcat) == 4) then
       parmi = 1. / parm(lcat)
       do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))         !k1,k2
          do ij =1, ij_final
             if (k>=k1_mod(ij,lcat).and.k<=k2_mod(ij,lcat)) then
                emb_mod(ij,k,lcat) = max(emb0(lcat), &
                     min(emb1(lcat),rx_mod(ij,k,lcat) * parmi))
                cx_mod(ij,k,lcat) = rx_mod(ij,k,lcat) / emb_mod(ij,k,lcat)
             endif
          enddo
       enddo
    elseif (jnmb(lcat) >= 5 .and. jflag == 1) then
       do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))         !k1,k2
          do ij =1, ij_final
             if (k>=k1_mod(ij,lcat).and. k<=k2_mod(ij,lcat))        then
                emb_mod(ij,k,lcat) = max(emb0(lcat), &
                     min(emb1(lcat),rx_mod(ij,k,lcat)/max(1.e-9,cx_mod(ij,k,lcat))))
                cx_mod(ij,k,lcat) = rx_mod(ij,k,lcat) / emb_mod(ij,k,lcat)
             endif
          enddo
       enddo
    elseif (jnmb(lcat) >= 5 .and. jflag == 2) then
       do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))         !k1,k2
          do ij =1, ij_final
             if (k>=k1_mod(ij,lcat).and.k<=k2_mod(ij,lcat)) then
                if (rx_mod(ij,k,lcat) >= 1.e-9) then
                   if (vap_mod(ij,k,lcat) < 0.) then
                      fracmass = min(1.,-vap_mod(ij,k,lcat) / rx_mod(ij,k,lcat))
                      cxloss = cx_mod(ij,k,lcat)*enmlttab(int(200.*fracmass)+1,&
                           jhcat_mod(ij,k,lcat))
                      cx_mod(ij,k,lcat) = cx_mod(ij,k,lcat) - cxloss
                   endif
                   emb_mod(ij,k,lcat) = max(emb0(lcat), &
                        min(emb1(lcat),rx_mod(ij,k,lcat)/max(1.e-9,cx_mod(ij,k,lcat))))
                   cx_mod(ij,k,lcat) = rx_mod(ij,k,lcat) / emb_mod(ij,k,lcat)
                endif
             endif
          enddo
       enddo
    endif

    return
  end subroutine enemb_opt

  !******************************************************************************

  subroutine diffprep_opt(ngr, lcat)

    use micphys, only: &
         ncat,         & !INTENT(IN)
         frefac1,      & !INTENT(IN)
         pwmasi,       & !INTENT(IN)
         frefac2,      & !INTENT(IN)
         cdp1,         & !INTENT(IN)
         pi4dt,        & !INTENT(IN)
         sl,           & !INTENT(IN)
         sj,           & !INTENT(IN)
         sc,           & !INTENT(IN)
         sk,           & !INTENT(IN)
         sm              !INTENT(INOUT)

    use grid_dims, only: nzpmax !INTENT(IN)


    implicit none

    ! Arguments:
    integer, intent(in)             :: ngr,lcat

    ! Local Variables:
    integer :: k,if1,if4,if6,if8, ij, cat !,lhcat
    ! For optimizations
    real :: fre_opt(ij_final), sb_opt(ij_final,nzpmax,ncat), &
         scdei_opt(ij_final), ttest_opt(ij_final,nzpmax,ncat)

    ! Putting Zero on local arrays
    do ij =1, ij_final
       fre_opt(ij)   = 0.
       scdei_opt(ij) = 0.
    enddo
    do cat = 1, ncat
       do k = 1, nzpmax
          do ij =1, ij_final
             sb_opt(ij,k,cat)    = 0.
             ttest_opt(ij,k,cat) = 0.
          enddo
       enddo
    enddo

    if (lcat .le. 2) then
       if1 = 1
       if4 = 4
       if6 = 6
       if8 = 8
    else
       if1 = 2
       if4 = 5
       if6 = 7
       if8 = 9
    endif

    do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))         !k1,k2

       do ij =1, ij_final

          if ((k>=k1_mod(ij,lcat)).and.(k<=k2_mod(ij,lcat))) then

             !lhcat = jhcat_mod(ij,k,lcat)

             if (rx_mod(ij,k,lcat) < 1.e-9) go to 229
             !if (rx_mod(ij,k,lcat) >= 1.e-9) then
             fre_opt(ij) = frefac1(jhcat_mod(ij,k,lcat)) * &
                  emb_mod(ij,k,lcat) ** pwmasi(jhcat_mod(ij,k,lcat)) + &
                  rdynvsci_mod(ij,k) * frefac2(jhcat_mod(ij,k,lcat)) * &
                  emb_mod(ij,k,lcat) ** cdp1(jhcat_mod(ij,k,lcat))
             sb_opt(ij,k,lcat) = cx_mod(ij,k,lcat) * dn0_opt_tile(ij,k) * &
                  fre_opt(ij) * pi4dt
             su_mod(ij,k,lcat) = vapdif_mod(ij,k) * sb_opt(ij,k,lcat)
             sd_mod(ij,k,lcat) = sh_mod(ij,k,lcat) * rx_mod(ij,k,lcat)
             se_mod(ij,k,lcat) = su_mod(ij,k,lcat) * sa_mod(ij,k,if6) + &
                  sb_opt(ij,k,lcat) * thrmcon_mod(ij,k)
             sf_mod(ij,k,lcat) = su_mod(ij,k,lcat) * sl(if1) - &
                  sb_opt(ij,k,lcat) * sa_mod(ij,k,2)
             sg_mod(ij,k,lcat) = su_mod(ij,k,lcat) * sa_mod(ij,k,if8) + &
                  sb_opt(ij,k,lcat)*sa_mod(ij,k,3) + sj(lcat)*qr_mod(ij,k,lcat)
             !     + lambda_j 
             !     [Joules/kg_air added by radiative heating this timestep]
             scdei_opt(ij) = 1./(sc(if1)*sd_mod(ij,k,lcat) + se_mod(ij,k,lcat))
             ss_mod(ij,k,lcat) = sf_mod(ij,k,lcat) * scdei_opt(ij)
             sw_mod(ij,k,lcat) = (sg_mod(ij,k,lcat)-sk(if1)*sd_mod(ij,k,lcat)) *&
                  scdei_opt(ij)
             ttest_opt(ij,k,lcat) = ss_mod(ij,k,lcat)*rv_opt_tile(ij,k) + &
                  sw_mod(ij,k,lcat)
             !endif
229          continue

             if (lcat >= 3 .and. lcat <= 5) then
                if (rx_mod(ij,k,lcat) < 1.e-9) go to 228
                if (ttest_opt(ij,k,lcat) >= 0.) then
                   sm_mod(ij,k,lcat) = 0.
                   sh_mod(ij,k,lcat) = 1.
                   sd_mod(ij,k,lcat) = sh_mod(ij,k,lcat) * rx_mod(ij,k,lcat)
                   scdei_opt(ij) = 1./(sc(if1)*sd_mod(ij,k,lcat) + se_mod(ij,k,lcat))
                   ss_mod(ij,k,lcat) = sf_mod(ij,k,lcat) * scdei_opt(ij)
                   sw_mod(ij,k,lcat) = (sg_mod(ij,k,lcat) - sk(if1) *        &
                        sd_mod(ij,k,lcat)) * scdei_opt(ij)
                else
                   sm_mod(ij,k,lcat) = 1.
                endif
             endif
228          continue

             if (lcat >= 6) then
                if (rx_mod(ij,k,lcat) < 1.e-9) go to 227
                if (ttest_opt(ij,k,lcat) >= 0.) then
                   sm_mod(ij,k,lcat) = 0.
                else
                   sm_mod(ij,k,lcat) = 1.
                endif
             endif
227          continue

             if (rx_mod(ij,k,lcat) < 1.e-9) go to 226
             sy_mod(ij,k,lcat) = rvsrefp_mod(ij,k,if1) * sm_mod(ij,k,lcat) * &
                  sw_mod(ij,k,lcat) - sa_mod(ij,k,if4)
             sz_mod(ij,k,lcat) = 1. - rvsrefp_mod(ij,k,if1) * ss_mod(ij,k,lcat) * &
                  sm_mod(ij,k,lcat)
             sumuy_mod(ij,k) = sumuy_mod(ij,k) + su_mod(ij,k,lcat) * sy_mod(ij,k,lcat)
             sumuz_mod(ij,k) = sumuz_mod(ij,k) + su_mod(ij,k,lcat) * sz_mod(ij,k,lcat)
226          continue

          endif

       enddo

    enddo

    return
  end subroutine diffprep_opt

  !******************************************************************************

  subroutine vapdiff_opt (ngr, lcat)

    implicit none

    !Arguments:
    integer, intent(in)              :: ngr, lcat

    !Local Variables:
    integer ::k, ij

    do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))        !kf1,kf2

       do ij =1, ij_final

          if ((k>=k1_mod(ij,lcat)).and.(k<=k2_mod(ij,lcat))) then

             rv_opt_tile(ij,k) = (rvstr_mod(ij,k) + sumuy_mod(ij,k)) /     &
                  (1.0 + sumuz_mod(ij,k))

          endif

       enddo

    enddo

    return
  end subroutine vapdiff_opt

  !**************************************************************************

  subroutine vapflux_opt(ngr, lcat)

    use micphys, only: &
         sc,           & !INTENT(IN)
         sk              !INTENT(IN)

    implicit none

    ! Arguments:
    integer, intent(in)                :: ngr,lcat

    ! Local Variables:

    integer :: if1,if4,k,ij
    real :: rxx(ij_final)

    if (lcat .le. 2) then
       if1 = 1
       if4 = 4
    else
       if1 = 2
       if4 = 5
    endif

    do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat)) !k1,k2

       do ij =1, ij_final

          if ((k>=k1_mod(ij,lcat)).and. (k<=k2_mod(ij,lcat))) then

             if (rx_mod(ij,k,lcat) .lt. 1.e-9) go to 229

             tx_mod(ij,k,lcat) = (ss_mod(ij,k,lcat) * rv_opt_tile(ij,k) + &
                  sw_mod(ij,k,lcat)) * sm_mod(ij,k,lcat)
             vap_mod(ij,k,lcat) = su_mod(ij,k,lcat) * (rv_opt_tile(ij,k) + &
                  sa_mod(ij,k,if4) - rvsrefp_mod(ij,k,if1) * tx_mod(ij,k,lcat))

             if (vap_mod(ij,k,lcat) > -rx_mod(ij,k,lcat)) then

                rxx(ij) = rx_mod(ij,k,lcat) + vap_mod(ij,k,lcat)

                if (sm_mod(ij,k,lcat) .gt. .5) then
                   qx_mod(ij,k,lcat) = sc(if1) * tx_mod(ij,k,lcat) + sk(if1)
                   qr_mod(ij,k,lcat) = qx_mod(ij,k,lcat) * rxx(ij)
                else
                   qx_mod(ij,k,lcat) = (rv_opt_tile(ij,k) * sf_mod(ij,k,lcat) + &
                        sg_mod(ij,k,lcat) - tx_mod(ij,k,lcat) * &
                        se_mod(ij,k,lcat)) / sd_mod(ij,k,lcat)
                   qx_mod(ij,k,lcat) = min(350000., &
                        max(-100000.,qx_mod(ij,k,lcat)))
                   qr_mod(ij,k,lcat) = qx_mod(ij,k,lcat) * rxx(ij)
                endif

             endif

             !bob Now also do the following section if pristine ice totally 
             ! melts: evaporate it too.

             if ((lcat==3 .and. qx_mod(ij,k,lcat)>330000.) .or. &
                  vap_mod(ij,k,lcat) <= -rx_mod(ij,k,lcat)) then

                sumuy_mod(ij,k) = sumuy_mod(ij,k) - su_mod(ij,k,lcat) * &
                     sy_mod(ij,k,lcat)
                sumuz_mod(ij,k) = sumuz_mod(ij,k) - su_mod(ij,k,lcat) * &
                     sz_mod(ij,k,lcat)
                sumvr_mod(ij,k) = sumvr_mod(ij,k) + rx_mod(ij,k,lcat)
                rv_opt_tile(ij,k) = (rvstr_mod(ij,k) + sumuy_mod(ij,k) + &
                     sumvr_mod(ij,k)) / (1.0 + sumuz_mod(ij,k))

                vap_mod(ij,k,lcat) = - rx_mod(ij,k,lcat)
                tx_mod(ij,k,lcat) = 0.
                rx_mod(ij,k,lcat) = 0.
                qx_mod(ij,k,lcat) = 0.
                qr_mod(ij,k,lcat) = 0.
             else
                rx_mod(ij,k,lcat) = rxx(ij)
             endif

229          continue

          endif

       enddo

    enddo
    return
  end subroutine vapflux_opt

  !***************************************************************************

  subroutine psxfer_opt()

    use micphys, only: &
         dnfac,        & !INTENT(IN)
         pwmasi,       & !INTENT(IN)
         gam,          & !INTENT(IN)
         dps2,         & !INTENT(IN)
         dps,          & !INTENT(IN)
         gnu,          & !INTENT(IN)
         gamn1,        & !INTENT(IN)
         pwmas,        & !INTENT(IN)
         dpsmi           !INTENT(IN)

    implicit none

    !Arguments:

    !Local Variables:
    integer :: k,lhcat,it, ij
    real :: embx,dn,xlim,dvap,dqr,dnum

    do k = min(minval(k1_mod(:,3)),                             &
         minval(k1_mod(:,4))) ,                                 &
         max(maxval(k2_mod(:,3)),                               &
         maxval(k2_mod(:,4)))
       !k1,k2 
       !min(k1(3),k1(4)),max(k2(3),k2(4))

       do ij =1, ij_final

          if ((k>=min(k1_mod(ij,3),k1_mod(ij,4))) .and. &
               (k<=max(k2_mod(ij,3),k2_mod(ij,4)))) then

             if ((vap_mod(ij,k,3) > 0.) .or. (vap_mod(ij,k,4) < 0.)) then

                if (vap_mod(ij,k,3) > 0.) then
                   lhcat = jhcat_mod(ij,k,3)
                   embx = max(1.e-9,rx_mod(ij,k,3)) / max(1.e-3,cx_mod(ij,k,3))
                   dn = dnfac(lhcat) * embx ** pwmasi(lhcat)
                   it = nint(dn * 1.e6)

                   xlim = gam(it,3) * dps2 * (dps / dn) ** (gnu(3) - 1.)  /  &
                        (gamn1(3) * pwmas(lhcat) * dn ** 2)

                   dvap = min(rx_mod(ij,k,3), vap_mod(ij,k,3) * &
                        (xlim + gam(it,1) / gamn1(3)))
                   dqr = dvap * qx_mod(ij,k,3)
                   dnum = dvap * min(dpsmi(lhcat), 1./embx)
                else
                   lhcat = jhcat_mod(ij,k,4)
                   embx = max(1.e-9,rx_mod(ij,k,4)) / max(1.e-3,cx_mod(ij,k,4))
                   dn = dnfac(lhcat) * embx ** pwmasi(lhcat)
                   it = nint(dn * 1.e6)

                   xlim = gam(it,3) * dps2 * (dps / dn) ** (gnu(4) - 1.)  /  &
                        (gamn1(4) * pwmas(lhcat) * dn ** 2)

                   dvap = max(-rx_mod(ij,k,4), vap_mod(ij,k,4) * xlim)
                   dqr = dvap * qx_mod(ij,k,4)
                   dnum = dvap * max(dpsmi(lhcat),1./embx)
                endif

                rx_mod(ij,k,3) = rx_mod(ij,k,3) - dvap
                cx_mod(ij,k,3) = cx_mod(ij,k,3) - dnum
                qr_mod(ij,k,3) = qr_mod(ij,k,3) - dqr
                rx_mod(ij,k,4) = rx_mod(ij,k,4) + dvap
                cx_mod(ij,k,4) = cx_mod(ij,k,4) + dnum
                qr_mod(ij,k,4) = qr_mod(ij,k,4) + dqr
             endif
          endif
       enddo
    enddo
    return
  end subroutine psxfer_opt

  !***************************************************************************

  subroutine getict_opt(lcat)

    use micphys, only: &
         dict,         & !INTENT(IN)
         emb0log,      & !INTENT(IN)
         rictmin,      & !INTENT(IN)
         rictmax         !INTENT(IN)

    implicit none

    ! Arguments:
    integer, intent(in) :: lcat

    ! Local Variables:
    integer :: k, ij
    real :: rict,rictmm

    do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))            !k1,k2

       do ij =1, ij_final

          if ((k>=k1_mod(ij,lcat)) .and. (k<=k2_mod(ij,lcat))) then   

             if (rx_mod(ij,k,lcat) >= 1.e-9) then

                rict = dict(lcat) * (log(emb_mod(ij,k,lcat)) - emb0log(lcat)) + 1.
                rictmm = max(rictmin,min(rictmax,rict))

                ict1_mod(ij,k,lcat) = int(rictmm)
                ict2_mod(ij,k,lcat) = ict1_mod(ij,k,lcat) + 1
                wct2_mod(ij,k,lcat) = rictmm - float(ict1_mod(ij,k,lcat))
                wct1_mod(ij,k,lcat) = 1.0 - wct2_mod(ij,k,lcat)

             endif

          endif

       enddo

    enddo
    return
  end subroutine getict_opt

  !***************************************************************************

  subroutine newtemp_opt(ngr)

    use rconstants, only: cp !INTENT(IN)

    implicit none

    ! Arguments:
    integer, intent(in)              :: ngr

    ! Local Variables:
    integer :: k, ij
    real    :: rslf,rsif

    do k = minval(k1_mod(:,10)), maxval(k2_mod(:,10))            !kf1,kf2

       do ij =1, ij_final

          if ((k>=k1_mod(ij,10)) .and. (k<=k2_mod(ij,10))) then

             tairc_mod(ij,k) = tairstrc_mod(ij,k) + sa_mod(ij,k,1) * &
                  (rvstr_mod(ij,k) - rv_opt_tile(ij,k))
             tair_mod(ij,k)  = tairc_mod(ij,k) + 273.16
             theta_opt_tile(ij,k) = tair_mod(ij,k) * cp / pitot_mod(ij,k)

             rvlsair_mod(ij,k) = rslf(press_mod(ij,k), tair_mod(ij,k))
             rvisair_mod(ij,k) = rsif(press_mod(ij,k), tair_mod(ij,k))
          endif
       enddo
    enddo

    return
  end subroutine newtemp_opt

  !***************************************************************************

  subroutine auto_accret_opt(ngr, dtlt)

    use micphys, only: &
         cfmas,        & !INTENT(IN)
         pwmas,        & !INTENT(IN)
         d1min,        & !INTENT(IN)
         d1max,        & !INTENT(IN)
         r2min,        & !INTENT(IN)
         r2max,        & !INTENT(IN)
         d2min,        & !INTENT(IN)
         d2max,        & !INTENT(IN)
         nd1cc,        & !INTENT(IN)
         d1ecr,        & !INTENT(IN)
         r2ecr,        & !INTENT(IN)
         nd2cr,        & !INTENT(IN)
         r2err,        & !INTENT(IN)
         nd2rr,        & !INTENT(IN)
         r1tabcc,      & !INTENT(IN)
         c1tabcc,      & !INTENT(IN)
         c2tabcc,      & !INTENT(IN)
         r1tabcr,      & !INTENT(IN)
         c1tabcr,      & !INTENT(IN)
         c2tabrr         !INTENT(IN)

    implicit none

    ! Arguments:
    integer, intent(in)             :: ngr
    real, intent(in)                :: dtlt

    ! Local Variables:
    integer :: k,id1cc,id1cr,id1crn,ir2cr,id2cr,ir2rr,id2rr,ij
    real :: dtlt3,dtlt6,dmb1cgs,dmb2cgs,r2cgs,en1cgs,ad1,ar2,d2minx,ad2  &
         ,bd1,br2,bd2,d2e,bd1cc,bd1cr,br2cr,bd2cr,br2rr,bd2rr,wd1cr  &
         ,wr2rr,wd2rr,tm1cc,tn1cc,tn2cc,tm1cr,tn1cr,tn2rr,en1cgs_2  &
         ,um1cc,un1cc,un2cc,um1cr,un1cr,un2rr,um2,un1,cfmasi1,cfmasi2  &
         ,pwmasi1,pwmasi2,wr2cr

    if (old_dtlt /= dtlt) then

       old_dtlt = dtlt

       dtlt3 = 1.e3 * dtlt
       dtlt6 = 1.e6 * dtlt
       cfmasi1 = 1. / cfmas(1)
       cfmasi2 = 1. / cfmas(2)
       pwmasi1 = 1. / pwmas(1)
       pwmasi2 = 1. / pwmas(2)

    endif

    do k = minval(k1_mod(:,1)), maxval(k2_mod(:,1))            !k1,k2

       do ij =1, ij_final

          if ((k>=k1_mod(ij,1)) .and. (k<=k2_mod(ij,1))) then

             if(rx_mod(ij,k,1) >= 1.e-9) then

                ! This subroutine works in cgs units, so convert inputs
                ! from mks

                dmb1cgs = 100. * (emb_mod(ij,k,1) * cfmasi1) ** pwmasi1
                dmb2cgs = 100. * (emb_mod(ij,k,2) * cfmasi2) ** pwmasi2
                r2cgs = 1.e-3 * rx_mod(ij,k,2) * dn0_opt_tile(ij,k)
                en1cgs = 1.e-6 * cx_mod(ij,k,1) * dn0_opt_tile(ij,k)

                ad1 = max(d1min,min(d1max,dmb1cgs))
                ar2 = max(r2min,min(r2max,r2cgs))
                d2minx = max(d2min,(r2cgs / (.1 * .5236)) ** pwmasi2)
                ad2 = max(d2minx,min(d2max,dmb2cgs))

                bd1 = alog10(ad1/d1min)
                br2 = alog10(ar2/r2min)
                bd2 = alog10(ad2/d2minx)
                d2e =  alog10(d2max / d2minx)

                bd1cc = float(nd1cc-1) * (ad1 - d1min) / (d1max - d1min) + 1.
                bd1cr = bd1 / d1ecr + 1.
                br2cr = br2 / r2ecr + 1.
                bd2cr = bd2 / d2e * float(nd2cr-1) + 1.
                br2rr = br2 / r2err + 1.
                bd2rr = bd2 / d2e * float(nd2rr-1) + 1.

                !         id1cc  =  int(bd1cc)
                id1cc  =  nint(bd1cc)
                id1cr  =  int(bd1cr)
                id1crn = nint(bd1cr)
                ir2cr  =  int(br2cr)
                id2cr  = nint(bd2cr)
                ir2rr  =  int(br2rr)
                id2rr  =  int(bd2rr)

                wd1cr = bd1cr - float(id1cr)
                wr2cr = br2cr - float(ir2cr)
                wr2rr = br2rr - float(ir2rr)
                wd2rr = bd2rr - float(id2rr)

                tm1cc = r1tabcc(id1cc)

                tn1cc = c1tabcc(id1cc)

                tn2cc = c2tabcc(id1cc)

                tm1cr = (1.-wd1cr) * ((1.-wr2cr) * r1tabcr(id1cr  ,ir2cr  ,id2cr) + &
                     wr2cr * r1tabcr(id1cr  ,ir2cr+1,id2cr)) + &
                     wd1cr * ((1.-wr2cr) * r1tabcr(id1cr+1,ir2cr  ,id2cr) + &
                     wr2cr * r1tabcr(id1cr+1,ir2cr+1,id2cr))

                tn1cr = (1.-wr2cr) * c1tabcr(id1crn,ir2cr  ,id2cr) + &
                     wr2cr * c1tabcr(id1crn,ir2cr+1,id2cr)

                tn2rr = (1.-wd2rr) * ((1.-wr2rr) * c2tabrr(ir2rr  ,id2rr  ) +&
                     wr2rr * c2tabrr(ir2rr+1,id2rr  )) + &
                     wd2rr * ((1.-wr2rr) * c2tabrr(ir2rr  ,id2rr+1) + &
                     wr2rr * c2tabrr(ir2rr+1,id2rr+1))

                en1cgs_2 = en1cgs ** 2

                um1cc = tm1cc * en1cgs_2 * dtlt3
                un1cc = tn1cc * en1cgs_2 * dtlt6
                un2cc = tn2cc * en1cgs_2 * dtlt6
                um1cr = 10. ** tm1cr * en1cgs * dtlt3
                un1cr = 10. ** tn1cr * en1cgs * dtlt6
                un2rr = 10. ** tn2rr * dtlt6

                ! The above values are amounts in kg/m^3 or #/m^3 converted
                ! in the present timestep, but must still be corrected for
                ! the effect of density on fall velocity.  Thus, they must be
                ! multiplied by (dn0i ** .5) which fall velocity is
                ! proportional to.  Also, since rxfer and enxfer are in units
                ! of kg/kg and #/kg, respectively, the above transfer amounts
                ! must also be multiplied by dn0i.  Together, these factors 
                ! make (dn0i ** 1.5).

                um2 = min(rx_mod(ij,k,1), (um1cc + um1cr)*dn0i_mod(ij,k))
                un1 = min(cx_mod(ij,k,1)*dn0_opt_tile(ij,k), (un1cc + un1cr))

                rxfer_mod(ij,k,1,2)  = rxfer_mod(ij,k,1,2) + um2
                qrxfer_mod(ij,k,1,2) = qrxfer_mod(ij,k,1,2) + um2 * qx_mod(ij,k,1)
                enxfer_mod(ij,k,1,1) = enxfer_mod(ij,k,1,1) + un1 - un2cc
                enxfer_mod(ij,k,1,2) = enxfer_mod(ij,k,1,2) + un2cc

                ! no collis breakup yet - do not use next line but use 
                ! col(2,2) in 3d micro

                !cc         enxfer(k,2,2) = enxfer(k,2,2) + un2rr

                ! aerosol loss here?

             endif
          endif
       enddo
    enddo
    return
  end subroutine auto_accret_opt

  !***************************************************************************

  subroutine effxy_opt(m1)

    use micphys, only: &
         jnmb            !INTENT(IN)

    implicit none

    ! Arguments:
    integer, intent(in)                :: m1

    ! Local Variables:
    integer :: k, ij
    !data ncall7/0/ ! use global variable in module

    !     1 = rp,rs,ra,rg,rh

    if (ncall7==0 .and. jnmb(2)>=1 .and. jnmb(3)>=1) then
       ncall7 = 7
       eff_mod(:,2:m1-1,1) = 1.0
    endif

    !     2 = cs,ca

    if (jnmb(2)>=1 .or. jnmb(3)>=1) then

       do k = minval(k1_mod(:,1)), maxval(k2_mod(:,1))           !k1(1),k2(1)

          do ij =1, ij_final

             if ((k>=k1_mod(ij,1)) .and. (k<=k2_mod(ij,1))) then

                ! Rough fit from Pruppacher and Klett Fig. 14-14 p. 496:
                ! close to curve for 404 microns.
                ! Replace with auto_accret eventually.

                if (emb_mod(ij,k,1)>9.e-13) then
                   eff_mod(ij,k,2) = min(1., 30.*(emb_mod(ij,k,1) - 9.e-13) ** .15)
                else
                   eff_mod(ij,k,2) = 0.
                endif

             endif

          enddo

       enddo
    endif

    !     3 = rr

    if (jnmb(2)>=1) then

       do k = minval(k1_mod(:,2)), maxval(k2_mod(:,2))           !k1(2),k2(2)

          do ij =1, ij_final

             if ((k>=k1_mod(ij,2)) .and. (k<=k2_mod(ij,2))) then

                if (rx_mod(ij,k,2) >= 1.e-9) then

                   ! rain breakup (old)

                   !            dmr = dn(k,2) * gnu2
                   !            if (dmr .lt. .0006) then
                   !               eff(k,3) = 1.0
                   !            elseif (dmr .gt. .001446) then
                   !               eff(k,3) = -5.0
                   !            else
                   !               eff(k,3) = exp(2300. * (dmr - .0006))
                   !            endif

                   ! rain breakup (new - temporary; 
                   ! eventually combine with autoconv/accret

                   if (emb_mod(ij,k,2)<.113e-6) then
                      eff_mod(ij,k,3) = 1.0
                   elseif (emb_mod(ij,k,2)>.158e-5) then
                      eff_mod(ij,k,3) = -5.0
                   else
                      eff_mod(ij,k,3) = 2. - exp(.1326e7*(emb_mod(ij,k,2) - .113e-6))
                   endif

                endif

             endif

          enddo

       enddo

    endif

    !     4 = pp,ps,pa

    if (jnmb(5)>=1) then
       do k = minval(k1_mod(:,3)), maxval(k2_mod(:,3))            !k1(3),k2(3)

          do ij =1, ij_final

             if ((k>=k1_mod(ij,3)) .and. (k<=k2_mod(ij,3))) then

                if (abs(tx_mod(ij,k,3)+14.)<=2.) then
                   eff_mod(ij,k,4) = 1.4
                else
                   eff_mod(ij,k,4) = min(0.2, 10.**(0.035 * tx_mod(ij,k,3) - 0.7))
                endif

             endif

          enddo

       enddo

       !     5 = ss,sa

       do k = minval(k1_mod(:,4)), maxval(k2_mod(:,4))            !k1(4),k2(4)

          do ij =1, ij_final

             if ((k>=k1_mod(ij,4)) .and. (k<=k2_mod(ij,4))) then

                if (abs(tx_mod(ij,k,4)+14.)<=2.) then
                   eff_mod(ij,k,5) = 1.4
                else
                   eff_mod(ij,k,5) = min(0.2, 10.**(0.035 * tx_mod(ij,k,4) - 0.7))
                endif

             endif

          enddo

       enddo

       !     6 = aa

       do k = minval(k1_mod(:,5)), maxval(k2_mod(:,5))            !k1(5),k2(5)

          do ij =1, ij_final

             if ((k>=k1_mod(ij,5)) .and. (k<=k2_mod(ij,5))) then

                if (rx_mod(ij,k,5)>=1.e-9) then

                   if (abs(tx_mod(ij,k,5)+14.)<=2.) then
                      eff_mod(ij,k,6) = 1.4
                   elseif (tx_mod(ij,k,5)>=-1.) then
                      eff_mod(ij,k,6) = 1.
                   else
                      eff_mod(ij,k,6) = min(0.2, 10.**(0.035 * tx_mod(ij,k,5) - 0.7))
                   endif

                endif

             endif

          enddo

       enddo
    endif

    !     7 = pg,sg,ag,gg,gh

    if (jnmb(6)>=1) then
       do k = minval(k1_mod(:,6)), maxval(k2_mod(:,6))            !k1(6),k2(6)

          do ij =1, ij_final

             if ((k>=k1_mod(ij,6)) .and. (k<=k2_mod(ij,6))) then

                if (qr_mod(ij,k,6)>0.) then
                   eff_mod(ij,k,7) = 1.0
                else
                   eff_mod(ij,k,7) = min(0.2, 10.**(0.035 * tx_mod(ij,k,6) - 0.7))
                endif

             endif

          enddo

       enddo
    endif

    !     8 = ph,sh,ah,gh

    if (jnmb(7)>=1) then
       do k = minval(k1_mod(:,7)), maxval(k2_mod(:,7))            !k1(7),k2(7)

          do ij =1, ij_final

             if ((k>=k1_mod(ij,7)) .and. (k<=k2_mod(ij,7))) then

                if (rx_mod(ij,k,7)>=1.e-9) then

                   if (qr_mod(ij,k,7)>0.) then
                      eff_mod(ij,k,8) = 1.0
                   else
                      eff_mod(ij,k,8) = min(0.2, 10.**(0.035 * tx_mod(ij,k,7) - 0.7))
                   endif

                endif

             endif

          enddo

       enddo
    endif

    !     9 = cg,ch

    if (jnmb(2)>=1 .or. jnmb(3)>=1) then

       do k = minval(k1_mod(:,1)), maxval(k2_mod(:,1))            !k1(1),k2(1)

          do ij =1, ij_final

             if ((k>=k1_mod(ij,1)) .and. (k<=k2_mod(ij,1))) then


                ! Rough fit from Pruppacher and Klett Fig. 14-11 p. 485:
                !  close to curves for 142 and 305 microns.  Replace with
                !  auto_accret eventually.

                if (emb_mod(ij,k,1)>3.4e-14) then
                   eff_mod(ij,k,9) = min(1., 1426.*(emb_mod(ij,k,1) - 3.4e-14)**.28)
                else
                   eff_mod(ij,k,9) = 0.
                endif

             endif

          enddo

       enddo
    endif

    !     10 = hh (trial)

    if (jnmb(7)>=1) then
       do k = minval(k1_mod(:,7)), maxval(k2_mod(:,7))            !k1(7),k2(7)

          do ij =1, ij_final

             if ((k>=k1_mod(ij,7)) .and. (k<=k2_mod(ij,7))) then

                eff_mod(ij,k,10) = max(0.,.1 + .005 * tx_mod(ij,k,7))

             endif

          enddo

       enddo
    endif

    return
  end subroutine effxy_opt

  !***************************************************************************

  subroutine cols_opt(mx, mc1, lcat)

    use micphys, only: &
         ipairc,       & !INTENT(IN)
         coltabc         !INTENT(IN)

    implicit none

    ! Arguments:
    integer, intent(in) :: mx,mc1,lcat

    ! Local Variables:
    integer :: ipc,k, ij
    real :: colnum,tabval

    do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))            !k1,k2

       do ij =1, ij_final

          if ((k>=k1_mod(ij,lcat)) .and. (k<=k2_mod(ij,lcat))) then

             if(rx_mod(ij,k,mx)>=1.e-9) then
                ipc = ipairc(jhcat_mod(ij,k,mx), jhcat_mod(ij,k,mx))

                tabval = wct1_mod(ij,k,mx)**2 * &
                     coltabc(ict1_mod(ij,k,mx), ict1_mod(ij,k,mx),ipc) + &
                     2.*wct1_mod(ij,k,mx)*wct2_mod(ij,k,mx) * &
                     coltabc(ict1_mod(ij,k,mx), ict2_mod(ij,k,mx),ipc) + &
                     wct2_mod(ij,k,mx)**2 * &
                     coltabc(ict2_mod(ij,k,mx), ict2_mod(ij,k,mx),ipc)

                colnum = colfacc_mod(ij,k) * eff_mod(ij,k,mc1) * &
                     cx_mod(ij,k,mx)**2 * 10.**(-tabval)
                enxfer_mod(ij,k,mx,mx) = enxfer_mod(ij,k,mx,mx) + &
                     min(0.5 * cx_mod(ij,k,mx), colnum)

             endif

          endif

       enddo

    enddo
    return
  end subroutine cols_opt

  !***************************************************************************

  subroutine col3344_opt(mx, mz, mc1, lcat)

    use micphys, only: &
         ipairr,       & !INTENT(IN)
         ipairc,       & !INTENT(IN)
         coltabr,      & !INTENT(IN)
         jnmb,         & !INTENT(IN)
         coltabc         !INTENT(IN)

    implicit none

    ! Arguments:
    integer, intent(in) :: mx,mz,mc1,lcat

    ! Local Variables:
    integer :: k,ip,ipc, ij
    real :: c1,tabvalx,colamt,tabvaln,colnum

    do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))            !k1,k2

       do ij =1, ij_final

          if ((k>=k1_mod(ij,lcat)) .and. (k<=k2_mod(ij,lcat))) then

             if(rx_mod(ij,k,mx)>=1.e-9) then

                ip = ipairr(jhcat_mod(ij,k,mx), jhcat_mod(ij,k,mx))
                ipc = ipairc(jhcat_mod(ij,k,mx), jhcat_mod(ij,k,mx))
                c1 = eff_mod(ij,k,mc1) * cx_mod(ij,k,mx)**2

                tabvalx = wct1_mod(ij,k,mx)**2 * &
                     coltabr(ict1_mod(ij,k,mx), ict1_mod(ij,k,mx),ip) + &
                     2. * wct1_mod(ij,k,mx) * wct2_mod(ij,k,mx) * &
                     coltabr(ict1_mod(ij,k,mx), ict2_mod(ij,k,mx),ip) + &
                     wct2_mod(ij,k,mx)**2 * &
                     coltabr(ict2_mod(ij,k,mx), ict2_mod(ij,k,mx),ip)

                colamt = min(rx_mod(ij,k,mx), colfacr2_mod(ij,k)*c1*10.**(-tabvalx))
                rxfer_mod(ij,k,mx,mz) = rxfer_mod(ij,k,mx,mz) + colamt
                qrxfer_mod(ij,k,mx,mz) = qrxfer_mod(ij,k,mx,mz) + colamt*qx_mod(ij,k,mx)

                if (jnmb(mz)>=5) then

                   tabvaln = wct1_mod(ij,k,mx)**2 * &
                        coltabc(ict1_mod(ij,k,mx), ict1_mod(ij,k,mx),ipc) + &
                        2. * wct1_mod(ij,k,mx) * wct2_mod(ij,k,mx) * &
                        coltabc(ict1_mod(ij,k,mx), ict2_mod(ij,k,mx),ipc) + &
                        wct2_mod(ij,k,mx)**2 * &
                        coltabc(ict2_mod(ij,k,mx), ict2_mod(ij,k,mx),ipc)

                   colnum = min(0.5*cx_mod(ij,k,mx), &
                        colfacc2_mod(ij,k)*c1*10.**(-tabvaln))
                   enxfer_mod(ij,k,mx,mz) = enxfer_mod(ij,k,mx,mz) + colnum
                   enxfer_mod(ij,k,mx,mx) = enxfer_mod(ij,k,mx,mx) + colnum

                endif
             endif

          endif

       enddo

    enddo
    return
  end subroutine col3344_opt

  !***************************************************************************

  subroutine col3443_opt(mx, my, mz)

    use micphys, only: &
         ipairr,       & !INTENT(IN)
         ipairc,       & !INTENT(IN)
         coltabr,      & !INTENT(IN)
         coltabc         !INTENT(IN)

    implicit none

    ! Arguments:
    integer, intent(in) :: mx,my,mz

    ! Local Variables:
    integer :: k,jhcatx,jhcaty,ipxy,ipyx,ipc, ij
    real :: c1,tabvalx,rcx,tabvaly,rcy,tabvaln,colnum

    do k = max(minval(k1_mod(:,3)), minval(k1_mod(:,4))), &
         min(maxval(k2_mod(:,3)), maxval(k2_mod(:,4)))            !k1,k2

       do ij =1, ij_final

          if ((k>=max(k1_mod(ij,3),k1_mod(ij,4))) .and. &
               (k<=min(k2_mod(ij,3),k2_mod(ij,4)))) then

             if(rx_mod(ij,k,mx)>=1.e-9 .and. rx_mod(ij,k,my)>=1.e-9) then

                jhcatx = jhcat_mod(ij,k,mx)
                jhcaty = jhcat_mod(ij,k,my)
                ipxy = ipairr(jhcatx,jhcaty)
                ipyx = ipairr(jhcaty,jhcatx)
                ipc  = ipairc(jhcatx,jhcaty)
                c1 = eff_mod(ij,k,4) * cx_mod(ij,k,mx) * cx_mod(ij,k,my)

                tabvalx = wct1_mod(ij,k,mx) * wct1_mod(ij,k,my) * &
                     coltabr(ict1_mod(ij,k,mx), ict1_mod(ij,k,my),ipxy) + &
                     wct2_mod(ij,k,mx) * wct1_mod(ij,k,my) * &
                     coltabr(ict2_mod(ij,k,mx), ict1_mod(ij,k,my),ipxy) + &
                     wct1_mod(ij,k,mx) * wct2_mod(ij,k,my) * &
                     coltabr(ict1_mod(ij,k,mx), ict2_mod(ij,k,my),ipxy) + &
                     wct2_mod(ij,k,mx) * wct2_mod(ij,k,my) * &
                     coltabr(ict2_mod(ij,k,mx), ict2_mod(ij,k,my),ipxy)
                rcx = min(rx_mod(ij,k,mx), c1*colfacr_mod(ij,k)*10.**(-tabvalx))

                tabvaly = wct1_mod(ij,k,my) * wct1_mod(ij,k,mx) * &
                     coltabr(ict1_mod(ij,k,my), ict1_mod(ij,k,mx),ipyx) + &
                     wct2_mod(ij,k,my) * wct1_mod(ij,k,mx) * &
                     coltabr(ict2_mod(ij,k,my), ict1_mod(ij,k,mx),ipyx) + &
                     wct1_mod(ij,k,my) * wct2_mod(ij,k,mx) * &
                     coltabr(ict1_mod(ij,k,my), ict2_mod(ij,k,mx),ipyx) + &
                     wct2_mod(ij,k,my) * wct2_mod(ij,k,mx) * &
                     coltabr(ict2_mod(ij,k,my), ict2_mod(ij,k,mx),ipyx)
                rcy = min(rx_mod(ij,k,my), c1*colfacr_mod(ij,k)*10.**(-tabvaly))

                rxfer_mod(ij,k,mx,mz) = rxfer_mod(ij,k,mx,mz) + rcx
                qrxfer_mod(ij,k,mx,mz) = qrxfer_mod(ij,k,mx,mz) + rcx * qx_mod(ij,k,mx)

                rxfer_mod(ij,k,my,mz) = rxfer_mod(ij,k,my,mz) + rcy
                qrxfer_mod(ij,k,my,mz) = qrxfer_mod(ij,k,my,mz) + rcy * qx_mod(ij,k,my)

                tabvaln = wct1_mod(ij,k,mx) * wct1_mod(ij,k,my) * &
                     coltabc(ict1_mod(ij,k,mx), ict1_mod(ij,k,my),ipc) + &
                     wct2_mod(ij,k,mx) * wct1_mod(ij,k,my) * &
                     coltabc(ict2_mod(ij,k,mx), ict1_mod(ij,k,my),ipc) + &
                     wct1_mod(ij,k,mx) * wct2_mod(ij,k,my) * &
                     coltabc(ict1_mod(ij,k,mx), ict2_mod(ij,k,my),ipc) + &
                     wct2_mod(ij,k,mx) * wct2_mod(ij,k,my) * &
                     coltabc(ict2_mod(ij,k,mx), ict2_mod(ij,k,my),ipc)
                colnum = c1*colfacc_mod(ij,k)*10.**(-tabvaln)

                if (cx_mod(ij,k,mx) > cx_mod(ij,k,my)) then
                   enxfer_mod(ij,k,my,mz) = min(cx_mod(ij,k,my),colnum)
                   enxfer_mod(ij,k,mx,mx) = min(cx_mod(ij,k,mx),colnum)
                else
                   enxfer_mod(ij,k,mx,mz) = min(cx_mod(ij,k,mx),colnum)
                   enxfer_mod(ij,k,my,my) = min(cx_mod(ij,k,my),colnum)
                endif

                ! also loss for aerosol

             endif

          endif

       enddo

    enddo
    return
  end subroutine col3443_opt

  !***************************************************************************

  subroutine col1_opt(mx, my, mz, mc4, lcat1, lcat2)

    use micphys, only: &
         ipairr,       & !INTENT(IN)
         ipairc,       & !INTENT(IN)
         coltabr,      & !INTENT(IN)
         jnmb,         & !INTENT(IN)
         coltabc         !INTENT(IN)

    implicit none

    ! Arguments:
    integer, intent(in) :: mx,my,mz,mc4,lcat1,lcat2

    ! Local Variables:
    integer :: k,ipxy,ipc, ij
    real :: c1,tabvalx,rcx,tabvaln,colnum

    do k = max(minval(k1_mod(:,lcat1)), minval(k1_mod(:,lcat2))), &
         min(maxval(k2_mod(:,lcat1)), maxval(k2_mod(:,lcat2))) !k1,k2

       do ij =1, ij_final

          if ((k>=max(k1_mod(ij,lcat1), k1_mod(ij,lcat2))) .and. &
               (k<=min(k2_mod(ij,lcat1), k2_mod(ij,lcat2)))) then

             if(rx_mod(ij,k,mx)>=1.e-9 .and. rx_mod(ij,k,my)>=1.e-9) then
                ipxy = ipairr(jhcat_mod(ij,k,mx), jhcat_mod(ij,k,my))
                ipc  = ipairc(jhcat_mod(ij,k,mx), jhcat_mod(ij,k,my))
                c1 = eff_mod(ij,k,mc4) * cx_mod(ij,k,mx) * cx_mod(ij,k,my)

                tabvalx = wct1_mod(ij,k,mx) * wct1_mod(ij,k,my) * &
                     coltabr(ict1_mod(ij,k,mx), ict1_mod(ij,k,my), ipxy) + &
                     wct2_mod(ij,k,mx) * wct1_mod(ij,k,my) * &
                     coltabr(ict2_mod(ij,k,mx), ict1_mod(ij,k,my), ipxy) + &
                     wct1_mod(ij,k,mx) * wct2_mod(ij,k,my) * &
                     coltabr(ict1_mod(ij,k,mx), ict2_mod(ij,k,my), ipxy) + &
                     wct2_mod(ij,k,mx) * wct2_mod(ij,k,my) * &
                     coltabr(ict2_mod(ij,k,mx), ict2_mod(ij,k,my),ipxy)

                rcx = min(rx_mod(ij,k,mx), c1*colfacr_mod(ij,k)*10.**(-tabvalx))
                rxfer_mod(ij,k,mx,mz) = rxfer_mod(ij,k,mx,mz) + rcx
                qrxfer_mod(ij,k,mx,mz) = qrxfer_mod(ij,k,mx,mz) + rcx*qx_mod(ij,k,mx)

                if (jnmb(mx)>=5) then
                   tabvaln = wct1_mod(ij,k,mx) * wct1_mod(ij,k,my) * &
                        coltabc(ict1_mod(ij,k,mx), ict1_mod(ij,k,my),ipc) + &
                        wct2_mod(ij,k,mx) * wct1_mod(ij,k,my) * &
                        coltabc(ict2_mod(ij,k,mx), ict1_mod(ij,k,my),ipc) + &
                        wct1_mod(ij,k,mx) * wct2_mod(ij,k,my) * &
                        coltabc(ict1_mod(ij,k,mx), ict2_mod(ij,k,my),ipc) + &
                        wct2_mod(ij,k,mx) * wct2_mod(ij,k,my) * &
                        coltabc(ict2_mod(ij,k,mx), ict2_mod(ij,k,my),ipc)

                   colnum = c1 * colfacc_mod(ij,k) * 10.**(-tabvaln)
                   enxfer_mod(ij,k,mx,mx) = enxfer_mod(ij,k,mx,mx) + &
                        min(colnum,cx_mod(ij,k,mx))

                   ! also loss for aerosol

                endif

             endif

          endif

       enddo

    enddo
    return
  end subroutine col1_opt

  !**************************************************************************

  subroutine col2_opt(ngr,mx,my,mz,mc2,dtlt,lcat1,lcat2)

    use micphys, only: &
         ipairr,       & !INTENT(IN)
         ipairc,       & !INTENT(IN)
         coltabr,      & !INTENT(IN)
         jnmb,         & !INTENT(IN)
         coltabc,      & !INTENT(IN)
         sipfac,       & !INTENT(IN)
         pwmasi,       & !INTENT(IN)
         emb1,         & !INTENT(IN)
         gamsip13,     & !INTENT(IN)
         gamsip24,     & !INTENT(IN)
         emb0            !INTENT(IN)

    implicit none

    ! Arguments:
    integer, intent(in)             :: ngr, mx, my, mz, mc2, lcat1, lcat2
    real, intent(in)                :: dtlt

    ! Local Variables:
    integer :: k, jhcatx, jhcaty, ipxy, ipyx, ipc, it, ij
    real :: c1(ij_final), c2(ij_final), tabvalx(ij_final), rcx(ij_final),   &
         tabvaly(ij_final), rcy(ij_final), tabvaln(ij_final),               &
         colnum0(ij_final), colnum(ij_final), rcoal(ij_final),              &
         qrcx(ij_final), qrcy(ij_final), qrcoal(ij_final), qcoal(ij_final), &
         fracliq(ij_final), tcoal(ij_final), coalliq(ij_final),             &
         area(ij_final), cn13(ij_final), cn24(ij_final), sip(ij_final),     &
         rsip(ij_final), qrsip(ij_final), rfinlz(ij_final), xtoz(ij_final)
    !coalice,

    real, dimension(15), parameter :: alpha = &
         !  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    (/00.,00.,00., 1., 1., 1., 1.,00.,00.,00.,00., 1., 1., 1., 1./)
    real, dimension(15), parameter :: beta  = &
         !  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    (/00.,00.,00.,1.5,1.1,0.0,0.0,00.,00.,00.,00.,1.2,1.1,1.1,1.3/)

    do k = max(minval(k1_mod(:,lcat1)), minval(k1_mod(:,lcat2))),&
         min(maxval(k2_mod(:,lcat1)), maxval(k2_mod(:,lcat2))) !k1,k2

       do ij =1, ij_final

          if ((k>=max(k1_mod(ij,lcat1), k1_mod(ij,lcat2))) .and. &
               (k<=min(k2_mod(ij,lcat1), k2_mod(ij,lcat2)))) then

             if(rx_mod(ij,k,mx)>=1.e-9 .and. rx_mod(ij,k,my)>=1.e-9) then

                jhcatx = jhcat_mod(ij,k,mx)
                jhcaty = jhcat_mod(ij,k,my)
                ipxy = ipairr(jhcatx,jhcaty)
                ipyx = ipairr(jhcaty,jhcatx)
                ipc  = ipairc(jhcatx,jhcaty)
                c2(ij) = cx_mod(ij,k,mx) * cx_mod(ij,k,my)
                c1(ij) = eff_mod(ij,k,mc2) * c2(ij)

                tabvalx(ij) = wct1_mod(ij,k,mx) * wct1_mod(ij,k,my) * &
                     coltabr(ict1_mod(ij,k,mx), ict1_mod(ij,k,my),ipxy) + &
                     wct2_mod(ij,k,mx) * wct1_mod(ij,k,my) * &
                     coltabr(ict2_mod(ij,k,mx), ict1_mod(ij,k,my),ipxy) + &
                     wct1_mod(ij,k,mx) * wct2_mod(ij,k,my) * &
                     coltabr(ict1_mod(ij,k,mx), ict2_mod(ij,k,my),ipxy) + &
                     wct2_mod(ij,k,mx) * wct2_mod(ij,k,my) * &
                     coltabr(ict2_mod(ij,k,mx), ict2_mod(ij,k,my),ipxy)

                rcx(ij) = min(rx_mod(ij,k,mx),c1(ij) *         &
                     colfacr_mod(ij,k) * 10. ** (-tabvalx(ij)))

                tabvaly(ij) = wct1_mod(ij,k,my) * wct1_mod(ij,k,mx) * &
                     coltabr(ict1_mod(ij,k,my), ict1_mod(ij,k,mx),ipyx) + &
                     wct2_mod(ij,k,my) * wct1_mod(ij,k,mx) * &
                     coltabr(ict2_mod(ij,k,my), ict1_mod(ij,k,mx),ipyx) + &
                     wct1_mod(ij,k,my) * wct2_mod(ij,k,mx) * &
                     coltabr(ict1_mod(ij,k,my), ict2_mod(ij,k,mx),ipyx) + &
                     wct2_mod(ij,k,my) * wct2_mod(ij,k,mx) * &
                     coltabr(ict2_mod(ij,k,my), ict2_mod(ij,k,mx),ipyx)

                rcy(ij) = min(rx_mod(ij,k,my), &
                     c1(ij)*colfacr_mod(ij,k)*10.**(-tabvaly(ij)))

                if (jnmb(mx)>=5 .or. jnmb(my)>=5) then

                   tabvaln(ij) = wct1_mod(ij,k,mx) * wct1_mod(ij,k,my) * &
                        coltabc(ict1_mod(ij,k,mx), ict1_mod(ij,k,my),ipc) + &
                        wct2_mod(ij,k,mx) * wct1_mod(ij,k,my) * &
                        coltabc(ict2_mod(ij,k,mx), ict1_mod(ij,k,my),ipc) + &
                        wct1_mod(ij,k,mx) * wct2_mod(ij,k,my) * &
                        coltabc(ict1_mod(ij,k,mx), ict2_mod(ij,k,my),ipc) + &
                        wct2_mod(ij,k,mx) * wct2_mod(ij,k,my) * &
                        coltabc(ict2_mod(ij,k,mx), ict2_mod(ij,k,my),ipc)

                   colnum0(ij) = c2(ij) * colfacc_mod(ij,k) * 10.**(-tabvaln(ij))
                   colnum(ij) = colnum0(ij) * eff_mod(ij,k,mc2)

                endif

                rcoal(ij) = rcx(ij) + rcy(ij)
                qrcx(ij) = rcx(ij) * qx_mod(ij,k,mx)
                qrcy(ij) = rcy(ij) * qx_mod(ij,k,my)
                qrcoal(ij) = qrcx(ij) + qrcy(ij)
                qcoal(ij) = qrcoal(ij) / (1.e-13 + rcoal(ij))

                call qtc(qcoal(ij),tcoal(ij),fracliq(ij))

                coalliq(ij) = rcoal(ij) * fracliq(ij)
                !coalice = rcoal(ij) - coalliq(ij)

                ! secondary ice production: cn24 is the number fraction of
                ! collected cloud droplets larger than 24 microns and is 
                ! obtained from an incomplete gamma function table.  cn13 is
                ! the fraction of collected cloud droplets smaller than 13
                ! microns.  area is cross section area of collecting ice per
                ! m^3 of atmospheric volume.

                if (tcoal(ij)>-8. .and. tcoal(ij)<-3.) then

                   area(ij) = cx_mod(ij,k,my) * dn0_opt_tile(ij,k)*sipfac(jhcaty)*  &
                        emb_mod(ij,k,my)**(2.*pwmasi(jhcaty))
                   it = nint(emb_mod(ij,k,mx) / emb1(1) * 5000.)
                   cn13(ij) = colnum(ij) * gamsip13(it) / (area(ij) * dtlt)
                   cn24(ij) = min(cx_mod(ij,k,mx)*dn0_opt_tile(ij,k),colnum0(ij)) * &
                        gamsip24(it)
                   sip(ij) = 9.1e-10 * cn24(ij) * cn13(ij) ** .93

                   if (tcoal(ij)<-5.) then
                      sip(ij) = 0.33333 * (tcoal(ij) + 8.) * sip(ij)
                   else
                      sip(ij) = -0.5 * (tcoal(ij) + 3.) * sip(ij)
                   endif

                   rsip(ij) = sip(ij) * emb0(3) * dn0i_mod(ij,k)

                   qrsip(ij) = qcoal(ij) * rsip(ij)

                   rcoal(ij) = rcoal(ij) - rsip(ij)
                   qrcoal(ij) = qrcoal(ij) - qrsip(ij)

                   enxfer_mod(ij,k,mx,3) = enxfer_mod(ij,k,mx,3) + sip(ij)

                   rxfer_mod(ij,k,mx,3) =  rxfer_mod(ij,k,mx,3) + rsip(ij)

                   qrxfer_mod(ij,k,mx,3) = qrxfer_mod(ij,k,mx,3) + qrsip(ij)

                endif

                ! ALWAYS NEED (ALPHA + BETA)>=1 but in the (rare) case that
                ! fracliq may be a little larger than fracx due to collected
                ! liquid being above 0C, need (ALPHA+BETA) to be at least 1.1
                ! or 1.2, or need ALPHA itself to be at least 1.0.

                rfinlz(ij) = min(rcoal(ij), &
                     alpha(jhcaty) * coalliq(ij) + beta(jhcaty) * rcx(ij))

                xtoz(ij) = min(rcx(ij),rfinlz(ij))

                rxfer_mod(ij,k,mx,mz) = rxfer_mod(ij,k,mx,mz) + xtoz(ij)
                rxfer_mod(ij,k,mx,my) = rxfer_mod(ij,k,mx,my) + rcx(ij) - xtoz(ij)

                if (my/=mz) &
                     rxfer_mod(ij,k,my,mz) = rxfer_mod(ij,k,my,mz)+rfinlz(ij)-xtoz(ij)

                qrxfer_mod(ij,k,mx,mz) = qrxfer_mod(ij,k,mx,mz)+qx_mod(ij,k,mx)*xtoz(ij)
                qrxfer_mod(ij,k,mx,my) = qrxfer_mod(ij,k,mx,my)+ &
                     qx_mod(ij,k,mx)*(rcx(ij) - xtoz(ij))

                if (my/=mz) &
                     qrxfer_mod(ij,k,my,mz) = qrxfer_mod(ij,k,my,mz) + &
                     qx_mod(ij,k,my) * (rfinlz(ij) - xtoz(ij))

                enxfer_mod(ij,k,mx,mx) = enxfer_mod(ij,k,mx,mx) + &
                     min(colnum(ij),cx_mod(ij,k,mx))

                if (my/=mz) &
                     enxfer_mod(ij,k,my,mz) =enxfer_mod(ij,k,my,mz) + &
                     (rfinlz(ij)-xtoz(ij))* &
                     min(colnum(ij),cx_mod(ij,k,my)) / (1.e-20 + rcy(ij))

                ! BUT NEED TO CHANGE THE ABOVE FOR 177 COLLECTION BECAUSE X=Y

                ! also include loss of aerosol

             endif

          endif

       enddo

    enddo
    return
  end subroutine col2_opt

  !**************************************************************************

  subroutine col3_opt(mx, my, mz, lcat1, lcat2)

    use micphys, only: &
         ipairr,       & !INTENT(IN)
         ipairc,       & !INTENT(IN)
         coltabr,      & !INTENT(IN)
         jnmb,         & !INTENT(IN)
         coltabc         !INTENT(IN)

    implicit none
    ! Arguments:
    integer, intent(in) :: mx, my, mz, lcat1, lcat2

    ! Local Variables:
    integer :: k, ipxy, ipyx, ipc, jhcaty, ij
    real :: c1(ij_final), tabvalx(ij_final), rcx(ij_final),                 &
         tabvaly(ij_final), rcy(ij_final), tabvaln(ij_final),               &
         colnum(ij_final), colnumx(ij_final), colnumy(ij_final),            &
         coalnum(ij_final), rcoal(ij_final), qrcx(ij_final), qrcy(ij_final),&
         qrcoal(ij_final), qcoal(ij_final), fracliq(ij_final),              &
         coalliq(ij_final), xtoz(ij_final),rfinlz(ij_final), tcoal,         &
         cfinlz(ij_final) !coalice,
    real, dimension(15),parameter :: alpha = &
         !  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    (/00.,00., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1./)
    real, dimension(15),parameter :: beta = &
         !  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    (/00.,00., 2., 2., 2., 1., 0., 2., 2., 2., 2., 2., 2., 2., 2./)

    do k = max(minval(k1_mod(:,lcat1)), minval(k1_mod(:,lcat2))), &
         min(maxval(k2_mod(:,lcat1)), maxval(k2_mod(:,lcat2)))    !k1,k2

       do ij =1, ij_final

          if ((k>=max(k1_mod(ij,lcat1), k1_mod(ij,lcat2))) .and. &
               (k<=min(k2_mod(ij,lcat1), k2_mod(ij,lcat2)))) then

             if(rx_mod(ij,k,mx)>=1.e-9 .and. rx_mod(ij,k,my) .ge. 1.e-9) then

                jhcaty = jhcat_mod(ij,k,my)
                ipxy = ipairr(jhcat_mod(ij,k,mx),jhcaty)
                ipyx = ipairr(jhcaty,jhcat_mod(ij,k,mx))
                ipc  = ipairc(jhcat_mod(ij,k,mx),jhcaty)
                c1(ij) = eff_mod(ij,k,1) * cx_mod(ij,k,mx) * cx_mod(ij,k,my)

                tabvalx(ij) = wct1_mod(ij,k,mx) * wct1_mod(ij,k,my) * &
                     coltabr(ict1_mod(ij,k,mx), ict1_mod(ij,k,my),ipxy) + &
                     wct2_mod(ij,k,mx) * wct1_mod(ij,k,my) * &
                     coltabr(ict2_mod(ij,k,mx), ict1_mod(ij,k,my),ipxy) + &
                     wct1_mod(ij,k,mx) * wct2_mod(ij,k,my) * &
                     coltabr(ict1_mod(ij,k,mx), ict2_mod(ij,k,my),ipxy) + &
                     wct2_mod(ij,k,mx) * wct2_mod(ij,k,my) * &
                     coltabr(ict2_mod(ij,k,mx), ict2_mod(ij,k,my),ipxy)

                rcx(ij) = min(rx_mod(ij,k,mx), c1(ij) * colfacr_mod(ij,k) * &
                     10. ** (-tabvalx(ij)))

                tabvaly(ij) = wct1_mod(ij,k,my) * wct1_mod(ij,k,mx) * &
                     coltabr(ict1_mod(ij,k,my), ict1_mod(ij,k,mx),ipyx) + &
                     wct2_mod(ij,k,my) * wct1_mod(ij,k,mx) * &
                     coltabr(ict2_mod(ij,k,my), ict1_mod(ij,k,mx),ipyx) + &
                     wct1_mod(ij,k,my) * wct2_mod(ij,k,mx) * &
                     coltabr(ict1_mod(ij,k,my), ict2_mod(ij,k,mx),ipyx) + &
                     wct2_mod(ij,k,my) * wct2_mod(ij,k,mx) * &
                     coltabr(ict2_mod(ij,k,my), ict2_mod(ij,k,mx),ipyx)

                rcy(ij) = min(rx_mod(ij,k,my), c1(ij) * colfacr_mod(ij,k) * &
                     10. ** (-tabvaly(ij)))

                if (jnmb(mx)>=5) then

                   tabvaln(ij) = wct1_mod(ij,k,mx) * wct1_mod(ij,k,my) * &
                        coltabc(ict1_mod(ij,k,mx), ict1_mod(ij,k,my),ipc) + &
                        wct2_mod(ij,k,mx) * wct1_mod(ij,k,my) * &
                        coltabc(ict2_mod(ij,k,mx), ict1_mod(ij,k,my),ipc) + &
                        wct1_mod(ij,k,mx) * wct2_mod(ij,k,my) * &
                        coltabc(ict1_mod(ij,k,mx), ict2_mod(ij,k,my),ipc) + &
                        wct2_mod(ij,k,mx) * wct2_mod(ij,k,my) * &
                        coltabc(ict2_mod(ij,k,mx), ict2_mod(ij,k,my),ipc)

                   colnum(ij) = c1(ij) * colfacc_mod(ij,k) * 10.**(-tabvaln(ij))
                   colnumx(ij) = min(cx_mod(ij,k,mx),colnum(ij))
                   colnumy(ij) = min(cx_mod(ij,k,my),colnum(ij))
                   coalnum(ij) = min(colnumx(ij),colnumy(ij))

                endif

                rcoal(ij) = rcx(ij) + rcy(ij)
                qrcx(ij) = rcx(ij) * qx_mod(ij,k,mx)
                qrcy(ij) = rcy(ij) * qx_mod(ij,k,my)
                qrcoal(ij) = qrcx(ij) + qrcy(ij)
                qcoal(ij) = qrcoal(ij) / (1.e-20 + rcoal(ij))

                call qtc(qcoal(ij),tcoal,fracliq(ij))

                coalliq(ij) = rcoal(ij) * fracliq(ij)
                !coalice = rcoal(ij) - coalliq(ij)

                if (fracliq(ij)>=.99) then

                   rxfer_mod(ij,k,my,mx)  = rxfer_mod(ij,k,my,mx) + rcy(ij)
                   qrxfer_mod(ij,k,my,mx) = qrxfer_mod(ij,k,my,mx) + qrcy(ij)

                   if (jnmb(mx)>=5)  enxfer_mod(ij,k,my,my) = &
                   enxfer_mod(ij,k,my,my) + colnumy(ij)

                else

                   rfinlz(ij) = min(rcoal(ij), &
                        alpha(jhcaty) * coalliq(ij) + beta(jhcaty) * rcx(ij))

                   xtoz(ij) = min(rcx(ij),rfinlz(ij))

                   rxfer_mod(ij,k,mx,mz) = rxfer_mod(ij,k,mx,mz) + xtoz(ij)
                   rxfer_mod(ij,k,mx,my) = rxfer_mod(ij,k,mx,my) + rcx(ij) - &
                        xtoz(ij)
                   if (my/=mz) rxfer_mod(ij,k,my,mz) = rxfer_mod(ij,k,my,mz) + &
                        rfinlz(ij) - xtoz(ij)

                   ! NEED TO USE QCOAL TO TRANSFER Q?

                   qrxfer_mod(ij,k,mx,mz) = qrxfer_mod(ij,k,mx,mz) + &
                        qx_mod(ij,k,mx) * xtoz(ij)
                   qrxfer_mod(ij,k,mx,my) = qrxfer_mod(ij,k,mx,my) + &
                        qx_mod(ij,k,mx) * (rcx(ij) - xtoz(ij))
                   if (my/=mz) qrxfer_mod(ij,k,my,mz) = qrxfer_mod(ij,k,my,mz) +&
                        qx_mod(ij,k,my) * (rfinlz(ij) - xtoz(ij))

                   if (jnmb(mx)>=5) then

                      if (my==mz) then

                         enxfer_mod(ij,k,mx,mx) = enxfer_mod(ij,k,mx,mx) + &
                              colnumx(ij)

                      elseif (colnumy(ij)>=colnumx(ij)) then

                         cfinlz(ij) = coalnum(ij)*rfinlz(ij)/(rcoal(ij)+1.e-20)
                         enxfer_mod(ij,k,mx,mz) = enxfer_mod(ij,k,mx,mz) + &
                              cfinlz(ij)
                         enxfer_mod(ij,k,mx,mx) = enxfer_mod(ij,k,mx,mx) + &
                              colnumx(ij) - cfinlz(ij)
                         enxfer_mod(ij,k,my,my) = enxfer_mod(ij,k,my,my) + &
                              colnumy(ij)

                      else

                         cfinlz(ij) = coalnum(ij)*rfinlz(ij)/(rcoal(ij)+1.e-20)
                         enxfer_mod(ij,k,my,mz) = enxfer_mod(ij,k,my,mz) + &
                              cfinlz(ij)
                         enxfer_mod(ij,k,mx,mx) = enxfer_mod(ij,k,mx,mx) + &
                              colnumx(ij)
                         enxfer_mod(ij,k,my,my) = enxfer_mod(ij,k,my,my) + &
                              colnumy(ij) - cfinlz(ij)

                      endif

                   endif

                endif

             endif

          endif

       enddo

    enddo

    ! also include loss of aerosol

    return

  end subroutine col3_opt

  !***************************************************************************

  subroutine colxfers_opt()

    ! rloss  => scrmic1 ! INTENT(INOUT)
    ! enloss => scrmic2 ! INTENT(INOUT)

    use micphys, only: &
         jnmb            !INTENT(IN)

    implicit none

    ! Local Variables:
    integer :: k,lcat,kd1,kd2,jcat, ij

    !  All rxfer values are nonnegative.

    do lcat = 1,7
       if (jnmb(lcat)>=1) then
          kd1 = minval(k1_mod(:,lcat))
          kd2 = maxval(k2_mod(:,lcat))

          do k = kd1,kd2
             do ij=1, ij_final
                if (k>=k1_mod(ij,lcat) .and. k<=k2_mod(ij,lcat)) then
                   scrmic1_mod(ij,k) = 0. !rloss(k)  = 0.
                   scrmic2_mod(ij,k) = 0. !enloss(k) = 0.
                endif
             enddo
          enddo

          do jcat = 1,7

             ! change this to include enxfer of the same categories
             if (jnmb(jcat)>=1) then
                if (lcat/=jcat) then
                   do k = kd1,kd2
                      do ij=1, ij_final
                         if (k>=k1_mod(ij,lcat) .and.  k<=k2_mod(ij,lcat)) then
                            scrmic1_mod(ij,k) = scrmic1_mod(ij,k) + &
                                 rxfer_mod(ij,k,lcat,jcat)
                         endif
                      enddo
                   enddo
                endif

                do k = kd1,kd2
                   do ij=1, ij_final
                      if (k>=k1_mod(ij,lcat) .and. k<=k2_mod(ij,lcat)) then
                         scrmic2_mod(ij,k) = scrmic2_mod(ij,k) + &
                              enxfer_mod(ij,k,lcat,jcat)
                      endif
                   enddo
                enddo
             endif
          enddo

          do k = kd1,kd2
             do ij=1, ij_final
                if (k>=k1_mod(ij,lcat) .and. k<=k2_mod(ij,lcat)) then
                   scrmic1_mod(ij,k) = min(1.,rx_mod(ij,k,lcat) / &
                        max(1.e-20,scrmic1_mod(ij,k)))
                   scrmic2_mod(ij,k) = min(1.,cx_mod(ij,k,lcat) / &
                        max(1.e-10,scrmic2_mod(ij,k)))
                endif
             enddo
          enddo

          do jcat = 1,7
             if (jnmb(jcat)>=1) then
                if (lcat/=jcat) then
                   do k = kd1,kd2
                      do ij=1, ij_final
                         if (k>=k1_mod(ij,lcat) .and. k<=k2_mod(ij,lcat)) then
                            rxfer_mod(ij,k,lcat,jcat) =       &
                                 rxfer_mod(ij,k,lcat,jcat) * scrmic1_mod(ij,k)
                            qrxfer_mod(ij,k,lcat,jcat) =      &
                                 qrxfer_mod(ij,k,lcat,jcat) * scrmic1_mod(ij,k)
                         endif
                      enddo
                   enddo
                endif

                do k = kd1,kd2
                   do ij=1, ij_final
                      if (k>=k1_mod(ij,lcat) .and. k<=k2_mod(ij,lcat)) then
                         enxfer_mod(ij,k,lcat,jcat) =         &
                              enxfer_mod(ij,k,lcat,jcat) * scrmic2_mod(ij,k)
                      endif
                   enddo
                enddo
             endif
          enddo
       endif
    enddo

    do lcat = 1,7

       if (jnmb(lcat)>=1) then

          kd1 = minval(k1_mod(:,lcat))
          kd2 = maxval(k2_mod(:,lcat))

          do jcat = 1,7
             if (jnmb(jcat)>=1 .and. lcat/=jcat) then
                do k = kd1,kd2
                   do ij=1, ij_final
                      if (k>=k1_mod(ij,lcat) .and. k<=k2_mod(ij,lcat)) then
                         rx_mod(ij,k,lcat) = rx_mod(ij,k,lcat) - &
                              rxfer_mod(ij,k,lcat,jcat)
                         rx_mod(ij,k,jcat) = rx_mod(ij,k,jcat) + &
                              rxfer_mod(ij,k,lcat,jcat)
                         qr_mod(ij,k,lcat) = qr_mod(ij,k,lcat) - &
                              qrxfer_mod(ij,k,lcat,jcat)
                         qr_mod(ij,k,jcat) = qr_mod(ij,k,jcat) + &
                              qrxfer_mod(ij,k,lcat,jcat)
                         cx_mod(ij,k,lcat) = cx_mod(ij,k,lcat) - &
                              enxfer_mod(ij,k,lcat,jcat)
                         cx_mod(ij,k,jcat) = cx_mod(ij,k,jcat) + &
                              enxfer_mod(ij,k,lcat,jcat)
                      endif
                   enddo
                enddo
             endif
          enddo

          if (jnmb(lcat)>=5) then
             do k = kd1,kd2
                do ij=1, ij_final
                   if (k>=k1_mod(ij,lcat) .and. k<=k2_mod(ij,lcat)) then
                      cx_mod(ij,k,lcat) = cx_mod(ij,k,lcat) - &
                           enxfer_mod(ij,k,lcat,lcat)
                   endif
                enddo
             enddo
          endif

       endif
    enddo
    return
  end subroutine colxfers_opt

  !***************************************************************************

  subroutine x02_opt(ngr, lcat) 

    use rconstants, only: alli   ! INTENT(IN)

    use micphys, only:  &
         vtfac,         & ! INTENT(IN)
         pwvtmasi,      & ! INTENT(IN)
         enmlttab,      & ! INTENT(IN)
         dnfac,         & ! INTENT(IN)
         pwmasi,        & ! INTENT(IN)
         gnu,           & ! INTENT(IN)
         shedtab          ! INTENT(IN)

    implicit none

    ! Arguments
    integer, intent(in)                :: ngr, lcat

    ! Local Variables
    integer :: k,jflag,lhcat,inc,idns, ij
    real :: rinv,closs,rxinv,rmelt,fracliq,cmelt,tcoal,ricetor6,rshed,rmltshed  &
         ,qrmltshed,shedmass,fracmloss,dn

    k1_mod(:,lcat) = k1_mod(:,10)
    k2_mod(:,lcat) = 1

    do k = minval(k1_mod(:,10)),maxval(k2_mod(:,10))

       do ij=1, ij_final

          if (k>=k1_mod(ij,10) .and. k<=k2_mod(ij,10)) then

             if (rx_mod(ij,k,lcat)>=1.e-9) k2_mod(ij,lcat) = k
             if (k2_mod(ij,lcat)==1 .and. rx_mod(ij,k,lcat)<1.e-9) &
                  k1_mod(ij,lcat) = k + 1

          endif

       enddo

    enddo

    if (lcat==2 .or. lcat>=4) then
       jflag = 1

       call enemb_opt(ngr, lcat, jflag)

       do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))

          do ij=1, ij_final

             if (k>=k1_mod(ij,lcat) .and. k<=k2_mod(ij,lcat)) then

                if (rx_mod(ij,k,lcat)>=1.e-9) then

                   lhcat = jhcat_mod(ij,k,lcat)
                   vterm_mod(ij,k,lcat) = -vtfac(lhcat) * emb_mod(ij,k,lcat) ** &
                        pwvtmasi(lhcat) * denfac_mod(ij,k)

                endif

             endif

          enddo

       enddo
    endif

    if (lcat==2) then

       do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat)) !k1(lcat),k2(lcat)

          do ij=1, ij_final

             if (k>=k1_mod(ij,lcat) .and. k<=k2_mod(ij,lcat)) then

                if (rx_mod(ij,k,lcat)>=1.e-9) then

                   rxinv = 1. / rx_mod(ij,k,lcat)
                   qx_mod(ij,k,lcat) = qr_mod(ij,k,lcat) * rxinv
                   ! limit rain to under 48C and over -80C
                   qx_mod(ij,k,lcat) = max(0., min(1.6*alli, qx_mod(ij,k,lcat)))

                endif

             endif

          enddo

       enddo

    elseif (lcat==3) then

       do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))

          do ij=1, ij_final

             if (k>=k1_mod(ij,lcat) .and. k<=k2_mod(ij,lcat)) then

                if (rx_mod(ij,k,lcat)>=1.e-9) then

                   rinv = 1. / rx_mod(ij,k,lcat)
                   qx_mod(ij,k,lcat) = qr_mod(ij,k,lcat) * rinv

                   call qtc(qx_mod(ij,k,lcat),tcoal,fracliq)

                   rmelt = rx_mod(ij,k,lcat) * fracliq
                   cmelt = cx_mod(ij,k,lcat) * fracliq

                   rx_mod(ij,k,lcat) = rx_mod(ij,k,lcat) - rmelt
                   rx_mod(ij,k,1) = rx_mod(ij,k,1) + rmelt
                   cx_mod(ij,k,lcat) = cx_mod(ij,k,lcat) - cmelt
                   cx_mod(ij,k,1) = cx_mod(ij,k,1) + cmelt

                endif

             endif

          enddo

       enddo
       !
       ! meyers - source for cloud aerosol number here?
       !
    elseif (lcat==4 .or. lcat==5) then

       do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))

          do ij=1, ij_final

             if (k>=k1_mod(ij,lcat) .and. k<=k2_mod(ij,lcat)) then

                if (rx_mod(ij,k,lcat)>=1.e-9) then

                   rinv = 1. / rx_mod(ij,k,lcat)
                   qx_mod(ij,k,lcat) = qr_mod(ij,k,lcat) * rinv
                   call qtc(qx_mod(ij,k,lcat),tcoal,fracliq)

                   if (fracliq>1.e-6) then
                      rmelt = rx_mod(ij,k,lcat) * fracliq

                      ! change this??? move to rain instead ???
                      ! look at melting decisions in col2

                      ricetor6 = min(rx_mod(ij,k,lcat)-rmelt,rmelt)
                      rx_mod(ij,k,lcat) = rx_mod(ij,k,lcat) - rmelt - ricetor6
                      rx_mod(ij,k,6)    = rx_mod(ij,k,6) + rmelt + ricetor6
                      qr_mod(ij,k,6)    = qr_mod(ij,k,6) + rmelt * alli
                      qx_mod(ij,k,lcat) = 0.

                      ! keep the above the same with ricetor6
                      ! meyers - use sa melt table here? yes
                      !
                      fracmloss = (rmelt + ricetor6) * rinv
                      closs = enmlttab(int(200. * fracmloss) + 1, &
                           jhcat_mod(ij,k,lcat)) * cx_mod(ij,k,lcat)
                      cx_mod(ij,k,lcat) = cx_mod(ij,k,lcat) - closs
                      cx_mod(ij,k,6)    = cx_mod(ij,k,6) + closs
                   endif

                endif

             endif

          enddo

       enddo

    elseif (lcat==6) then

       do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))

          do ij=1, ij_final

             if (k>=k1_mod(ij,lcat) .and. k<=k2_mod(ij,lcat)) then

                if (rx_mod(ij,k,lcat)>=1.e-9) then

                   rxinv = 1. / rx_mod(ij,k,lcat)
                   qx_mod(ij,k,lcat) = qr_mod(ij,k,lcat) * rxinv
                   call qtc(qx_mod(ij,k,lcat),tcoal,fracliq)

                   if (fracliq>0.95) then
                      rx_mod(ij,k,2) = rx_mod(ij,k,2) + rx_mod(ij,k,6)
                      qr_mod(ij,k,2) = qr_mod(ij,k,2) + rx_mod(ij,k,6) * alli
                      cx_mod(ij,k,2) = cx_mod(ij,k,2) + cx_mod(ij,k,6)
                      rx_mod(ij,k,6) = 0.
                      qr_mod(ij,k,6) = 0.
                      cx_mod(ij,k,6) = 0.
                   endif

                endif

             endif

          enddo

       enddo

    elseif (lcat==7) then

       shedmass = 5.236e-7
       do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat))

          do ij=1, ij_final

             if (k>=k1_mod(ij,lcat) .and. k<=k2_mod(ij,lcat)) then

                if (rx_mod(ij,k,lcat)>=1.e-9) then

                   rxinv = 1. / rx_mod(ij,k,lcat)
                   qx_mod(ij,k,lcat) = qr_mod(ij,k,lcat) * rxinv
                   !c          qx(k,lcat) = max(-50.,qx(k,lcat))
                   call qtc(qx_mod(ij,k,lcat),tcoal,fracliq)

                   if (fracliq>0.95) then
                      rx_mod(ij,k,2) = rx_mod(ij,k,2) + rx_mod(ij,k,7)
                      qr_mod(ij,k,2) = qr_mod(ij,k,2) + rx_mod(ij,k,7) * alli
                      cx_mod(ij,k,2) = cx_mod(ij,k,2) + cx_mod(ij,k,7)
                      rx_mod(ij,k,7) = 0.
                      qr_mod(ij,k,7) = 0.
                      cx_mod(ij,k,7) = 0.
                      !
                      !  take out following IF statement?
                      !

                   elseif (fracliq>0.3) then

                      lhcat = jhcat_mod(ij,k,lcat)
                      inc = nint(200. * fracliq) + 1
                      dn = dnfac(lhcat) * emb_mod(ij,k,lcat) ** pwmasi(lhcat)
                      idns = max(1,nint(1.e3 * dn * gnu(lcat)))
                      rshed = rx_mod(ij,k,lcat) * shedtab(inc,idns)
                      !cc            rmltshed = rx(k,lcat) * rmlttab(inc) + rshed
                      rmltshed = rshed
                      qrmltshed = rmltshed * alli

                      rx_mod(ij,k,2)    = rx_mod(ij,k,2) + rmltshed
                      qr_mod(ij,k,2)    = qr_mod(ij,k,2) + qrmltshed
                      rx_mod(ij,k,lcat) = rx_mod(ij,k,lcat) - rmltshed
                      qr_mod(ij,k,lcat) = qr_mod(ij,k,lcat) - qrmltshed
                      !               closs = cx(k,lcat) * enmlttab(inc,lhcat)
                      !               cx(k,lcat) = cx(k,lcat) - closs
                      !               cx(k,2) = cx(k,2) + closs + rshed/shedmass
                      cx_mod(ij,k,2) = cx_mod(ij,k,2) + rshed / shedmass

                   endif

                endif

             endif

          enddo

       enddo

    endif
    return
  end subroutine x02_opt

  ! *****************************************************************************

  subroutine cldnuc_opt(ngr, m1, k1cnuc, k2cnuc)

    use micphys, only: &
         jnmb,         & ! INTENT(IN)
         parm,         & ! INTENT(IN)
         emb0,         & ! INTENT(IN)
         cparm,        & ! INTENT(IN)
         icloud          ! INTENT(IN)

    implicit none

    ! Arguments:
    integer, intent(in)                       :: ngr, m1
    integer, intent(out), dimension(ij_final) :: k1cnuc, k2cnuc

    ! Local Variables:
    integer :: k, jtemp(ij_final), jw(ij_final), jconcen(ij_final), ij,        &
         stoping_flag=0
    real :: rnuc, excessrv(ij_final), rcnew,tab(ij_final),                     &
         concen_tab(ij_final), cxadd(ij_final), tairc_nuc(ij_final),           &
         w_nuc(ij_final), rjw(ij_final), wtw2(ij_final), wtw1(ij_final),       &
         concen_nuc(ij_final), rjconcen(ij_final),wtconcen2(ij_final),         &
         wtconcen1(ij_final)

    real, dimension(9,7,7), parameter :: cldnuctab = reshape ( &
         (/.307,  .520,  .753,  .919,  .990,  .990,  .990,  .990,  .990,  &
         .230,  .426,  .643,  .860,  .969,  .990,  .990,  .990,  .990,  &
         .164,  .336,  .552,  .777,  .940,  .990,  .990,  .990,  .990,  &
         .098,  .254,  .457,  .701,  .892,  .979,  .990,  .990,  .990,  &
         .045,  .145,  .336,  .614,  .822,  .957,  .990,  .990,  .990,  &
         .018,  .073,  .206,  .426,  .672,  .877,  .969,  .990,  .990,  &
         .008,  .027,  .085,  .206,  .280,  .336,  .336,  .336,  .906,  &
         !
    .230,  .426,  .643,  .860,  .969,  .990,  .990,  .990,  .990,  &
         .164,  .336,  .552,  .777,  .930,  .990,  .990,  .990,  .990,  &
         .112,  .254,  .457,  .701,  .877,  .974,  .990,  .990,  .990,  &
         .073,  .184,  .365,  .583,  .822,  .949,  .990,  .990,  .990,  &
         .038,  .112,  .254,  .489,  .727,  .906,  .982,  .990,  .990,  &
         .015,  .054,  .145,  .365,  .614,  .841,  .957,  .990,  .990,  &
         .005,  .018,  .073,  .184,  .395,  .614,  .800,  .940,  .990,  &
         !
    .164,  .336,  .552,  .800,  .949,  .990,  .990,  .990,  .990,  &
         .128,  .254,  .457,  .701,  .892,  .979,  .990,  .990,  .990,  &
         .085,  .184,  .365,  .583,  .822,  .949,  .990,  .990,  .990,  &
         .054,  .128,  .280,  .489,  .727,  .906,  .982,  .990,  .990,  &
         .027,  .085,  .206,  .395,  .643,  .841,  .963,  .990,  .990,  &
         .012,  .038,  .112,  .280,  .520,  .777,  .930,  .990,  .990,  &
         .004,  .015,  .054,  .145,  .365,  .614,  .822,  .949,  .990,  &
         !
    .145,  .280,  .489,  .727,  .919,  .990,  .990,  .990,  .990,  &
         .098,  .206,  .395,  .614,  .841,  .963,  .990,  .990,  .990,  &
         .063,  .145,  .307,  .520,  .753,  .919,  .990,  .990,  .990,  &
         .038,  .098,  .230,  .426,  .643,  .860,  .963,  .990,  .990,  &
         .022,  .063,  .164,  .336,  .552,  .777,  .930,  .990,  .990,  &
         .010,  .027,  .085,  .230,  .426,  .701,  .877,  .974,  .990,  &
         .003,  .012,  .038,  .112,  .280,  .552,  .777,  .940,  .990,  &
         !
    .112,  .230,  .457,  .701,  .892,  .982,  .990,  .990,  .990,  &
         .073,  .164,  .336,  .552,  .800,  .940,  .990,  .990,  .990,  &
         .054,  .112,  .254,  .457,  .672,  .877,  .979,  .990,  .990,  &
         .032,  .085,  .184,  .365,  .583,  .800,  .940,  .990,  .990,  &
         .018,  .045,  .128,  .254,  .457,  .701,  .892,  .979,  .990,  &
         .008,  .022,  .073,  .184,  .365,  .614,  .822,  .949,  .990,  &
         .003,  .010,  .032,  .098,  .230,  .489,  .727,  .906,  .979,  &
         !
    .098,  .206,  .395,  .643,  .860,  .974,  .990,  .990,  .990,  &
         .063,  .145,  .307,  .520,  .753,  .930,  .990,  .990,  .990,  &
         .045,  .098,  .206,  .395,  .643,  .841,  .963,  .990,  .990,  &
         .027,  .063,  .145,  .307,  .520,  .753,  .919,  .990,  .990,  &
         .015,  .038,  .098,  .230,  .426,  .643,  .841,  .963,  .990,  &
         .007,  .018,  .054,  .145,  .307,  .552,  .777,  .919,  .990,  &
         .003,  .008,  .027,  .073,  .206,  .395,  .672,  .860,  .969,  &
         !
    .098,  .206,  .365,  .614,  .841,  .969,  .990,  .990,  .990,  &
         .054,  .128,  .280,  .489,  .727,  .906,  .990,  .990,  .990,  &
         .038,  .085,  .184,  .365,  .583,  .822,  .957,  .990,  .990,  &
         .022,  .063,  .128,  .280,  .457,  .701,  .892,  .982,  .990,  &
         .012,  .038,  .085,  .184,  .365,  .583,  .822,  .949,  .990,  &
         .005,  .018,  .045,  .128,  .280,  .489,  .727,  .892,  .979,  &
         .002,  .007,  .022,  .063,  .164,  .365,  .614,  .822,  .949/),&
         (/ 9, 7, 7 /) ) 

    do ij=1, ij_final
       k1cnuc(ij) = lpw_opt(ij)
       k2cnuc(ij) = 1
    enddo

    if (jnmb(1)==1 .or. jnmb(1)==4) then

       ! cloud number specified in parm(1)   

       rnuc = parm(1) * emb0(1)

       do k = minval(lpw_opt(:)), m1-1          !lpw,m1-1

          do ij=1, ij_final

             if (k>=lpw_opt(ij)) then

                excessrv(ij) = rv_opt_tile(ij,k) - 1.0001 * rvlsair_mod(ij,k)
                rcnew = 0.

                if (excessrv(ij) > 0.) then

                   rcnew = min(rnuc,.5*excessrv(ij))
                   rx_mod(ij,k,1) = rx_mod(ij,k,1) + rcnew
                   rv_opt_tile(ij,k) = rv_opt_tile(ij,k) - rcnew
                   k2cnuc(ij) = k
                   cx_mod(ij,k,1) = min(parm(1), rx_mod(ij,k,1) / emb0(1))

                elseif (k2cnuc(ij)==1) then

                   k1cnuc(ij) = k + 1

                endif

             endif

          enddo

       enddo

    elseif (jnmb(1)>=5) then

       ! cloud number predicted from ccn field

       do k = minval(lpw_opt(:)), m1-1          !lpw,m1-1

          do ij=1, ij_final

             if (k>=lpw_opt(ij)) then           

                excessrv(ij) = rv_opt_tile(ij,k) - 1.0001 * rvlsair_mod(ij,k)

                if (excessrv(ij) > 0.) then

                   tairc_nuc(ij) = tairc_mod(ij,k)

                   if (tairc_nuc(ij) < -30.) then
                      tairc_nuc(ij) = -30.
                   elseif (tairc_nuc(ij) > 30.) then
                      tairc_nuc(ij) = 30.
                   endif

                   jtemp(ij) = nint(.1 * (tairc_nuc(ij) + 30.)) + 1

                   w_nuc(ij) = wp_opt_tile(ij,k)

                   if (w_nuc(ij) < .010001) then 
                      w_nuc(ij) = .010001
                   elseif (w_nuc(ij) > 99.99) then
                      w_nuc(ij) = 99.99
                   endif

                   rjw(ij) = 2. * log10(100. * w_nuc(ij)) + 1.
                   jw(ij) = int(rjw(ij))
                   wtw2(ij) = rjw(ij) - float(jw(ij))
                   wtw1(ij) = 1. - wtw2(ij)

                   if (jnmb(1)==5) then

                      concen_nuc(ij) = cparm
                      ! ccn concen const value specified in cparm

                   elseif (jnmb(1)==6) then

                      ! concen_nuc(ij) = prof(k)
                      ! ccn concen specified vertical profile

                   elseif (jnmb(1)==7) then

                      concen_nuc(ij) = cccnx_mod(ij,k)
                      ! ccn concen predicted in cccnp

                   else

                      !print*, 'icloud set to value greater than 7: ',icloud
                      !print*, 'stopping model '
                      !stop 'icloud'
                      stoping_flag=1

                   endif

                   if (concen_nuc(ij) < 10.001e6) then
                      concen_nuc(ij) = 10.001e6
                   elseif (concen_nuc(ij) > 9999.e6) then
                      concen_nuc(ij) = 9999.e6
                   endif

                   rjconcen(ij) = 2. * log10(1.e-7 * concen_nuc(ij)) + 1.
                   jconcen(ij) = int(rjconcen(ij))
                   wtconcen2(ij) = rjconcen(ij) - float(jconcen(ij))
                   wtconcen1(ij) = 1. - wtconcen2(ij)

                   tab(ij) = wtconcen1(ij) * (wtw1(ij) *                        &
                        cldnuctab(jw(ij),jconcen(ij),jtemp(ij)) + wtw2(ij) *    &
                        cldnuctab(jw(ij)+1,jconcen(ij),jtemp(ij))) +            &
                        wtconcen2(ij) * (wtw1(ij) *                             &
                        cldnuctab(jw(ij),jconcen(ij)+1,jtemp(ij)) +             &
                        wtw2(ij) * cldnuctab(jw(ij)+1,jconcen(ij)+1,jtemp(ij)))

                   concen_tab(ij) = concen_nuc(ij) * tab(ij)

                   ! Nucleate cloud droplets only if 
                   ! concen_tab > existing cloud concentration

                   if (concen_tab(ij) > cx_mod(ij,k,1)) then

                      cxadd(ij) = concen_tab(ij) - cx_mod(ij,k,1)

                      if (cxadd(ij) > excessrv(ij) / emb0(1))                   &
                           cxadd(ij) = excessrv(ij) / emb0(1)

                      cx_mod(ij,k,1) = cx_mod(ij,k,1) + cxadd(ij)
                      rx_mod(ij,k,1) = rx_mod(ij, k,1) + excessrv(ij)
                      k2cnuc(ij) = k

                   endif

                elseif (k2cnuc(ij)==1) then
                   k1cnuc(ij) = k + 1
                endif

             endif

          enddo

       enddo

    else
       print*, 'icloud not allowed to be 2 or 3'
       print*, 'stopping model '
       stop 'icloud'
    endif

    if (stoping_flag==1) then
       print*, 'icloud set to value greater than 7: ',icloud
       print*, 'stopping model '
       stop 'icloud'
    endif

    return
  end subroutine cldnuc_opt

  !******************************************************************************

  subroutine c03_opt(ngr, lcat)

    use micphys, only: jnmb  ! INTENT(IN)

    implicit none

    ! Arguments
    integer, intent(in)             :: ngr, lcat

    ! Local Variables
    integer :: jflag

    jflag = 1
    if (jnmb(lcat)>=3) call enemb_opt(ngr, lcat, jflag)

    return
  end subroutine c03_opt

  !******************************************************************************

  subroutine icenuc_opt(ngr,m1,k1pnuc,k2pnuc,dtlt)

    use micphys, only: &
         dnfac,        & ! INTENT(IN)
         pwmasi,       & ! INTENT(IN)
         ndnc,         & ! INTENT(IN)
         ddnc,         & ! INTENT(IN)
         dtc,          & ! INTENT(IN)
         fracc,        & ! INTENT(IN)
         drhhz,        & ! INTENT(IN)
         dthz,         & ! INTENT(IN)
         frachz,       & ! INTENT(IN)
         ipris,        & ! INTENT(IN)
         emb0            ! INTENT(IN)

    implicit none

    ! Arguments:
    integer, intent(in)                       :: ngr, m1
    integer, intent(out), dimension(ij_final) :: k1pnuc, k2pnuc
    real, intent(in)                          :: dtlt

    ! Local Variables:
    integer :: k, idnc(ij_final), itc(ij_final), irhhz(ij_final),               &
         ithz(ij_final), ij
    real :: dn1(ij_final), fraccld(ij_final), ridnc(ij_final), wdnc2(ij_final), &
         tc(ij_final), ritc(ij_final), wtc2(ij_final), pbvi, ptvi, pdvi,        &
         ptotvi(ij_final), fracifn(ij_final), cldnuc(ij_final),                 &
         cldnucr(ij_final), rhhz(ij_final), haznuc(ij_final), rirhhz(ij_final), &
         wrhhz2(ij_final), thz(ij_final), rithz(ij_final), wthz2(ij_final),     &
         frachaz(ij_final), ssi(ij_final), diagni(ij_final), vapnuc(ij_final),  &
         vapnucr(ij_final), availvap(ij_final)

    ! Define ssi0 to be maximum supersaturation with respect to ice for
    ! determining total number of IFN that can nucleate in Meyers' formula
    real, parameter :: ssi0=0.40

    !
    ! implement paul's immersion freezing of rain here.  This would
    ! replace mike's homogeneous freezing of rain which was in h03.
    !

    do k = minval(k1_mod(:,1)), maxval(k2_mod(:,1))  !kc1,kc2

       do ij=1, ij_final

          if (k>=k1_mod(ij,1) .and. k<=k2_mod(ij,1)) then

             !  Homogeneous ice nucleation of cloud droplets

             ! define dn locally from emb

             dn1(ij) = dnfac(1) * emb_mod(ij,k,1) ** pwmasi(1)

             fraccld(ij) = 0.

             if (rx_mod(ij,k,1)>1.e-10 .and. tairc_mod(ij,k)<=-30.01) then

                ridnc(ij) = max(1., min(float(ndnc-1), dn1(ij)/ddnc))
                idnc(ij) = int(ridnc(ij))
                wdnc2(ij) = ridnc(ij) - float(idnc(ij))

                tc(ij) = max(-49.99, tairc_mod(ij,k))
                ritc(ij) = (tc(ij) + 50.00)/dtc + 1.0
                itc(ij) = int(ritc(ij))
                wtc2(ij) = ritc(ij) - float(itc(ij))
                fraccld(ij) = (1.-wdnc2(ij)) * (1.-wtc2(ij)) *                 &
                     fracc(idnc(ij),itc(ij), ngr) + wdnc2(ij) * (1.-wtc2(ij)) *&
                     fracc(idnc(ij)+1,itc(ij), ngr) + (1. - wdnc2(ij)) *       &
                     wtc2(ij) * fracc(idnc(ij), itc(ij)+1, ngr) + wdnc2(ij) *  &
                     wtc2(ij) * fracc(idnc(ij)+1, itc(ij)+1, ngr)

             endif

             ! Heterogeneous contact ice nucleation of cloud droplets by diffusio
             ! phoresis, thermophoresis, and Brownian motion (transport of IN)

             call contnuc_opt(rx_mod(ij,k,1), cx_mod(ij,k,1), tx_mod(ij,k,1), &
                  vap_mod(ij,k,1), press_mod(ij,k), dynvisc_mod(ij,k), &
                  thrmcon_mod(ij,k), tair_mod(ij,k), &
                  tairc_mod(ij,k), pbvi, ptvi, pdvi, ptotvi(ij), dn1(ij), dtlt)

             ! progIFN: Scale ptotvi returned from contnuc by prognosed 
             ! IFN fraction

             !::later   ptotvi = ptotvi * fracifn

             ! MIKE ADDED THIS COMMENTED ccinp(k)=ccinp(k)-ptotvi, but
             ! probably do not want sink of ccinp here.

             cldnuc(ij) = ptotvi(ij) + &
                  max(0., fraccld(ij) * cx_mod(ij,k,1) - cx_mod(ij,k,3))
             !         cldnucr = cldnuc * emb(k,1)
             cldnucr(ij) = min(rx_mod(ij,k,1), ptotvi(ij) * emb_mod(ij,k,1) + &
                  fraccld(ij) * rx_mod(ij,k,1))

             rx_mod(ij,k,3) = rx_mod(ij,k,3) + cldnucr(ij)
             rx_mod(ij,k,1) = rx_mod(ij,k,1) - cldnucr(ij)
             cx_mod(ij,k,3) = cx_mod(ij,k,3) + cldnuc(ij)
             cx_mod(ij,k,1) = cx_mod(ij,k,1) - cldnuc(ij)

          endif

       enddo

    enddo

    ! DEMOTT'S NEW SCHEME: In 4.3 and beyond, assume that it gives #/KG

    !  Homogeneous nucleation of haze

    k1pnuc(:) = 2
    k2pnuc(:) = 1

    do k = minval(lpw_opt(:)), m1-1 !lpw,m1-1

       do ij=1, ij_final

          if (k>=lpw_opt(ij)) then

             rhhz(ij) = rv_opt_tile(ij,k) / rvlsair_mod(ij,k)
             haznuc(ij) = 0.
             if (rhhz(ij)>0.82 .and. tairc_mod(ij,k)<=-35.01) then
                rirhhz(ij) = min(0.1799, rhhz(ij)-0.82)/drhhz + 1.0
                irhhz(ij) = int(rirhhz(ij))
                wrhhz2(ij) = rirhhz(ij) - float(irhhz(ij))
                thz(ij) = max(-59.99, tairc_mod(ij,k))
                rithz(ij) = (thz(ij) + 60.00)/dthz + 1.0
                ithz(ij) = int(rithz(ij))
                wthz2(ij) = rithz(ij) - float(ithz(ij))
                frachaz(ij) = (1. - wrhhz2(ij)) * (1. - wthz2(ij)) *           &
                     frachz(irhhz(ij), ithz(ij)) + wrhhz2(ij) *                &
                     (1. - wthz2(ij)) * frachz(irhhz(ij)+1, ithz(ij)) +        &
                     (1. - wrhhz2(ij)) * wthz2(ij) *                           &
                     frachz(irhhz(ij), ithz(ij)+1) + wrhhz2(ij) * wthz2(ij) *  &
                     frachz(irhhz(ij)+1, ithz(ij)+1)
                frachaz(ij) = 1. - exp(-frachaz(ij)*dtlt)
                ! OPTION 1
                haznuc(ij) = frachaz(ij) * 300.e6
                ! OPTION 2
                !           haznuc = frachaz * caero(k)
             endif

             ! meyers -  no cld aerosol source or sink here

             !  Heterogeneous nucleation by deposition condensation freezing
             !  with deposition nuclei. In 4.3 and beyond, assume that it gives
             !  #/kg.

             ssi(ij) = min(ssi0, rv_opt_tile(ij,k) / rvisair_mod(ij,k) - 1.)

             if (ssi(ij)>0. .and. tairc_mod(ij,k)<-5.) then
                fracifn(ij) = exp(12.96 * (ssi(ij) - ssi0))
             else
                fracifn(ij) = 0.
             endif

             ! Diagnose maximum number of IFN to activate based on ipris

             if (ipris==5) then
                diagni(ij) = fracifn(ij) * 1.e5
             elseif (ipris==6) then
                diagni(ij) = fracifn(ij) * dn0_opt_tile(ij,k) ** 5.4 * 1.e5
             elseif (ipris==7) then
                diagni(ij) = fracifn(ij) * cifnx_mod(ij,k)
             endif

             ! orig Meyers formula:     +      diagni = exp(6.269 + 12.96 * ssi)

             !  Combine nucleation types, and limit amounts
             ! vapnuc is #/kg_air and vapnucr is kg/kg_air

             ! BEGIN MIKE'S SECTION FOR LIMITING NUMBER OF CRYSTALS NUCLEATED
             ! BY NUMBER OF ICE CRYSTALS PRESENT ALREADY

             vapnuc(ij) = max(0., haznuc(ij) + diagni(ij) - cx_mod(ij,k,3))
             vapnucr(ij) = vapnuc(ij) * emb0(3)
             if (vapnucr(ij)>0.) then
                availvap(ij) = .5 * (rv_opt_tile(ij,k) - rvisair_mod(ij,k))
                if (vapnucr(ij)>availvap(ij)) then
                   vapnucr(ij) = min(vapnucr(ij), max(0., availvap(ij)))
                endif
             endif
             vapnuc(ij) = vapnucr(ij)/emb0(3)

             rx_mod(ij,k,3) = rx_mod(ij,k,3) + vapnucr(ij)
             cx_mod(ij,k,3) = cx_mod(ij,k,3) + vapnuc(ij)

             if (rx_mod(ij,k,3)>1.e-9) k2pnuc(ij) = k
             if (k2pnuc(ij)==1 .and. rx_mod(ij,k,3)<1.e-9) k1pnuc(ij) = k + 1

          endif

       enddo

    enddo

    ! here mike has the habit diagnosis. option 1 is to use habit
    ! at cloud top, option 2 is to use new habit at each level.
    ! need to consider other options.  how about method of formation?
    ! my question about how much of habit is due to existing ice
    ! structure, and how much is due to current growth environment
    ! (temp and supsat). relevant supsat is wrt liquid?

    return
  end subroutine icenuc_opt

  !******************************************************************************

  subroutine contnuc_opt (rx,cx,tx,vap,press,  &
       dynvisc,thrmcon,tair,tairc,pbvi,ptvi,pdvi,ptotvi,dn1,dtlt) !i

    implicit none

    ! Arguments:
    real, intent(in)    :: rx,cx,tx,vap,press,dynvisc,thrmcon,tair,tairc,dn1,dtlt
    real, intent(out)   :: pbvi,ptvi,pdvi,ptotvi

    ! Local Variables:
    real :: ana,akn,dfar,f1,f2,ft
    real, parameter :: aka = 5.39e-3, raros = 3.e-7

    !  Heterogeneous contact ice nucleation of cloud droplets by diffusio-
    !  phoresis, thermophoresis, and Brownian motion (transport of IN)
    !
    !  ana   = # IN per kg available for contact freezing 
    !            (from Meyers et al. 1992
    !          where ana was interpreted as # per m^3)
    !  akn   = Knudsen number (Walko et al. 1995, Eq. 58)
    !          [2.28e-5 = mfp * p00 / 293.15]
    !  raros = aerosol radius = 3.e-7 m from Cotton et al. (1986)
    !  dfar  = aerosol diffusivity (Pruppacher and Klett Eq. 12-15)
    !          [7.32e-25 = Boltzmann constant / (6 pi)]
    !  f1    = "function 1" (Walko et al. 1995 Eq. 55) multiplied by delta t
    !           but now cld concen in #/kg_air so (pvbi, ptvi, pdvi) all per 
    !           kg_air
    !  f2    = "function 2" (Walko et al. 1995 Eq. 56)
    !  ft    = "function ft" (Walko et al. 1995 Eq. 57)
    !  pbvi  = Brownian motion nucleation amount this timestep [#/kg_air]
    !  ptvi  = Thermophoretic nucleation amount this timestep [#/kg_air]
    !  pdvi  = Diffusiophoretic nucleation amount this timestep [#/kg_air],
    !          reformulated to use vapor diffusion directly.  Factor of 1.2
    !          is (1+sigma_va x_a) from Pruppacher and Klett Eq. 12-102
    !          divided by .622, the molecular weight ratio between water and air.

    ptotvi = 0.

    if (tx<=-2. .and. rx>1.e-10) then

       ana    = exp(4.11 - 0.262*tx)
       akn    = 2.28e-5*tair/(press*raros)
       dfar   = 7.32e-25*tair*(1. + akn)/(raros*dynvisc)
       f1     = 6.28318*dn1*cx*ana*dtlt
       f2     = thrmcon*(tairc - tx)/press
       ft     = 0.4*(1. + 1.45*akn + 0.4*akn*exp(-1./akn)) * &
            (thrmcon+2.5*akn*aka)/((1.+3.*akn)*(2.*thrmcon+5.*aka*akn+aka))
       pbvi   = f1*dfar
       ptvi   = f1*f2*ft
       pdvi   = 1.2*ana*vap
       ptotvi = max(0., pbvi + ptvi + pdvi)

    endif
    return
  end subroutine contnuc_opt

  !******************************************************************************

  subroutine pc03_opt(ngr, lcat)

    use micphys, only: &
         jnmb,         & ! INTENT(IN)
         jhcat,        & ! INTENT(IN)
         vtfac,        & ! INTENT(IN)
         pwvtmasi        ! INTENT(IN)

    implicit none

    ! Arguments
    integer, intent(in)             :: ngr, lcat

    ! Local Variables
    integer :: k, lhcat, jflag, ij

    jflag = 1
    if (jnmb(lcat)>=3) call enemb_opt(ngr, lcat, jflag)

    do k = minval(k1_mod(:,lcat)), maxval(k2_mod(:,lcat)) !k1,k2

       do ij =1, ij_final

          if ((k>=k1_mod(ij,lcat)).and.(k<=k2_mod(ij,lcat))) then

             if (rx_mod(ij,k,lcat)>=1.e-9) then

                lhcat = jhcat_mod(ij,k,lcat)
                vterm_mod(ij,k,lcat) = -vtfac(lhcat) *         &
                     emb_mod(ij,k,lcat) ** pwvtmasi(lhcat) * denfac_mod(ij,k)

             endif

          endif

       enddo

    enddo

    return
  end subroutine pc03_opt

  !*****************************************************************************

  subroutine sedim_opt(m1, lcat, ngr, nembfall, maxkfall, alphasfc, dtlti, &
       pcpfillc, pcpfillr, sfcpcp, dtlt, local_ch1, ij_final, pcprx_mod,  &
       k1_mod, k2_mod, scrmic1_mod, scrmic2_mod, scrmic3_mod, jhcat_mod, &
       rx_mod, cx_mod, qx_mod, emb_mod, dn0i_mod, accpx_mod, tairc_mod, &
       tair_mod, dn0_opt_tile, rtgt_opt_tile, pcpg_opt_tile, qpcpg_opt_tile, &
       dpcpg_opt_tile, rtp_opt_tile, thp_opt_tile, theta_opt_tile)

    use rconstants, only: cpi  ! INTENT(IN)

    use micphys, only: &
         nhcat,        & ! INTENT(IN)
         ch1,          & ! INTENT(INOUT)
         cfmas,        & ! INTENT(IN)
         ch3,          & ! INTENT(IN)
         ch2,          & ! INTENT(IN)
         dispemb0,     & ! INTENT(IN)
         cfvt            ! INTENT(IN) - ALF

    implicit none

    ! Arguments
    integer, intent(in) :: m1
    integer, intent(in) :: lcat
    integer, intent(in) :: ngr
    integer, intent(in) :: nembfall
    integer, intent(in) :: maxkfall
    real,    intent(in) :: alphasfc
    real,    intent(in) :: dtlti
    real,    intent(in) :: dtlt
    real,    intent(in) :: pcpfillr(m1,maxkfall,nembfall,nhcat)
    real,    intent(in) :: pcpfillc(m1,maxkfall,nembfall,nhcat)
    real,    intent(in) :: sfcpcp(maxkfall,nembfall,nhcat)
    real,    intent(in) :: local_ch1(:,:)

    ! novos argumentos:
    integer, intent(in) :: ij_final
    ! global
    real,    intent(inout) :: pcprx_mod(:,:)
    integer, intent(in   ) :: k1_mod(:,:)
    integer, intent(in   ) :: k2_mod(:,:)
    real,    intent(inout) :: scrmic1_mod(:,:)
    real,    intent(inout) :: scrmic2_mod(:,:)
    real,    intent(inout) :: scrmic3_mod(:,:)
    integer, intent(in   ) :: jhcat_mod(:,:,:)
    real,    intent(inout) :: rx_mod(:,:,:)
    real,    intent(inout) :: cx_mod(:,:,:)
    real,    intent(inout) :: qx_mod(:,:,:)
    real,    intent(in   ) :: emb_mod(:,:,:)
    real,    intent(in   ) :: dn0i_mod(:,:)
    real,    intent(inout) :: accpx_mod(:,:)
    real,    intent(inout) :: tairc_mod(:,:)
    real,    intent(in   ) :: tair_mod(:,:)


    ! from micro_g_opt(ngr)
    real,    intent(in   ) :: dn0_opt_tile(:,:)
    real,    intent(in   ) :: rtgt_opt_tile(:)
    real,    intent(inout) :: pcpg_opt_tile(:)
    real,    intent(inout) :: qpcpg_opt_tile(:)
    real,    intent(inout) :: dpcpg_opt_tile(:)
    real,    intent(inout) :: rtp_opt_tile(:,:)
    real,    intent(in   ) :: thp_opt_tile(:,:)
    real,    intent(in   ) :: theta_opt_tile(:,:)

    ! Local Variables
    integer :: k
    integer :: k1Min
    integer :: k2Max
    integer :: kkf
    integer :: kk
    integer :: ij

    ! intermediario para atualizacao de qpcpg_opt_tile e pcprx_mod
    real    :: psfc    
    ! intermediario para calculo de scrmic1_mod
    real    :: colddn0
    ! intermediario para calculo de iemb
    real    :: riemb
    ! indice 
    integer :: lhcat
    ! intermediario para calculo de scrmic2_mod e psfc
    real    :: rolddn0
    ! intermediario para calculo de scrmic2_mod e scrmic3_mod
    real    :: qrolddn0
    ! intermediario para calculo de iemb
    real    :: dispemb
    ! indice dependendo de ij e k
    integer, allocatable :: iemb(:,:)
    ! mascara
    logical, allocatable :: mask1(:,:)

    k1Min = minval(k1_mod(:,lcat))
    k2Max = maxval(k2_mod(:,lcat))
    pcprx_mod(:,lcat) = 0.0

    allocate (iemb(ij_final,k1Min:k2Max))
    allocate (mask1(ij_final,k1Min:k2Max))

    do k = k1Min, k2Max
       do ij =1, ij_final
          mask1(ij,k) =  &
               (k>=k1_mod(ij,lcat)) .and. &
               (k<=k2_mod(ij,lcat)) .and. &
               (rx_mod(ij,k,lcat)>1.e-9)
       end do
    end do

    do k = 2, k2Max !2,k2
       do ij =1, ij_final
          if (k<=k2_mod(ij,lcat)) then
             scrmic1_mod(ij,k) = 0. !rnew(k) = 0.
             scrmic2_mod(ij,k) = 0. !cnew(k) = 0.
             scrmic3_mod(ij,k) = 0. !qrnew(k) = 0.
          end if
       end do
    end do

    do k = k1Min, k2Max
       do ij =1, ij_final
          if (mask1(ij,k)) then
             lhcat = jhcat_mod(ij,k,lcat)

             ! ALF
             ch1(lhcat) = dtlt * cfvt(lhcat) / rtgt_opt_tile(ij)
             !

             dispemb = local_ch1(ij,lhcat)*(emb_mod(ij,k,lcat)/ &
                  cfmas(lhcat)) ** ch3(lhcat) * sqrt(dn0i_mod(ij,k))
             !
             riemb = 1. + ch2(lhcat, ngr) *                         &
                  log10(dispemb / dispemb0(lhcat,ngr))

             !Bob (10/24/00):  Now, limiting iemb to max of nembfall

             iemb(ij,k) = min(nint(riemb),nembfall)

             if (k <= maxkfall) then
                psfc = rx_mod(ij,k,lcat) * dn0_opt_tile(ij,k) * &
                     sfcpcp(k,iemb(ij,k),lhcat)
                qpcpg_opt_tile(ij) = qpcpg_opt_tile(ij) + psfc*qx_mod(ij,k,lcat)
                pcprx_mod(ij,lcat) = pcprx_mod(ij,lcat) + psfc
             end if
          end if
       end do
    end do

    do k = k1Min, k2Max
       do kkf = 1, min(maxkfall,k-1)
          do ij =1, ij_final
             if (mask1(ij,k)) then
                kk = k + 1 - kkf
                lhcat = jhcat_mod(ij,k,lcat)
                colddn0 = cx_mod(ij,k,lcat) * dn0_opt_tile(ij,k)
                rolddn0 = rx_mod(ij,k,lcat) * dn0_opt_tile(ij,k)
                qrolddn0 = qx_mod(ij,k,lcat) * rolddn0
                !cnew
                scrmic1_mod(ij,kk) = scrmic1_mod(ij,kk) + colddn0 *  &
                     dn0i_mod(ij,kk) * pcpfillc(k,kkf,iemb(ij,k),lhcat)
                !rnew
                scrmic2_mod(ij,kk) = scrmic2_mod(ij,kk) + rolddn0 *  &
                     dn0i_mod(ij,kk) * pcpfillr(k,kkf,iemb(ij,k),lhcat)
                !qrnew
                scrmic3_mod(ij,kk) = scrmic3_mod(ij,kk) + qrolddn0 * &
                     dn0i_mod(ij,kk) * pcpfillr(k,kkf,iemb(ij,k),lhcat)
             end if
          end do
       end do
    end do

    pcpg_opt_tile(:) = pcpg_opt_tile(:) + pcprx_mod(:,lcat)
    accpx_mod(:,lcat) = pcprx_mod(:,lcat)
    dpcpg_opt_tile(:) = dpcpg_opt_tile(:) + pcprx_mod(:,lcat) * alphasfc
    pcprx_mod(:,lcat) = pcprx_mod(:,lcat) * dtlti

    do k = 2, k2Max !2,k2
       do ij =1, ij_final
          if (k<=k2_mod(ij,lcat)) then
             rtp_opt_tile(ij,k) = rtp_opt_tile(ij,k) + &
                  scrmic2_mod(ij,k) - rx_mod(ij,k,lcat)
             tairc_mod(ij,k) = tairc_mod(ij,k) - thp_opt_tile(ij,k) * &
                  thp_opt_tile(ij,k) * (2820. * (scrmic2_mod(ij,k) - &
                  rx_mod(ij,k,lcat)) - cpi * (scrmic3_mod(ij,k) - &
                  qx_mod(ij,k,lcat) * rx_mod(ij,k,lcat))) / &
                  (max(tair_mod(ij,k), 253.) * theta_opt_tile(ij,k))
             rx_mod(ij,k,lcat) = scrmic2_mod(ij,k)
             cx_mod(ij,k,lcat) = scrmic1_mod(ij,k)
             qx_mod(ij,k,lcat) = scrmic3_mod(ij,k)/max(1.e-20,scrmic2_mod(ij,k))
             if (rx_mod(ij,k,lcat)<1.e-9) then
                rx_mod(ij,k,lcat) = 0.
                cx_mod(ij,k,lcat) = 0.
                qx_mod(ij,k,lcat) = 0.
             end if
          end if
       end do
    end do
    deallocate (iemb)
    deallocate (mask1)
  end subroutine sedim_opt

  !***************************************************************************

  subroutine final_copy_mic_block(ngr, m1, micro, basic, rtgt, flpw)

    use mem_micro, only:  &
         micro_vars            ! INTENT(OUT)

    use mem_basic, only:  &
         basic_vars            ! INTENT(OUT)

    use micphys, only:  &
         level,         &      ! INTENT(IN)
         irain,         &      ! INTENT(IN)
         ipris,         &      ! INTENT(IN)
         isnow,         &      ! INTENT(IN)
         iaggr,         &      ! INTENT(IN)
         igraup,        &      ! INTENT(IN)
         icloud,        &      ! INTENT(IN)
         ihail                 ! INTENT(IN)

    implicit none

    ! Arguments
    type (micro_vars), intent(out)    :: micro 
    type (basic_vars), intent(out)    :: basic
    integer, intent(in)               :: ngr, m1
    real, intent(out)              :: flpw(:,:)
    real, intent(out)                 :: rtgt(:,:)

    ! Local Variables
    integer :: k, ij

    ! *** 2D Arrays:

    do ij=1,ij_final

       flpw(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))  = real(lpw_opt(ij))
       rtgt(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) =              &
            rtgt_opt_tile(ij)

       if (level >= 3) then
          if(irain >= 1)  then
             micro%accpr(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) = &
                  accpr_tile(ij)
             micro%pcprr(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) = &
                  pcprr_tile(ij)
          endif
          if(ipris >= 1)  then
             micro%accpp(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) = &
                  accpp_tile(ij)
             micro%pcprp(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) = &
                  pcprp_tile(ij)
          endif
          if(isnow >= 1)  then
             micro%accps(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) = &
                  accps_tile(ij)
             micro%pcprs(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) = &
                  pcprs_tile(ij)
          endif
          if(iaggr >= 1)  then
             micro%accpa(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) = &
                  accpa_tile(ij)
             micro%pcpra(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) = &
                  pcpra_tile(ij)
          endif
          if(igraup >= 1) then
             micro%accpg(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) = &
                  accpg_tile(ij)
             micro%pcprg(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) = &
                  pcprg_tile(ij)
          endif
          if(ihail >= 1)  then
             micro%accph(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) = &
                  accph_tile(ij)
             micro%pcprh(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) = &
                  pcprh_tile(ij)
          endif
          micro%pcpg(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))  =    &
               pcpg_opt_tile(ij)
          micro%qpcpg(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) =    &
               qpcpg_opt_tile(ij)
          micro%dpcpg(indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) =    &
               dpcpg_opt_tile(ij)
       endif

    enddo

    ! *** 3D Arrays:

    do k = 1,m1

       do ij=1,ij_final

          ! Input data used in thrmstr_opt
          basic%pp(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))   =   &
               pp_opt_tile(ij,k)
          basic%thp(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))   =  &
               thp_opt_tile(ij,k)
          basic%rtp(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))   =  &
               rtp_opt_tile(ij,k)
          basic%theta(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) =  &
               theta_opt_tile(ij,k)
          basic%rv(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))    =  &
               rv_opt_tile(ij,k)
          ! each_column
          basic%dn0(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))   =  &
               dn0_opt_tile(ij,k)
          if (level >= 2 ) then
             micro%rcp(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr)) = &
                  rcp_tile(ij,k)
          endif
          if (level >= 3) then
             if(irain >= 1) then
                micro%rrp(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                     rrp_tile(ij,k)
                micro%q2(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                     q2_tile(ij,k)
             endif
             if(ipris >= 1)  then
                micro%rpp(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                     rpp_tile(ij,k)
             endif
             if(isnow >= 1)  then
                micro%rsp(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                     rsp_tile(ij,k)
             endif
             if(iaggr >= 1)  then
                micro%rap(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                     rap_tile(ij,k)
             endif
             if(igraup >= 1) then
                micro%rgp(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                     rgp_tile(ij,k)
                micro%q6(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                     q6_tile(ij,k)
                micro%rhp(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                     rhp_tile(ij,k)
                micro%q7(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                     q7_tile(ij,k)
             endif
             if(icloud == 5) then
                micro%ccp(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                     ccp_tile(ij,k)
             endif
             if(irain == 5) then
                micro%crp(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                     crp_tile(ij,k)
             endif
             if(ipris == 5) then
                micro%cpp(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                     cpp_tile(ij,k)
             endif
             if(isnow == 5) then
                micro%csp(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                     csp_tile(ij,k)
             endif
             if(iaggr == 5) then
                micro%cap(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                     cap_tile(ij,k)
             endif
             if(igraup == 5) then
                micro%cgp(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                     cgp_tile(ij,k)
             endif
             if(ihail == 5) then
                micro%chp(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                     chp_tile(ij,k)
             endif
             micro%cifnp(k,indicei(ij,indextile,ngr),indicej(ij,indextile,ngr))=&
                  cifnp_tile(ij,k)
          endif

       enddo

    enddo

  end subroutine final_copy_mic_block

  !****************************************************************************

  subroutine copy_ijk2k(ka, kz, ij, source_ijk, dest_k)

    implicit none

    ! Arguments
    integer, intent(in) :: ka, kz, ij
    real, intent(in)    :: source_ijk(:,:)
    real, intent(out)   :: dest_k(:)

    ! Local Variables
    integer :: k

    do k = ka,kz

       dest_k(k) = source_ijk(ij,k)

    enddo

  end subroutine copy_ijk2k

  !**************************************************************************

end Module mem_micro_opt
