module mem_scratch1


  type scratch_vars1

     real, pointer, dimension(:) :: &
          vtu,vtv,vtw,vtp,vtscalar

  end type scratch_vars1

  type (scratch_vars1) :: scratch1

contains
  !---------------------------------------------------------------

  subroutine alloc_scratch1(nodebounds,maxgrds,ngrids,mmzp,mynum)

    use var_tables

    implicit none

    integer :: maxgrds, ngrids, mynum, nodebounds(maxgrds,8)
    integer, dimension (*) :: mmzp
    integer :: nu, nv, nw, np, ns, ng
    integer :: mu, mv, mw, mp, ms
    integer :: ifm, icm, ilf, ilc, jlf, jlc, zlf, zlc


    !         Find the maximum number of grid points needed for any grid.

    nu = 0
    nv = 0
    nw = 0
    np = 0
    ns = 0
    do ng=1,ngrids-1
       icm = ng
       ifm = ng + 1
       ilf = nodebounds(ifm,2) - nodebounds(ifm,1) + 1
       jlf = nodebounds(ifm,4) - nodebounds(ifm,3) + 1
       ilc = nodebounds(icm,6) - nodebounds(icm,5) + 1
       jlc = nodebounds(icm,8) - nodebounds(icm,7) + 1
       zlf = mmzp(ifm) - 2
       zlc = mmzp(icm) - 2
       mu = ilc * jlf * zlf
       mv = ilf * jlc * zlf
       mw = ilf * jlf * zlc
       mp = ilf * jlf * zlf
       ms = mp * num_scalar(ifm)
       nu = max(nu,mu)
       nv = max(nv,mv)
       nw = max(nw,mw)
       np = max(np,mp)
       ns = max(ns,ms)
    enddo

    allocate (scratch1%vtu (nu))
    allocate (scratch1%vtv (nv))
    allocate (scratch1%vtw (nw))
    allocate (scratch1%vtp (np))
    allocate (scratch1%vtscalar (ns))

    return
  end subroutine alloc_scratch1

  !---------------------------------------------------------------

  subroutine nullify_scratch1()

    implicit none

    ! Deallocate all scratch arrays

    if (associated(scratch1%vtu ))  nullify (scratch1%vtu )
    if (associated(scratch1%vtv ))  nullify (scratch1%vtv )
    if (associated(scratch1%vtw ))  nullify (scratch1%vtw )
    if (associated(scratch1%vtp ))  nullify (scratch1%vtp )
    if (associated(scratch1%vtscalar ))  nullify (scratch1%vtscalar )


    return
  end subroutine nullify_scratch1

end module mem_scratch1
