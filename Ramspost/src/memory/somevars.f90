module somevars
  integer :: myngrids,myn1,myn2,myn3,myjdim,myihtran,mynbig,myistar
  integer, allocatable, dimension(:) :: mynnxp,mynnyp,mynnzp
  real   , allocatable, dimension (:)   :: myplatn,myplonn,mydeltaxn,mydeltayn,mydeltazn
  real   , allocatable, dimension(:,:)  :: mydzmn,mydztn
  real   , allocatable, dimension(:,:)  :: myxmn,myxtn,myymn,myytn,myzmn,myztn
  real   , allocatable, dimension(:,:)  :: myu01dn,myv01dn,mypi01dn,myth01dn,mydn01dn,myrt01dn

  contains 
  
  subroutine alloc_somevars(ngrids,n1,n2,n3)
     implicit none
     integer, intent(in) :: ngrids,n1,n2,n3
     if (.not. allocated(myplatn   )) allocate(myplatn      (ngrids))
     if (.not. allocated(myplonn   )) allocate(myplonn      (ngrids))
     if (.not. allocated(mydeltaxn )) allocate(mydeltaxn    (ngrids))
     if (.not. allocated(mydeltayn )) allocate(mydeltayn    (ngrids))
     if (.not. allocated(mydeltazn )) allocate(mydeltazn    (ngrids))
     if (.not. allocated(mynnxp    )) allocate(mynnxp       (ngrids))
     if (.not. allocated(mynnyp    )) allocate(mynnyp       (ngrids))
     if (.not. allocated(mynnzp    )) allocate(mynnzp       (ngrids))

     if (.not. allocated(myxmn     )) allocate(myxmn     (n1,ngrids))
     if (.not. allocated(myxtn     )) allocate(myxtn     (n1,ngrids))

     if (.not. allocated(myymn     )) allocate(myymn     (n2,ngrids))
     if (.not. allocated(myytn     )) allocate(myytn     (n2,ngrids))

     if (.not. allocated(mydzmn    )) allocate(mydzmn    (n3,ngrids))
     if (.not. allocated(mydztn    )) allocate(mydztn    (n3,ngrids))
     if (.not. allocated(myzmn     )) allocate(myzmn     (n3,ngrids))
     if (.not. allocated(myztn     )) allocate(myztn     (n3,ngrids))
     if (.not. allocated(myu01dn   )) allocate(myu01dn   (n3,ngrids))
     if (.not. allocated(myv01dn   )) allocate(myv01dn   (n3,ngrids))
     if (.not. allocated(mypi01dn  )) allocate(mypi01dn  (n3,ngrids))
     if (.not. allocated(myth01dn  )) allocate(myth01dn  (n3,ngrids))
     if (.not. allocated(mydn01dn  )) allocate(mydn01dn  (n3,ngrids))
     if (.not. allocated(myrt01dn  )) allocate(myrt01dn  (n3,ngrids))
     return
  end subroutine alloc_somevars
end module somevars

