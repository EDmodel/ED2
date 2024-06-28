AC_INIT([ED2],[m4_esyscmd_s([sed -n 's/.*ED_VERSION=\(.*\)/\1/p' ED/build/make/paths.mk])])
AM_INIT_AUTOMAKE

MAIN="ed-${PACKAGE_VERSION}"
os_name=$(uname | sed "s/-.*//")
SHLIB="rED2.so"

CFLAGS="${CFLAGS} -IED/src/include"

PKG_PROG_PKG_CONFIG()
AS_IF([test "x$os_name" != "xDarwin" -a "x$os_name" != "xLinux"],
  [
    AC_DEFINE_UNQUOTED([PKG_CONFIG_SYSTEM_INCLUDE_PATH], ["$PKG_CONFIG_SYSTEM_INCLUDE_PATH"], [MINGW defined pkg-config default system path])
    AS_IF([test "x$PKG_CONFIG_SYSTEM_INCLUDE_PATH" != "x"],
    [
      CFLAGS="${CFLAGS} -I\"${PKG_CONFIG_SYSTEM_INCLUDE_PATH}\" ${FCFLAGS} -D_WIN32"
    ])
    MAIN="${MAIN}.exe"
    SHLIB="rED2.dll"
  ]
)



AC_PROG_FC([ifx ifort mpif90 gfortran])
AC_PROG_CC
AC_PROG_RANLIB
AC_LANG([Fortran])

# Check for hdf5_fortran
PKG_CHECK_MODULES([hdf5_fortran], [hdf5_fortran], [
    LIBS="${LIBS} ${hdf5_fortran_LIBS}"
    CFLAGS="${CFLAGS} ${hdf5_fortran_CFLAGS}"
  ], [
    # If hdf5_fortran is not found, check for hdf5
    PKG_CHECK_MODULES([hdf5], [hdf5], [
      # Add hdf5 flags
      LIBS="${LIBS} ${hdf5_LIBS} -lhdf5_fortran"
      CFLAGS="${CFLAGS} ${hdf5_CFLAGS}"
    ], [
      # Try to link directly
      LIBS="${LIBS} -lhdf5_fortran -lhdf5"
    ])
])

# Link check for HDF5 support in Fortran
FCFLAGS="${CFLAGS}"
FFLAGS="${CFLAGS}"
AC_FC_SRCEXT(f90)

AC_MSG_CHECKING([if HDF5 library can link with Fortran code])
AC_LINK_IFELSE([AC_LANG_PROGRAM(, [[
  use hdf5
  integer :: hdferr
    call h5open_f(hdferr)
    call h5close_f(hdferr)
  ]])], [
    AC_MSG_RESULT([yes])
    have_hdf5=yes
  ], [
    AC_MSG_RESULT([no])
    have_hdf5=no
  ]
)

# Case where HDF5 library is not found    
AS_IF([test "x$have_hdf5" = "xno"],[
      AC_MSG_ERROR([
HDF5 library not found!

Please install HDF5 library and try again:

    Debian: sudo apt-get install libhdf5-dev
    Fedora: sudo dnf install hdf5-devel
    CentOS: sudo yum install hdf5-devel
    MacOS: brew install hdf5
    Windows MinGW: pacman -S mingw-w64-ucrt-x86_64-hdf5
      ])
])

# Check for MPI library if mpif90 was found
AS_IF([test "x$FC" = "xmpif90"],
  [
    ORIGINAL_CFLAGS="${CFLAGS}"
    ORIGINAL_LIBS="${LIBS}"
    
    AS_IF([test "x$os_name" != "xDarwin" -a "x$os_name" != "xLinux"],
      [
        PKG_CHECK_MODULES([MPI], [msmpi], [have_mpi=yes], [
          LIBS="${LIBS} -lmsmpi"
        ])
      ],
      [
        PKG_CHECK_MODULES([MPI], [mpi-fort], [have_mpi=yes], [
          LIBS="${LIBS} -lmpi_usempif08 -lmpi_mpifh -lmpi"
        ])
      ]
    )

    AS_IF([test "x$have_mpi" = "xyes"], [
        CFLAGS="${CFLAGS} ${MPI_CFLAGS} -DRAMS_MPI -fallow-invalid-boz -fallow-argument-mismatch"
        LIBS="${LIBS} ${MPI_LIBS}"
      ], [
        # If pkg-config is not found try to link directly
        CFLAGS="${CFLAGS} -DRAMS_MPI -fallow-invalid-boz -fallow-argument-mismatch"

      ]
    )

    FCFLAGS="${CFLAGS}"

    AC_MSG_CHECKING([for MPI can link in Fortran])
    AC_LINK_IFELSE([
      AC_LANG_PROGRAM([], [[
        program test_mpi
            use mpi
            integer :: ierr, rank
            call MPI_Init(ierr)
            call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
            call MPI_Finalize(ierr)
        end program test_mpi
        ]], [[]]
      )], [
          AC_MSG_RESULT([yes])
      ], [
          AC_MSG_RESULT([no])
          CFLAGS="${ORIGINAL_CFLAGS}"
          LIBS="${ORIGINAL_LIBS}"
          AC_PROG_FC([gfortran])
      ]
    )
  ]
)

AC_SUBST([CC])
AC_SUBST([FC])
AC_SUBST([SHLIB])
AC_SUBST([MAIN])
AC_SUBST([LIBS])
AC_SUBST([CFLAGS])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT