AC_INIT([ED2],[m4_esyscmd_s([sed -n 's/.*ED_VERSION=\(.*\)/\1/p' ED/build/make/paths.mk])])
AM_INIT_AUTOMAKE

MAIN="ed-${PACKAGE_VERSION}"
os_name=$(uname | sed "s/-.*//")
SHLIB="rED2.so"

CFLAGS="${CFLAGS} -I\"ED/src/include\""

PKG_PROG_PKG_CONFIG()
AS_IF([test "x$os_name" != "xDarwin" -a "x$os_name" != "xLinux"],
  [
    AC_DEFINE_UNQUOTED([PKG_CONFIG_SYSTEM_INCLUDE_PATH], ["$PKG_CONFIG_SYSTEM_INCLUDE_PATH"], [MINGW defined pkg-config default system path])
    CFLAGS="${CFLAGS} -I\"${PKG_CONFIG_SYSTEM_INCLUDE_PATH}\" ${FCFLAGS} -D_WIN32"
    CPPFLAGS="${CFLAGS} ${FCFLAGS} -D_WIN32"
    MAIN="${MAIN}.exe"
    SHLIB="rED2.dll"
  ]
)


AC_PROG_FC([ifx ifort mpif90 gfortran])
AC_PROG_CC
AC_PROG_RANLIB
AC_LANG([Fortran])


PKG_CHECK_MODULES([hdf5_fortran], [hdf5_fortran], [have_hdf5=yes], [have_hdf5=no])

AS_IF([test "x$FC" = "xmpif90"],
  [
    AS_IF([test "x$os_name" != "xDarwin" -a "x$os_name" != "xLinux"],
      [
        PKG_CHECK_MODULES([MPI], [msmpi], [have_mpi=yes], [have_mpi=no])
      ],
      [
        PKG_CHECK_MODULES([MPI], [mpi-fort], [have_mpi=yes], [have_mpi=no])
      ]
    )

    AS_IF([test "x$have_mpi" = "xyes"],
      [
        CFLAGS="${CFLAGS} ${MPI_CFLAGS} -DRAMS_MPI -fallow-invalid-boz -fallow-argument-mismatch"
        LIBS="${LIBS} ${MPI_LIBS}"
      ],
      [
        AC_PROG_FC([gfortran])
      ]
    )
  ]
)


AS_IF([test "x$have_hdf5" = "xyes"],
  [
    LIBS="${LIBS} ${hdf5_fortran_LIBS} ${MPILIBS}"
    CFLAGS="${CFLAGS} ${hdf5_fortran_CFLAGS}"
  ],
  [
    PKG_CHECK_MODULES([HDF5], [hdf5], [have_hdf5=yes], [have_hdf5=no])
    
    AS_IF([test "x$have_hdf5" = "xyes"],
      [
        LIBS="${LIBS} ${HDF5_LIBS} -lhdf5 -lhdf5_fortran"
        FCFLAGS="${CFLAGS} ${HDF5_CFLAGS}"

        AC_MSG_CHECKING([for HDF5 support in Fortran])
        AC_LANG_PUSH([Fortran])

        AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
              use hdf5
              integer :: hdferr
              call h5open_f(hdferr)
              call h5close_f(hdferr)
        ]])],
        [AC_MSG_RESULT([yes])
        HDF5_OK=yes],
        [AC_MSG_RESULT([no])
        HDF5_OK=no])

        AC_LANG_POP([Fortran])
      ],
      [have_hdf5=no])
  ]
)

AS_IF([test "x$have_hdf5" = "xno"],[
      AC_MSG_ERROR([
HDF5 library not found!

Please install HDF5 library and try again:

    Debian: sudo apt-get install libhdf5-dev
    Fedora: sudo dnf install hdf5-devel
    CentOS: sudo yum install hdf5-devel
    MacOS: brew install hdf5
    Windows MinGW: pacman -S mingw-w64-x86_64-hdf5
      ])
])


AC_SUBST([CC])
AC_SUBST([FC])
AC_SUBST([SHLIB])
AC_SUBST([MAIN])
AC_SUBST([LIBS])
AC_SUBST([CFLAGS])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT