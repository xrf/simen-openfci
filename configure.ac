AC_INIT(src/manybody/Slater.hpp)
AM_INIT_AUTOMAKE(openfci,0.6.1)
dnl AC_CONFIG_HEADER(config.h)
AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_PROG_INSTALL
AC_PROG_LIBTOOL

dnl
dnl Process LAPACK options.
dnl 
USER_DEFINED_LAPACK=no
AC_ARG_WITH(lapack,
    [ --with-lapack=<flags> Linker flags to include LAPACK library ],
    [
    AC_MSG_WARN(You are using a custom LAPACK library option.)
    LIBS="$LIBS $with_lapack"
    USER_DEFINED_LAPACK=yes
    ]
)
if test X${USER_DEFINED_LAPACK} = Xno; then

  AC_FC_FUNC(dsyev)
  AC_CHECK_LIB(lapack, $dsyev, [LIBS="-llapack -lblas $LIBS"],
    [AC_MSG_ERROR(You need the LAPACK library. See README.)])

fi



dnl
dnl Process ARPACK options.
dnl
USER_DEFINED_ARPACK=no
AC_ARG_WITH(arpack,
    [ --with-arpack=<flags> Linker flags to include ARPACK library ],
    [
    AC_MSG_WARN(You are using a custom ARPACK library option.)
    dnl ARPACK_FLAGS=$with_arpack
    LIBS="$LIBS $with_arpack"
    USER_DEFINED_ARPACK=yes
    ]
)
if test X${USER_DEFINED_ARPACK} = Xno; then
  AC_FC_FUNC(dseigt)
  AC_CHECK_LIB(arpack, $dseigt, [LIBS="-larpack $LIBS"],
    AC_MSG_ERROR(You need the ARPACK library. See README.))
fi

dnl
dnl Check for ARPACK++ library
dnl
dnl AC_LANG_PUSH(C++)
dnl sk_REQUIRE_HEADER([arpack++/arerror.h], [ArpackError x;])
dnl AC_LANG_POP



dnl
dnl Check for LPP library
dnl
AC_LANG_PUSH(C++)
AC_CHECK_HEADER([lpp/lapack.hh], [],
  AC_MSG_ERROR(You need the LAPACK Plus Plus (LPP) library. See README.))
AC_LANG_POP


dnl
dnl Static linking or not.
dnl
AC_MSG_CHECKING(whether to link statically)
AC_ARG_WITH(static-linking,
[  --with-static-linking   Link binaries statically],
[ case "$withval" in
  no)
     AC_MSG_RESULT(no)
     ;;
  *)
    AC_MSG_RESULT(yes)
    LDFLAGS="-static $LDFLAGS"
    ;;
  esac],
AC_MSG_RESULT(no)
)

dnl
dnl Debugging or not
dnl
AC_MSG_CHECKING(if debbuging should be turned off)
AC_ARG_ENABLE(nodebug,
[  --enable-nodebug    Turn off debugging],
[ 
  AC_DEFINE(NDEBUG,1)
  CXXFLAGS="-O3 $CXXFLAGS" 
  AC_MSG_RESULT(yes)
],
[ 
  CXXFLAGS="-g -O2 $CXXFLAGS" 
  AC_MSG_RESULT(no)
])


dnl
dnl Check for SLATER_BITS
dnl
AC_MSG_CHECKING(number of bits in a Slater determinant)
AC_ARG_WITH(slater-bits, 
  [  --with-slater-bits  Number of bits to use in Slater dets],
  [
    AC_DEFINE_UNQUOTED([SLATER_BITS],[$withval])
    slater_bits="$withval"
  ],
  [
    AC_DEFINE(SLATER_BITS,512)
    slater_bits="512"
  ])
AC_MSG_RESULT([$slater_bits])

AC_SUBST(INCLUDES)

dnl
dnl Output Makefiles
dnl
AC_OUTPUT(Makefile src/Makefile src/ConfigFile/Makefile src/manybody/Makefile src/quantumdot/Makefile src/qdot/Makefile src/examples/Makefile src/arpack/Makefile)

