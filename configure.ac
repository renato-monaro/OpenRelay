#
AC_INIT([LibOpenRelay], [0.1], [monaro@usp.br],[openrelay],[http://lsee.sel.eesc.usp.br/orelay/])
#libboost1.46
#libboost-thread1.46
#PKG_CHECK_MODULES([DEPS], [boost >= 1.38 libboost-thread >= 1.38 libc >= 6.0 libncurses >= 5.7 libcomedi >= 0.8])

PKG_CHECK_MODULES([DEPS], [fann >= 0.26])
#AC_SUBST(DEPS_CFLAGS)
#AC_SUBST(DEPS_LIBS)
CPPFLAGS="$CPPFLAGS $DEPS_CFLAGS" 
LIBS="$DEPS_LIBS $LIBS" 
LDFLAGS  = "${LDFLAGS} -shared"


AC_PREREQ([2.59])
AM_INIT_AUTOMAKE([1.10 -Wall no-define])

AC_CONFIG_HEADERS([config.h])
AC_PROG_CXX

#CXXFLAGS="-O1 -g3"
CXXFLAGS="-O3"
#LT_INIT([disable-static])
#LT_INIT([enable-shared])
LT_INIT

AC_SUBST([OPENRELAY_SO_VERSION], [1:2:0])
AC_SUBST([OPENRELAY_API_VERSION], [1.0])

#AC_SUBST([FUZZYLITE_SO_VERSION], [1:2:0])
#AC_SUBST([FUZZYLITE_API_VERSION], [1.0])

AC_ARG_WITH([rtai],
[AS_HELP_STRING([--with-rtai],[RTAI Real-Time Support])],
[],[with_rtai=no])

AC_ARG_WITH([comedi],
[AS_HELP_STRING([--with-comedi],[Comedi Support])],
[],[with_comedi=no])

AC_ARG_WITH([rtnet],
[AS_HELP_STRING([--with-rtnet],[Rtnet Support])],
[],[with_rtnet=no])

AC_ARG_WITH([fuzzy],
[AS_HELP_STRING([--with-fuzzy],[Fuzzylite Support])],
[],[with_fuzzy=no])

AC_ARG_WITH([rna],
[AS_HELP_STRING([--with-rna],[FANN Support])],
[],[with_rna=no])

AC_ARG_WITH([iec],
[AS_HELP_STRING([--with-iec],[IEC61850 Support])],
[],[with_iec=no])

AC_ARG_WITH([ncurses],
[AS_HELP_STRING([--with-ncurses],[NCurses Support])],
[],[with_ncurses=no])

AS_IF(
	[test "x$with_rtai" != xno],
	[
		#Set up RTAI
		CPPFLAGS="-DRTAI -I/usr/lib/realtime/include -I/usr/src/linux/include -Wall -pipe $CPPFLAGS"
	])
AS_IF(	[test "x$with_comedi" != xno],
	[
		#Set up Comedi
		CPPFLAGS="-DCOMEDI $CPPFLAGS"
	])
AS_IF(	[test "x$with_rtnet" != xno],
	[
		#Set up Rtnet
		CPPFLAGS="-DRTNET -I/usr/local/src/rtnet/trunk/rtnet/stack/include -I/usr/local/rtnet/include $CPPFLAGS"
	]
	)
AS_IF(	[test "x$with_fuzzy" != xno],
	[
		#Set up Rtnet
		CPPFLAGS="-DFUZZY $CPPFLAGS"
	]
	)
	
AS_IF(	[test "x$with_fuzzy" != xno],
	[
		#Set up Rtnet
		CPPFLAGS="-DFUZZY $CPPFLAGS"
	]
	)
	
AS_IF(	[test "x$with_rna" != xno],
	[
		#Set up Rtnet
		CPPFLAGS="-DRNA $CPPFLAGS"
	]
	)

AS_IF(	[test "x$with_iec" != xno],
	[
		#Set up Rtnet
		CPPFLAGS="-DIEC61850 -I/usr/local/include/goosesv-1.0 -I/usr/local/lib/goosesv-1.0/include -I/usr/local/include/goosesv-1.0/codecs -I/usr/local/lib/goosesv-1.0/include/codecs  -L/usr/local/lib -lgoosesv-1.0 $CPPFLAGS"
	]
	)
	
AS_IF(	[test "x$with_ncurses" != xno],
	[
		#Set up Rtnet
		CPPFLAGS="-DNCURSES  $CPPFLAGS"
	]
	)
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi
AM_CONDITIONAL([HAVE_DOXYGEN], 
[test -n "$DOXYGEN"])AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([doc/Doxyfile])])

AC_CONFIG_FILES([Makefile
		src/Makefile
		doc/Makefile
                openrelay-${OPENRELAY_API_VERSION}.pc:openrelay.pc.in])
AC_OUTPUT
