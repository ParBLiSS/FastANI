AC_INIT(src/cgi/core_genome_identity.cpp)

AC_ARG_WITH(gsl, [  --with-gsl=<path/to/gsl>     GNU Scientific Library install dir (default: /usr/local/)])
AC_ARG_ENABLE(static-gsl, [  --enable-static-gsl])

if test "$with_gsl" == ""
then
	with_gsl=/usr/local/
fi

AC_ARG_WITH(boost, [  --with-boost=<path/to/boost>     Boost Library install dir (will be used instead of GSL)])

AC_LANG(C++)

AC_OPENMP

AC_CHECK_HEADER(zlib.h, [result=1], [result=0])

if test $result == 0
then
	AC_MSG_ERROR([Zlib not found.])
fi


CPPFLAGS="-std=c++11"

if test "$with_boost" == ""
then
    CPPFLAGS="-I$with_gsl/include"
    AC_CHECK_HEADER(gsl/gsl_cdf.h, [result=1], [result=0])

    if test $result == 0
    then
        AC_MSG_ERROR([GNU Scientific Library headers not found.])
    fi
    
    if test "x$enable_static_gsl" == "xyes"
    then
        AC_SUBST(mathlib, "$with_gsl/lib/libgsl.a $with_gsl/lib/libgslcblas.a")
    else
        AC_SUBST(mathlib, "-L$with_gsl/lib -lgsl -lgslcblas")
    fi
    
    AC_SUBST(mathinc, $with_gsl/include)
else
    CPPFLAGS="-I$with_boost/include"
    AC_CHECK_HEADER(boost/math/distributions/binomial.hpp, [result=1], [result=0])

    if test $result == 0
    then
        AC_MSG_ERROR([Boost Library headers not found.])
    fi
    
    AC_SUBST(mathlib, $with_boost/lib/libboost_math_c99.a)
    AC_SUBST(mathinc, $with_boost/include)
    AC_SUBST(amcppflags, "-DUSE_BOOST")
fi

AC_OUTPUT(Makefile)
