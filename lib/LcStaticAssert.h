/**
 * @file	LcMathLib.h
 *
 * @brief	Math library for use in Lucee.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_STATIC_ASSERT_H
#define LC_STATIC_ASSERT_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

namespace Lucee
{
/** A help class to help with compile-time asserts */
  template<bool> class compile_time_check
  {
    public:
      compile_time_check(...) {}
  };

  template<> class compile_time_check<false>
  {
  };
}

/*
 * LcStaticAssert is only in operation when LC_DEBUG is defined. It
 * will test its first argument at compile time and on failure report
 * the error message of the second argument, which must be a valid C++
 * classname. i.e. no spaces, punctuation or reserved keywords.
*/
#ifdef LC_DEBUG
#   define LcStaticAssert(test, errormsg)                       \
    do {                                                        \
      struct ERROR_##errormsg {};                               \
      typedef Lucee::compile_time_check< (test) != 0 > tmplimpl;\
      tmplimpl aTemp = tmplimpl(ERROR_##errormsg());            \
      size_t __lc_x = sizeof(aTemp);                             \
      __lc_x += 1;                                               \
    } while (0)
#else
#   define LcStaticAssert(test, errormsg)
#endif

#endif // LC_STATIC_ASSERT_H
