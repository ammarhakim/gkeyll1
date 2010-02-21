/**
 * @file	LcField.cpp
 *
 * @brief	Fields hold multiple values per index location.
 *
 * @version	$Id: LcField.cpp 222 2009-11-17 04:46:22Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>

namespace Lucee
{

// instantiations
  template class Lucee::Field<1, int>;
  template class Lucee::Field<1, float>;
  template class Lucee::Field<1, double>;

  template class Lucee::Field<2, int>;
  template class Lucee::Field<2, float>;
  template class Lucee::Field<2, double>;

  template class Lucee::Field<3, int>;
  template class Lucee::Field<3, float>;
  template class Lucee::Field<3, double>;
}
