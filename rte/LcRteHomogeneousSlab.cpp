/**
 * @file	LcRteHomogeneousSlab.cpp
 *
 * @brief	Radiative transfer equation in homogeneous slab.
 *
 * @version	$Id: LcRteHomogeneousSlab.cpp 331 2010-03-10 03:55:02Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcRteHomogeneousSlab.h>

namespace Lucee
{
  const char *RteHomogeneousSlab::id = "RteHomogeneousSlab";

  RteHomogeneousSlab::RteHomogeneousSlab()
    : SolverIfc(RteHomogeneousSlab::id)
  {
  }

  RteHomogeneousSlab::~RteHomogeneousSlab()
  {
  }

  void
  RteHomogeneousSlab::readInput(Lucee::LuaTable& tbl)
  {
  }

  void 
  RteHomogeneousSlab::buildData()
  {
  }

  void 
  RteHomogeneousSlab::buildAlgorithms()
  {
  }

  void
  RteHomogeneousSlab::initialize()
  {
  }

  int
  RteHomogeneousSlab::advance(double t)
  {

    return 0;
  }

  void
  RteHomogeneousSlab::writeToFile(const std::string& baseName) const
  {
  }

  void
  RteHomogeneousSlab::restoreFromFile(const std::string& baseName)
  {
  }

  void 
  RteHomogeneousSlab::finalize()
  {
  }
}
