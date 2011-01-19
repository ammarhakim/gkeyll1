/**
 * @file	lctabledescripton.cxx
 *
 * @brief	Unit tests table description object.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcTableDescription.h>
#include <LcTest.h>

void
test_1()
{
  Lucee::TableDescription tbl("euler");

  double gamma;
  tbl.addValue<double>("gas_gamma", 1.4)
    .setHelp("Gas adiabatic constant")
    .setMinValue(0.0)
    .setVar(&gamma);

  LC_RAISES("Testing if non-existent value can be fetched", 
    tbl.getValue<double>("XXX"), Lucee::Except);

  const Lucee::ValueDescription<double>& gasGamma
    = tbl.getValue<double>("gas_gamma");

  Lucee::TableDescription grd("grid");
  std::vector<int> cells;
  grd.addVector<int>("cells")
    .setHelp("Number of cells in each direction")
    .setLength(3)
    .setMinValue(1)
    .setVar(&cells);

  std::vector<double> lower;
  grd.addVector<double>("lower")
    .setHelp("Coordinates of lower corner of box")
    .setLength(3)
    .setVar(&lower);

  std::vector<double> upper;
  grd.addVector<double>("upper")
    .setHelp("Coordinates of upper corner of box")
    .setLength(3)
    .setVar(&upper);
}

int
main(void)
{
  LC_BEGIN_TESTS("lctabledescription");
  test_1();
  LC_END_TESTS;
}
