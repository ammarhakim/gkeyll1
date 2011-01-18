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
  tbl.addValue<double>("gas_gamma", 1.4)
    .setHelp("Gas adiabatic constant")
    .setMinValue(0.0);

  std::vector<std::string> eosType(2);
  eosType[0] = "idealGas";
  eosType[1] = "twoTermGas";
  tbl.addValue<std::string>("eos")
    .setHelp("Equation of state type")
    .setOneOf(eosType);

  LC_RAISES("Testing if non-existent value can be fetched", 
    tbl.getValue<double>("XXX"), Lucee::Except);

  const Lucee::ValueDescription<double>& gasGamma
    = tbl.getValue<double>("gas_gamma");
}

int
main(void)
{
  LC_BEGIN_TESTS("lctabledescription");
  test_1();
  LC_END_TESTS;
}
