/**
 * @file	lcgridode.cxx
 *
 * @brief	Test for grid ODE solver.
 */

// lucee includes
#include <LcLuaModuleRegistry.h>
#include <LcLuaState.h>
#include <LcLuaTable.h>
#include <LcObjRegistry.h>
#include <LcPointSourceIfc.h>
#include <LcTest.h>

class HarmonicOde : public Lucee::PointSourceIfc
{
  public:
/** Class id: this is used by registration system */
    static const char *id;

    HarmonicOde()
      : Lucee::PointSourceIfc(2, 2) 
    {
      gamma = 5.0;
    }

    ~HarmonicOde()
    {
    }

    void
    getSource(const double loc[3], std::vector<double>& src)
    {
      src[0] = this->getData(1);
      src[1] = -gamma*this->getData(0);
    }

  private:
    double gamma;
};
const char *HarmonicOde::id = "Harmonic";

int
main(void)
{
  LC_BEGIN_TESTS("lcgridode");

  Lucee::LuaState L;
// register stuff
  new Lucee::ObjRegistry<Lucee::PointSourceIfc, HarmonicOde>;
  Lucee::LuaModuleRegistry<Lucee::PointSourceIfc>::registerModule(L);

// create Lua table for source
  std::string tblStr = "top = {"
    "harmonic = PointSource.Harmonic { inpComponents = {0, 1}, outComponents = {0, 1} }"
    "}";
// evaluate string as Lua code
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "top");
// construct LuaTable object
  Lucee::LuaTable top(L, "top");

// get hold of harmonic object
  Lucee::PointSourceIfc& point = top.getObjectAsBase<Lucee::PointSourceIfc>("harmonic");

  double loc[3] = {0.0, 0.0, 0.0};
  double inp[2] = {1.0, 2.0};
  double src[2];
// compute source
  point.calcSource(loc, inp, src);
// test it
  LC_ASSERT("Testing source", src[0] == inp[1]);
  LC_ASSERT("Testing source", src[1] == -5*inp[0]);

  LC_END_TESTS;
}
