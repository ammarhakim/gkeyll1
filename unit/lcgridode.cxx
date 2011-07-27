/**
 * @file	lcgridode.cxx
 *
 * @brief	Test for grid ODE solver.
 */

// lucee includes
#include <LcGridOdePointIntegrator.h>
#include <LcLuaModuleRegistry.h>
#include <LcLuaState.h>
#include <LcLuaTable.h>
#include <LcMathPhysConstants.h>
#include <LcObjRegistry.h>
#include <LcPointSourceIfc.h>
#include <LcRectCartGrid.h>
#include <LcStructGridField.h>
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
      src[1] = -gamma*gamma*this->getData(0);
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
    "harmonic = PointSource.Harmonic { inpComponents = {0, 1}, outComponents = {0, 1} },"
    "}"
    "gridOde = { terms = { top.harmonic } }";
// evaluate string as Lua code
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "top");
// construct LuaTable object
  Lucee::LuaTable top(L, "top");

// fetch table and put on top of stack
  lua_getglobal(L, "gridOde");
// construct LuaTable object
  Lucee::LuaTable gridOde(L, "gridOde");

// get hold of harmonic object
  Lucee::PointSourceIfc& point = top.getObjectAsBase<Lucee::PointSourceIfc>("harmonic");

  double loc[3] = {0.0, 0.0, 0.0};
  double inp[2] = {1.0, 2.0};
  double src[2];
// compute source
  point.calcSource(loc, inp, src);
// test it
  LC_ASSERT("Testing source", src[0] == inp[1]);
  LC_ASSERT("Testing source", src[1] == -25*inp[0]);

// create grid
  int lo[2] = {0};
  int up[2] = {5};
  Lucee::Region<1, int> localBox(lo, up);

  double plo[2] = {0.0};
  double pup[2] = {1.0};
  Lucee::Region<1, double> physBox(plo, pup);
  Lucee::RectCartGrid<1> grid(localBox, localBox, physBox);

// create field
  int lg[1] = {0}, ug[1] = {0};
  Lucee::StructGridField<1, double> fld(&grid, 2, lg, ug);

// set field at t=0 for initial condition
  for (int i=fld.getLower(0); i<fld.getUpper(0); ++i)
  {
    fld(i,0) = 10.0;
    fld(i,1) = 0.0;
  }

// now create a point-integrator class
  Lucee::GridOdePointIntegrator<1> integrator(grid);
  integrator.readInput(gridOde);

  unsigned nstep = 100;
  double tend = 2*Lucee::PI;
// integrate with this time-step
  double dt = tend/nstep, tcurr = 0.0;
  for (unsigned i=0; i<nstep; ++i)
  {
    integrator.integrate(tcurr, tcurr+dt, fld);
    std::cout << tcurr+dt << " " << fld(0,0) << std::endl;
    tcurr = tcurr+dt;
  }

  LC_END_TESTS;
}
