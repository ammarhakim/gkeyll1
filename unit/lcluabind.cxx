/**
 * @file	lcluabind.cxx
 *
 * @brief	Unit tests for using LUABIND into Lucee
 */

// lucee includes
#include <LcCmdLineArgs.h>
#include <LcLuaState.h>
#include <LcLuaTable.h>

// luabind
#include <luabind/iterator_policy.hpp>
#include <luabind/luabind.hpp>

// std includes
#include <cstdlib>
#include <iostream>
#include <vector>

void greet()
{
  std::cout << "Hello World" << std::endl;
}

double
square(double x)
{
  return x*x;
}

class LcTest
{
  public:

    double square(double x)
    {
      return x*x;
    }
};

struct Particle
{
    double x, y, z;
    double vx, vy, vz;
};

// particle list
struct ParticleList
{
    std::vector<Particle> p;

    ParticleList(unsigned np)
      : p(np)
    {}
    
    unsigned numPtcls() const
    { return p.size(); }
};

int
init(lua_State *L)
{
  using namespace luabind;
  luabind::open(L);
// register module
  module (L, "Lucee") 
    [
// free-functions
      def("greet", &greet),
      def("square", &square),
// class LcTest
      class_<LcTest>("LcTest")
      .def(constructor<>())
      .def("square", &LcTest::square),
// particle
      class_<Particle>("Particle")
      .def(constructor<>())
      .def_readwrite("x", &Particle::x)
      .def_readwrite("y", &Particle::y)
      .def_readwrite("z", &Particle::z)
      .def_readwrite("vx", &Particle::vx)
      .def_readwrite("vy", &Particle::vy)
      .def_readwrite("vz", &Particle::vz),
// particle list
      class_<ParticleList>("ParticleList")
      .def(constructor<unsigned>())
      .def_readwrite("p", &ParticleList::p, return_stl_iterator)
      .def("numPtcls", &ParticleList::numPtcls)
    ];
}

int
main(int argc, char *argv[])
{
  Lucee::CmdLineArgs cmd("lucee");
  cmd.addArg("i", "INPUT", "Input file");
// parse command line
  cmd.parse(argc, argv);

// check for input file
  if (cmd.hasArg("i") == false)
  {
    cmd.showHelp();
    exit(1);
  }
  std::string inpFile = cmd.getArg("i");

  Lucee::LuaState L;
  init(L);

// read lua script file
  if (luaL_loadfile(L, inpFile.c_str()) || lua_pcall(L, 0, 0, 0))
  {
    std::cerr << "Error parsing input file " << inpFile << std::endl;
    std::string err(lua_tostring(L, -1));
    lua_pop(L, 1);
    std::cerr << err << std::endl;
    exit(1);
  }

  return 0;
}
