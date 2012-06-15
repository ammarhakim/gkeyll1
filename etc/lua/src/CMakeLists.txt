# Includes for this level
include_directories(${PROJECT_SOURCE_DIR}/etc/etc/lua/src)

add_definitions(-DLUA_USE_LINUX)

# Create the library
add_library(lclua
  lapi.c
  lauxlib.c
  lbaselib.c
  lcode.c
  ldblib.c
  ldebug.c
  ldo.c
  ldump.c
  lfunc.c
  lgc.c
  linit.c
  liolib.c
  llex.c
  lmathlib.c
  lmem.c
  loadlib.c
  lobject.c
  lopcodes.c
  loslib.c
  lparser.c
  lstate.c
  lstring.c
  lstrlib.c
  ltable.c
  ltablib.c
  ltm.c
  lundump.c
  lvm.c
  lzio.c
  print.c)

# Add executables
#add_executable(lua lua.c)
#target_link_libraries(lua lclua)

#add_executable(luac luac.c)
#target_link_libraries(luac lclua)

