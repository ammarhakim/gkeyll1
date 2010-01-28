# Includes for this level
include_directories(${PROJECT_SOURCE_DIR}/etc/etc/lua/src)

add_definitions(-DLUA_USE_LINUX)

# Create the library
add_library(lclua
  lapi.c
  ldump.c
  lmem.c
  lstring.c
  lundump.c
  lauxlib.c
  lfunc.c
  loadlib.c
  lstrlib.c
  lvm.c
  lbaselib.c
  lgc.c
  lobject.c
  ltable.c
  lzio.c
  lcode.c
  linit.c
  lopcodes.c
  ltablib.c
  print.c
  ldblib.c
  liolib.c
  loslib.c
  ltm.c
  ldebug.c
  llex.c
  lparser.c
  lua.c
  ldo.c
  lmathlib.c
  lstate.c
  luac.c)

