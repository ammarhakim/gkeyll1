## -*- python -*-
##
# Top level build file: controls how lucee is built.
#
# Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
# Licence version 1.0.
##

import os
import re
import sys
import string

# setup initial construction environment
env = Environment(ENV=os.environ)
#env = Environment()
# construct options: values are stored in options.cache
# defaults are read from config.py
opts = Options(['options.cache', 'config.py'])
# add various options used to build WarpX
opts.AddOptions(
    BoolOption('debug', 'Set to yes to compile for debugging', 'no'),
    BoolOption('parallel', 'Set to yes to compile parallel version', 'no'),
    BoolOption('usepetsc', 'Set to yes to compile with PETSc', 'no'),
)

# update environment with options
opts.Update(env)
opts.Save('options.cache', env) # save stuff to cache

# when on OS X use -Framework Accelerate to pull in LAPACK and BLAS
# symbols
env.AppendUnique(FRAMEWORKS=['Accelerate'])

# generate help messages
Help(opts.GenerateHelpText(env))

# configure compiler flags
import scripts.config_flags as config_flags
config_flags.configFlags(env)

# determine build-directory
import scripts.config_builddir as config_builddir
buildin = config_builddir.configBuildDir(env)

# clone the environment
myEnv = env.Clone()

# CONFIGURATION OF INDIVIDUAL DEPENDENCIES GO BELOW

# add the build directory to include path to get hold of config.h header
buildDir = '#%s' % buildin
env.Append(CPPPATH = buildDir)

# add flag to indicate we have config.h header
if os.path.exists('%s/config.h' % buildDir):
    env.Append(CCFLAGS = '-DHAVE_CONFIG_H')

# create fresh clones to pass to our sub-builds
env = myEnv.Clone()

# export environments to sub-builds
Export('env')

##
# build core library
##
build_dir = os.path.join(buildin, 'lib')
SConscript('lib/SConscript', build_dir=build_dir, duplicate=0)

##
# build core library
##
build_dir = os.path.join(buildin, 'unit')
SConscript('unit/SConscript', build_dir=build_dir, duplicate=0)
