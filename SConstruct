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
    ('EXTRA_CCFLAGS', 'Extra compiler flags to pass to build', ''),
    ('EXTRA_LINKFLAGS', 'Extra link flags to pass to build', ''),
)

# update environment with options
opts.Update(env)
opts.Save('options.cache', env) # save stuff to cache

# generate help messages
Help(opts.GenerateHelpText(env))

# configure compiler flags
import scripts.config_flags as config_flags
config_flags.configFlags(env)

# determine build-directory
import scripts.config_builddir as config_builddir
buildin = config_builddir.configBuildDir(env)

# add flag to indicate we have config.h header
env.Append(CCFLAGS = '-DHAVE_CONFIG_H')

# add the build directory to include path to get hold of config.h header
buildDir = '#%s' % buildin
env.Append(CPPPATH = buildDir)

# clone the environment
myEnv = env.Clone()

# create fresh clones to pass to our sub-builds
env = myEnv.Clone()
parenv = env.Clone()

# add extra compiler and link flags if needed
if env['EXTRA_CCFLAGS'] != '':
    env.Append(CCFLAGS = env['EXTRA_CCFLAGS'])
    parenv.Append(CCFLAGS = env['EXTRA_CCFLAGS'])
if env['EXTRA_LINKFLAGS'] != '':
    env.Append(LINKFLAGS = env['EXTRA_LINKFLAGS'])
    parenv.Append(LINKFLAGS = env['EXTRA_LINKFLAGS'])

# export environments to sub-builds
Export('env')
Export('parenv')

# list of various object files that do registration
reg_objs = []
