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
#    ('CC', 'The C compiler to use', 'gcc'),
#    ('CXX', 'The C++ compiler to use', 'g++'),
    ('EXTRA_CCFLAGS', 'Extra compiler flags to pass to build', ''),
    ('EXTRA_LINKFLAGS', 'Extra link flags to pass to build', ''),
)

# HDF5 serial library setup
import scripts.config_hdf5s as config_hdf5s
config_hdf5s.bc.begin(opts)

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

if env['parallel']:
    # configure MPI if needed    
    if config_mpi.bc.conf(env):
        config_mpi.bc.finish(myEnv, addIncs=False, addLibs=False)
        # set C++ compiler to mpicxx
        mpicxx = config_mpi.bc.getBinWithPath()
        myEnv['CXX'] = mpicxx
        env['CXX'] = mpicxx        
        # determine location of mpicc file
        sp = mpicxx.split('/')
        sp[-1] = 'mpicc'
        mpicc = string.join(sp, '/')
        # set C compiler to mpicc
        myEnv['CC'] = mpicc
        env['CC'] = mpicc

        # add flag to work around MPICH bug        
        myEnv.Append(CCFLAGS = '-DMPICH_IGNORE_CXX_SEEK')
    else:
        print "Parallel build requested, but mpi not found"
        Exit(1)

# configure HDF5: we need different configurations depending on serial
# or parallel
if env['parallel']:
    if config_hdf5p.bc.conf(env):
        config_hdf5p.bc.finish(myEnv)
        if config_hdf5p.testNewIfc(config_hdf5p.bc.getIncPath()):
            # for some strange reason HDF5 interface for version great
            # than 1.6.4 is different
            myEnv.Append(CCFLAGS = '-DH5_HAVE_PARALLEL')
            myEnv.Append(CCFLAGS = '-DNEW_H5S_SELECT_HYPERSLAB_IFC')
            myEnv.Append(LIBS = 'z') # I DO NOT KNOW HOW THIS WORKS
    else:
        print "Parallel build needs parallel HDF5, which was not found"
        Exit(1)
else:
    if config_hdf5s.bc.conf(env):
        config_hdf5s.bc.finish(myEnv)
        if config_hdf5s.testNewIfc(config_hdf5s.bc.getIncPath()):
            # for some strange reason HDF5 interface for version great
            # than 1.6.4 is different
            myEnv.Append(CCFLAGS = '-DNEW_H5S_SELECT_HYPERSLAB_IFC')
            myEnv.Append(LIBS = 'z') # I DO NOT KNOW HOW THIS WORKS            
    else:
        print "Serial build needs serial HDF5, which was not found"
        Exit(1)

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

##
# build loki library
##
build_dir = os.path.join(buildin, 'etc/loki')
SConscript('etc/loki/src/SConscript', build_dir=build_dir, duplicate=0)

##
# build mup library
##
build_dir = os.path.join(buildin, 'etc/mup')
SConscript('etc/mup/SConscript', build_dir=build_dir, duplicate=0)

##
# build core library
##
build_dir = os.path.join(buildin, 'lib')
SConscript('lib/SConscript', build_dir=build_dir, duplicate=0)

##
# build test code
##
build_dir = os.path.join(buildin, 'unit')
SConscript('unit/SConscript', build_dir=build_dir, duplicate=0)
