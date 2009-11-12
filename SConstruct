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
opts = Variables(['options.cache', 'config.py'])
# add various options used to build WarpX
opts.AddVariables(
    BoolVariable('debug', 'Set to yes to compile for debugging', 'no'),
    BoolVariable('parallel', 'Set to yes to compile parallel version', 'no'),
    BoolVariable('usegsl', 'Set to yes to compile with GSL', 'no'),
    ('EXTRA_CCFLAGS', 'Extra compiler flags to pass to build', ''),
    ('EXTRA_LINKFLAGS', 'Extra link flags to pass to build', ''),
)

# MPI library setup
import scripts.config_mpi as config_mpi
config_mpi.bc.begin(opts)

# HDF5 serial library setup
import scripts.config_hdf5s as config_hdf5s
config_hdf5s.bc.begin(opts)

# ZLIB library setup
import scripts.config_z as config_z
config_z.bc.begin(opts)

# SZIP library setup
import scripts.config_szip as config_szip
config_szip.bc.begin(opts)

# HDF5 parallel library setup
import scripts.config_hdf5p as config_hdf5p
config_hdf5p.bc.begin(opts)

# GSL library setup
import scripts.config_gsl as config_gsl
config_gsl.bc.begin(opts)

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
parenv = env.Clone()

# CONFIGURATION OF INDIVIDUAL DEPENDENCIES GO BELOW

# configure the zlib package
if config_z.bc.conf(env):
    config_z.bc.finish(myEnv)
else:
    print "Zlib is required for HDF5, but was not found"
    Exit(1)

# configure the szip package
if config_szip.bc.conf(env):
    config_szip.bc.finish(myEnv)
else:
    print "SZIP library not found. Continuing anyway ..."

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
    else:
        print "Serial build needs serial HDF5, which was not found"
        Exit(1)

# configure GSL library if requested
if env['usegsl']:
    if config_gsl.bc.conf(env):
        config_gsl.bc.finish(myEnv)
    else:
        print "GSL requested but not found"
        Exit(1)

# add the build directory to include path to get hold of config.h header
buildDir = '#%s' % buildin
env.Append(CPPPATH = buildDir)

# add flag to indicate we have config.h header
if os.path.exists('%s/config.h' % buildDir):
    env.Append(CCFLAGS = '-DHAVE_CONFIG_H')

# create fresh clones to pass to our sub-builds
env = myEnv.Clone()

# add extra compiler and link flags if needed
if env['EXTRA_CCFLAGS'] != '':
    env.Append(CCFLAGS = env['EXTRA_CCFLAGS'])
if env['EXTRA_LINKFLAGS'] != '':
    env.Append(LINKFLAGS = env['EXTRA_LINKFLAGS'])

# export environments to sub-builds
Export('env')
Export('parenv')

##
# build loki library
##
build_dir = os.path.join(buildin, 'etc/loki')
SConscript('etc/loki/src/SConscript', build_dir=build_dir, duplicate=0)

##
# build lua library
##
build_dir = os.path.join(buildin, 'etc/lua')
SConscript('etc/lua/src/SConscript', build_dir=build_dir, duplicate=0)

##
# build core library
##
build_dir = os.path.join(buildin, 'lib')
SConscript('lib/SConscript', build_dir=build_dir, duplicate=0)

##
# build unit tests
##
build_dir = os.path.join(buildin, 'unit')
SConscript('unit/SConscript', build_dir=build_dir, duplicate=0)
