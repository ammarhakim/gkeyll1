##
# Class to perform Lucee build configuration.
#
# Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
# Licence version 1.0.
##

import os

class BuildConf:
    r"""BuildConf(pkg : str, base_paths = [], incs : [] str, inc_path : []
    str, libs : [] str, lib_path: [] str)
    
    Class to test existence of libraries and include headers for
    package `pkg`. Here `base_paths` is the list of base directories,
    `incs` is the list of include headers to look for, `inc_path` is
    the list of include paths to look in, `libs` is the list of
    libraries to look for, `lib_path` is the list of library paths to
    look in. If `inc_paths` and `lib_paths` are relative paths
    (i.e. missing a front-slash in the begining) they will be appended
    by the `base_paths` to create the search paths. Thus if
    `base_paths = [/usr/local]` and `inc_paths = ['include',
    '/myincs']`, then the search for include headers will be
    performend in `/usr/local/include` and `/myincs` directories."""
    
    def __init__(self, pkg, base_paths,
                 incs, inc_paths,
                 libs, lib_paths,
                 sub_inc_paths = [], sub_lib_paths = []):
        self.pkg = pkg.upper()
        self.incs = incs
        self.libs = libs
        self.base_paths = base_paths
        self.orig_inc_paths = inc_paths
        self.orig_lib_paths = lib_paths
        self.orig_bin_paths = []
        self.sub_inc_paths = sub_inc_paths
        self.sub_lib_paths = sub_lib_paths

        # construct inc_paths
        self.inc_paths = []
        for ip in inc_paths:
            if ip[0] == '/':
                # this is absolute path so just add it to inc_paths
                p = os.path.expandvars(ip)
                self.inc_paths.append(p)
            else:
                # this is relative path so append base_path
                for bp in base_paths:
                    p = os.path.expandvars(os.path.join(bp,ip))
                    self.inc_paths.append(p)

        # construct lib_paths
        self.lib_paths = []
        for lp in lib_paths:
            if lp[0] == '/':                
                # this is absolute path so just add it to lib_paths
                p = os.path.expandvars(lp)
                self.lib_paths.append(p)
            else:
                # this is relative path so append base_path
                for bp in base_paths:
                    p = os.path.expandvars(os.path.join(bp,lp))
                    self.lib_paths.append(p)

        # for storing results of configuration
        self.good_inc_path = ''
        self.good_lib_path = ''

        self.bin = None
        self.bin_paths = []
        self.bin_with_path = None

    def setBin(self, bin, bin_paths):
        self.bin = bin
        self.orig_bin_paths = bin_paths
        for b in bin_paths:
            if b[0] == '/':                
                # this is absolute path so just add it to bin_paths
                p = os.path.expandvars(b)
                self.bin_paths.append(p)
            else:
                # this is relative path so append base_path
                for bp in self.base_paths:
                    p = os.path.expandvars(os.path.join(bp,b))
                    self.bin_paths.append(p)

    def begin(self, sc_opts):
        r"""begin(sc_opts) -> None

        Adds options to set the base directory, include and link
        directories.
        """

        pkg = self.pkg.lower()
        sc_opts.AddOptions(
            (pkg, 'Directory for %s package' % self.pkg, './'),
            (pkg + '_incdir', 'Include directory for %s package' % self.pkg, './'),
            (pkg + '_libdir', 'Library directory for %s package' % self.pkg, './'),
            )
        if self.bin:
            sc_opts.AddOptions(
                (pkg + '_bindir', 'Executables directory for %s package' % self.pkg, './')
                )

    def getIncPath(self):
        return self.good_inc_path

    def getLibPath(self):
        return self.good_lib_path

    def getBinWithPath(self):
        return self.bin_with_path

    def conf(self, sc_env):
        r"""con(sc_env) -> bool

        Use `sc_env` to configure the package. Returns True if the
        configuration was successful, False otherwise. Once the
        configuration is successful, the include and link paths can be
        gotten by doing getIncPath() and getLibPath() functions.
        """
        # extract paths from environment
        pkg = self.pkg.lower()
        extra_base_dir = sc_env[pkg]

        # extra includes
        extra_inc_path = sc_env[pkg + '_incdir']
        if extra_inc_path != './':
            # include path was specified
            self.orig_inc_paths.append(extra_inc_path)

        # extra library path
        extra_lib_path = sc_env[pkg + '_libdir']
        if extra_lib_path != './':
            # link path was specified
            self.orig_lib_paths.append(extra_lib_path)

        # extra binary path
        if self.bin:
            extra_bin_path = sc_env[pkg + '_bindir']
        else:
            extra_bin_path = './'
        if extra_bin_path != './':
            # bin path was specified
            self.orig_bin_paths.append(extra_bin_path)

        # extend include path
        for ip in self.orig_inc_paths:
            if ip[0] == '/':
                # this is absolute path so just add it to inc_paths
                p = os.path.expandvars(ip)
                if not p in self.inc_paths:
                    self.inc_paths.insert(0, p)
            else:
                # this is relative path so append base_path
                if extra_base_dir != './':
                    p = os.path.expandvars(os.path.join(extra_base_dir,ip))
                    self.inc_paths.insert(0, p)

        # extend library path
        for lp in self.orig_lib_paths:
            if lp[0] == '/':
                # this is absolute path so just add it to lib_paths
                p = os.path.expandvars(lp)
                if not p in self.lib_paths:
                    self.lib_paths.insert(0, p)
            else:
                # this is relative path so append base_path
                if extra_base_dir != './':
                    p = os.path.expandvars(os.path.join(extra_base_dir,lp))
                    self.lib_paths.insert(0, p)

        # extend binary path
        for bp in self.orig_bin_paths:
            if bp[0] == '/':
                # this is absolute path so just add it to bin_paths
                p = os.path.expandvars(bp)
                if not p in self.bin_paths:
                    self.bin_paths.insert(0, p)
            else:
                # this is relative path so append base_path
                if extra_base_dir != './':
                    p = os.path.expandvars(os.path.join(extra_base_dir,bp))
                    self.bin_paths.insert(0, p)
        
        # configure include headers
        allConf = self.confIncs(sc_env)
        # configure libraries
        allConf = allConf and self.confLibs(sc_env)
        # configure executables
        allConf = allConf and self.confBin(sc_env)

        return allConf

    def finish(self, sc_env, addIncs = True, addLibs = True):
        r"""finish(sc_env) -> None

        Complete the configuration by extending the enviornment
        `sc_env` with the paths and libraries.
        """
        if addIncs:
            # add if requested
            sc_env.Append( CPPPATH = [self.getIncPath()] )
        if addLibs:
            # add if requested
            sc_env.Append( LIBPATH = [self.getLibPath()] )
            sc_env.Append( LIBS = self.libs )
        
        sc_env['HAVE_%s' % self.pkg] = True
        if sc_env['parallel']:
            if sc_env['debug']:
                conf = sc_env.Configure(config_h = 'build-par-deb/config.h')
            else:
                conf = sc_env.Configure(config_h = 'build-par/config.h')
        else:
            if sc_env['debug']:
                conf = sc_env.Configure(config_h = 'build-deb/config.h')
            else:
                conf = sc_env.Configure(config_h = 'build/config.h')
            
        conf.Define('HAVE_%s' % self.pkg)
        conf.Finish()

    def confIncs(self, sc_env):
        # test for include headers
        flg = True
        for p in self.inc_paths:
            print "*** BuildConf (%s): Checking headers in %s ...." % (self.pkg, p)
            myEnv = sc_env.Clone()
            myConf = myEnv.Configure()
            myEnv.Append(CPPPATH = [p])
            # add all sub_inc paths before testing
            for si in self.sub_inc_paths:
                sip = os.path.join(p, si)
                myEnv.Append(CPPPATH = [sip])
            # now test if headers exists
            foundAll = True
            for h in self.incs:
                foundAll = foundAll and myConf.CheckCHeader(h)
            myConf.Finish()                    
            if foundAll:
                # found a path which works, so return it
                self.good_inc_path = p
                return True
            else:
                flg = False
        return flg

    def confLibs(self, sc_env):
        # test for libraries
        flg = True
        for p in self.lib_paths:
            print "*** BuildConf (%s): Checking libraries in %s ...." % (self.pkg, p)
            myEnv = sc_env.Clone()
            myConf = myEnv.Configure()
            myEnv.Append(LIBPATH = p)
            # now test if library exists
            foundAll = True
            for l in self.libs:
                foundAll = foundAll and myConf.CheckLib(library=l)
            myConf.Finish()                    
            if foundAll:
                # found a path which works, so return it
                self.good_lib_path = p
                return True
            else:
                flg = False
        return flg

    def confBin(self, sc_env):
        # test for executable
        flg = True
        for p in self.bin_paths:
            print "*** BuildConf (%s): Checking executable in %s ...." % (self.pkg, p)
            # now test if library exists
            ep = sc_env.WhereIs(self.bin, p)
            if ep:
                print "Found binary %s ...." % ep
                self.bin_with_path = ep
                return True
            else:
                flg = False
        return flg
