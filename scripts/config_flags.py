##
# Configure the flags for C++ compiler
##

def configFlags(sc_env):
    # check if we are building for debug
    if sc_env['debug']:
        # add debug only flags
        sc_env.Append(CCFLAGS = '-g -Wall')
        sc_env.Append(CCFLAGS = '-D_DO_RANGE_CHECK_')
        conf = sc_env.Configure(config_h = 'config.h')
        conf.Define('_DO_RANGE_CHECK_')
        conf.Define('HAVE_DEBUG')
        conf.Finish()
    else:
        # add optimize flags
        sc_env.Append(CCFLAGS = '-O3 -funroll-loops -Wall')
