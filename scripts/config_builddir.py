##
# Configure the name of the build directory
##

def configBuildDir(sc_env):    
    if sc_env['debug']:
        # debug
        if sc_env['parallel']:
            # mpi build
            buildin = 'build-par-deb'
        else:
            # serial build
            buildin = 'build-deb'
    else:
        # optimize
        if sc_env['parallel']:
            # mpi build
            buildin = 'build-par'
        else:
            # serial build
            buildin = 'build'

    return buildin
