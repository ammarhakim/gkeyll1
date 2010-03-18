-- Input file for RTE solution in homogeneous slab

-- top-level simulation object
simulation = Solver.RteHomogeneousSlab {
   -- number of phase function coefficients (degree of anisotropy)
   L = 128,
   -- number of quadrature points in each hemisphere
   N = 8,
   -- cosine of incident angle
   mu0 = 0.5,
   -- beam flux: downward irradiance is mu0*pi*flux
   flux = 1.0,
   -- optical depth of slab
   tau0 = 1.0,
   -- albedo of single scattering
   albedo = 0.9,
   -- number of azimuthal modes
   numModes = 10,
   -- phase function
   phaseFunction = RtePhaseFunction.HG { g = 0.8 },
   -- irradiance moments to compute
   irradOut = {0, 1},
   -- optical depths at which outputs are to be computed
   tauOut = {0.0, 1.0},
}
