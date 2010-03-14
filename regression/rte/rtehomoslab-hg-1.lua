-- Input file for RTE solution in homogeneous slab

-- top-level simulation object
simulation = Solver.RteHomogeneousSlab {
   -- number of phase function coefficients (degree of anisotropy)
   L = 128,
   -- number of quadrature points in each hemisphere
   N = 32,
   -- cosine of incident angle
   mu0 = 0.5,
   -- beam flux: downward irradiance is mu0*pi*flux
   flux = 1.0,
   -- optical depth of slab
   tau0 = 50.0,
   -- albedo of single scattering
   albedo = 0.9,
   -- number of azimuthal modes
   numModes = 10,
   -- phase function
   phaseFunction = RtePhaseFunction.HG { g = 0.9 },
}
