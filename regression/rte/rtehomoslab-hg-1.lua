-- Input file for RTE solution in homogeneous slab

-- start time
tStart = 0.0
-- end time
tEnd = 1.0

-- top-level simulation object
simulation = Solver.RteHomogeneousSlab {
-- number of phase function coefficients
   L = 128,
-- number of quadrature points in each hemisphere
   N = 32,
-- cosine of incident angle
   mu0 = 0.5,
-- optical depth of slab
   tau0 = 50.0,
-- phase function: Henyey-Greenstein
--   phaseFunction = RtePhaseFunction.HG { g = 0.90 },
}