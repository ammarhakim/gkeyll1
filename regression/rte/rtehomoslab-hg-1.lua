-- Input file for RTE solution in homogenous slav

-- start time
tStart = 0.0
-- end time
tEnd = 1.0

simulation = Solver.RteHomogenousSlab {
-- number of phase function coefficients
   L = 300,
-- Number of quadrature points in each hemisphere
   N = 100,
-- cosine of incident angle
   mu0 = 0.5,
-- optical depth of slab
   tau0 = 50.0,
}