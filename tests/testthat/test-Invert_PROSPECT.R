test_that("PROSPECT-D inversion produces accurate biophysical assessment", {
  bp_init <- list(n_struct = 1.3, chl = 20, car = 5,
                  ant = 0.2, ewt = 0.015, lma = 0.008)
  lrt <- prospect(n_struct = bp_init$n_struct,
                  chl = bp_init$chl,
                  car = bp_init$car,
                  ant = bp_init$ant,
                  ewt = bp_init$ewt,
                  lma = bp_init$lma)
  leaf_bp <- invert_prospect(refl = lrt$reflectance,
                             tran = lrt$transmittance)

  expect_true(abs(leaf_bp$n_struct - bp_init$n_struct)<1e-5)
  expect_true(abs(leaf_bp$chl - bp_init$chl)<1e-3)
  expect_true(abs(leaf_bp$car - bp_init$car)<1e-3)
  expect_true(abs(leaf_bp$ant - bp_init$ant)<1e-4)
  expect_true(abs(leaf_bp$ewt - bp_init$ewt)<1e-6)
  expect_true(abs(leaf_bp$lma - bp_init$lma)<1e-6)

  lrt$reflectance <- lrt$reflectance*(1+rnorm(length(lrt$reflectance),
                                              0,0.01))
  lrt$transmittance <- lrt$transmittance*(1+rnorm(length(lrt$transmittance),
                                                  0,0.01))
  leaf_bp_noise <- invert_prospect(refl = lrt$reflectance,
                                   tran = lrt$transmittance)
  expect_true(abs(leaf_bp_noise$n_struct - bp_init$n_struct)<1e-3)
  expect_true(abs(leaf_bp_noise$chl - bp_init$chl)<1e-1)
  expect_true(abs(leaf_bp_noise$car - bp_init$car)<1e-1)
  expect_true(abs(leaf_bp_noise$ant - bp_init$ant)<1e-1)
  expect_true(abs(leaf_bp_noise$ewt - bp_init$ewt)<1e-4)
  expect_true(abs(leaf_bp_noise$lma - bp_init$lma)<1e-4)
})
