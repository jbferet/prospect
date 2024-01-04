test_that("PROSPECT-D inversion produces accurate biophysical assessment", {
  BPinit <- list(N = 1.3, CHL = 20, CAR = 5, ANT = 0.2,
                 EWT = 0.015, LMA = 0.008)
  lrt <- PROSPECT(N = BPinit$N, CHL = BPinit$CHL, CAR = BPinit$CAR,
                  ANT = BPinit$ANT, EWT = BPinit$EWT, LMA = BPinit$LMA)
  leafBP <- Invert_PROSPECT(Refl = lrt$Reflectance,
                            Tran = lrt$Transmittance)
  expect_true(abs(leafBP$N-BPinit$N)<1e-5)
  expect_true(abs(leafBP$CHL-BPinit$CHL)<1e-3)
  expect_true(abs(leafBP$CAR-BPinit$CAR)<1e-4)
  expect_true(abs(leafBP$ANT-BPinit$ANT)<1e-4)
  expect_true(abs(leafBP$EWT-BPinit$EWT)<1e-6)
  expect_true(abs(leafBP$LMA-BPinit$LMA)<1e-6)

  lrt$Reflectance <- lrt$Reflectance*(1+rnorm(length(lrt$Reflectance),
                                              0,0.01))
  lrt$Transmittance <- lrt$Transmittance*(1+rnorm(length(lrt$Transmittance),
                                                  0,0.01))

  leafBP_noise <- Invert_PROSPECT(Refl = lrt$Reflectance,
                                  Tran = lrt$Transmittance)
  expect_true(abs(leafBP_noise$N-BPinit$N)<1e-3)
  expect_true(abs(leafBP_noise$CHL-BPinit$CHL)<1e-1)
  expect_true(abs(leafBP_noise$CAR-BPinit$CAR)<1e-1)
  expect_true(abs(leafBP_noise$ANT-BPinit$ANT)<1e-1)
  expect_true(abs(leafBP_noise$EWT-BPinit$EWT)<1e-4)
  expect_true(abs(leafBP_noise$LMA-BPinit$LMA)<1e-4)
})
