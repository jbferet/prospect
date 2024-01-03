test_that("PROSPECT-D inversion over optimal domains produces accurate biophysical assessment", {
  BPinit <- list(N = 1.3, CHL = 20, CAR = 5, ANT = 0.2,
                 EWT = 0.015, LMA = 0.008)
  lrt <- PROSPECT(N = BPinit$N, CHL = BPinit$CHL, CAR = BPinit$CAR,
                  ANT = BPinit$ANT, EWT = BPinit$EWT, LMA = BPinit$LMA)
  leafBP <- Invert_PROSPECT_OPT(lambda = lrt$wvl,
                                Refl = lrt$Reflectance,
                                Tran = lrt$Transmittance,
                                Parms2Estimate = c('CHL', 'CAR', 'EWT', 'LMA'))
  expect_true(abs(leafBP$CHL-BPinit$CHL)<1e-0)
  expect_true(abs(leafBP$CAR-BPinit$CAR)<1e-1)
  expect_true(abs(leafBP$EWT-BPinit$EWT)<1e-7)
  expect_true(abs(leafBP$LMA-BPinit$LMA)<1e-7)

  lrt$Reflectance <- lrt$Reflectance*(1+rnorm(length(lrt$Reflectance),0,0.01))
  lrt$Transmittance <- lrt$Transmittance*(1+rnorm(length(lrt$Transmittance),0,0.01))
  leafBPopt <- Invert_PROSPECT_OPT(lambda = lrt$wvl,
                                   Refl = lrt$Reflectance,
                                   Tran = lrt$Transmittance,
                                   Parms2Estimate = c('CHL', 'CAR', 'EWT', 'LMA'))
  expect_true(abs(leafBPopt$CHL-BPinit$CHL)<1e-0)
  expect_true(abs(leafBPopt$CAR-BPinit$CAR)<1e-0)
  expect_true(abs(leafBPopt$EWT-BPinit$EWT)<1e-4)
  expect_true(abs(leafBPopt$LMA-BPinit$LMA)<1e-4)
})
