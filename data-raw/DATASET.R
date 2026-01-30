## code to prepare `DATASET` dataset goes here

chl <- 100*runif(10)
car <- 25*runif(10)
ant <- 2*runif(10)
ewt <- 0.04*runif(10)
lma <- 0.02*runif(10)
n_struct <- 1+2*runif(10)
input_prospect <- data.frame('chl' = chl, 'car' = car, 'ant' = ant,
                             'ewt' = ewt, 'lma' = lma, 'n_struct' = n_struct)

usethis::use_data(input_prospect, internal = FALSE, overwrite = TRUE)
