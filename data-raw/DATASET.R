## code to prepare `DATASET` dataset goes here
cdat <- matrix(round(rnorm(19 * 15, 1000, 100)), nrow = 19)

pdat <- matrix(round(rnorm(19 * (15 + 15), 1e5, 5000)), nrow = 19)

usethis::use_data(cdat, pdat, overwrite = TRUE)
