library(pqsfinder)
library(plumber)
r <- plumb("../pqsfinder-backend/plumber.R")
r$run(port=8000)
