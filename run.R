library(plumber)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("run.R <instance>")
}
instance_port <- 8000 + as.integer(args[[1]])

r <- plumb("../pqsfinder-backend/plumber.R")
r$run(port = instance_port)
