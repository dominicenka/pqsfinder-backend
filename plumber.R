# plumber.R
library(pqsfinder)
library(digest)

#* @filter cors
cors <- function(req, res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  
  if (req$REQUEST_METHOD == "OPTIONS") {
    res$setHeader("Access-Control-Allow-Methods", "*")
    res$setHeader("Access-Control-Allow-Headers",
                  req$HTTP_ACCESS_CONTROL_REQUEST_HEADERS)
    res$status <- 200
    return(list())
  } else {
    plumber::forward()
  }
}


# Get range limits of pqsfinder options
#
pqsfinder_limits <- function() {
  return(list(
    max_sequence_len = c(0,1e5),
    max_len = c(10,100),
    min_score = c(1,NA),
    max_defects = c(0,3),
    loop_min_len = c(0,NA),
    loop_max_len = c(0,NA),
    max_bulges = c(0,3),
    max_mismatches = c(0,3),
    run_min_len = c(2,NA),
    run_max_len = c(2,NA)
  ))
}

# Get job root directory
#
get_job_root_path <- function() {
  return("../results")
}

# Get job directory
#
get_job_dir_path <- function(id) {
  return(sprintf("%s/%s", get_job_root_path(), substring(id, 1, 2)))
}

# Get job file
#
get_job_file_path <- function(id) {
  return(sprintf("%s/%s.fa", get_job_dir_path(id), id))
}

# Generate short unique id in a concurrent safe way (extremely low chance of conflicts)
#
get_unique_job_id <- function() {
  while (TRUE) {
    unique_hash <- digest(list(Sys.time(), Sys.getpid()), "sha256")
    job_id <- substring(unique_hash, 1, 6)
    job_file <- get_job_file_path(job_id)
    if (!file.exists(job_file)) {
      return(job_id)
    }
  }
}

# Validate ranges of pqsfinder options
#
validate_pqsfinder_opt_range <- function(opts, limits) {
  validated_opts <- list()
  for (opt in names(opts)) {
    if (opt %in% names(limits)) {
      value = as.numeric(opts[[opt]])
      min_value <- limits[[opt]][1]
      max_value <- limits[[opt]][2]
      
      if (!is.na(min_value) && value < min_value) {
        stop(sprintf("pqsfinder option %s is lower than minimal value %s", opt, min_value))
      }
      if (!is.na(max_value) && value > max_value) {
        stop(sprintf("pqsfinder option %s is greater than maximal value %s", opt, max_value))
      }
      validated_opts[[opt]] <- value
    }
  }
  return(validated_opts)
}

#* Return the results of pqsfinder
#* @get /job/<id>
function(id, res) {
  res_file <- get_job_file_path(id)
  
  if (!file.exists(res_file)) {
    return(0)
  }
  include_file(res_file, res)
}

#* Return default pqsfinder parameters
#* @serializer unboxedJSON
#* @get /formals
function() {
  return(formals(pqsfinder)[-1])
}

#* Return pqsfinder limits
#* @serializer unboxedJSON
#* @get /limits
function() {
  return(pqsfinder_limits())
}

#* Return current version of pqsfinder
#* @serializer unboxedJSON
#* @get /version
function() {
  return(as.character(packageVersion("pqsfinder")))
}

#* Get options and dna strings, return job ID ... example:
#* @param opts Options to configure algorithm
#* @param sequences The data on which to look for quadruplexes
#* @post /analyze
function(opts, sequences) {
  limits <- pqsfinder_limits()
  call_args <- validate_pqsfinder_opt_range(opts, limits)
  call_args$strand <- opts[['strand']]
  
  job_id <- get_unique_job_id()
  job_dir <- get_job_dir_path(job_id)
  dir.create(job_dir)
  
  for (i in seq_len(nrow(sequences))) {
    seq_string <- sequences$seq_string[[i]]
    call_args$subject <- tryCatch(DNAString(seq_string), error = function(e) RNAString(seq_string))
    pqs <- do.call(pqsfinder, call_args)
    
    job_file <- get_job_file_path(job_id)
    
    if (length(pqs) == 0) {
      info <- sprintf("%s", sequences$seq_description[[i]])
      write(c("*", seq_string, length(pqs), info), file = job_file, append = TRUE)
    }
    else {
      if (class(subject(pqs)) == "DNAString") {
        seq_set <- as(pqs, "DNAStringSet")
      } else {
        seq_set <- as(pqs, "RNAStringSet")
      }
      names(seq_set) <- sprintf("%s;%s", sequences$seq_description[[i]], names(seq_set))
      write(c("*", seq_string, length(pqs)), file = job_file, append = TRUE)
      writeXStringSet(seq_set, file = job_file, format = "fasta", append = TRUE)
    }
  }
  return(job_id)
}
