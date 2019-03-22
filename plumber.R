# plumber.R

#* @filter cors
cors <- function(req, res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  
  if (req$REQUEST_METHOD == "OPTIONS") {
    res$setHeader("Access-Control-Allow-Methods","*")
    res$setHeader("Access-Control-Allow-Headers", req$HTTP_ACCESS_CONTROL_REQUEST_HEADERS)
    res$status <- 200 
    return(list())
  } else {
    plumber::forward()
  }
}

#* Return the result of pqsfinder() ... example: http://localhost:8000/pqs?id=1
#* @get /job/<id>
function(id, res){
  if (!file.exists(paste("../results/", id, ".fa", sep=""))) return(0);
  include_file(paste("../results/", id, ".fa", sep=""), res)
}

#* Get options and dna strings, return job ID ... example: 
#* @param opts Options to configure algorithm
#* @param sequences The data on which to look for quadruplexes
#* @post /analyze
function(opts, sequences) {
  maxLength <- as.integer(opts['maxLength'])
  minScore <- as.integer(opts['minScore'])
  minLL <- as.integer(opts['minLL'])
  maxLL <- as.integer(opts['maxLL'])
  minRL <- as.integer(opts['minRL'])
  maxRL <- as.integer(opts['maxRL'])
  maxNB <- as.integer(opts['maxNB'])
  maxNM <- as.integer(opts['maxNM'])
  maxND <- as.integer(opts['maxND'])
  strand <- as.character(opts['strand'])
  my_data <- read.delim("../pqsfinder-backend/config.txt")
  id <- my_data[,1]
  i<-0
  for(seqDnaString in sequences['dnaString'][,1]) {
    i<-i+1
    pqs <- pqsfinder(DNAString(seqDnaString), strand = strand ,max_len = maxLength,
                     min_score = minScore, loop_min_len = minLL, loop_max_len = maxLL,
                     run_min_len = minRL, run_max_len = maxRL,
                     max_bulges = maxNB, max_mismatches = maxNM, max_defects = maxND
                     );
    dnaset <- as(pqs, "DNAStringSet")
    names(dnaset) <- sprintf("%s;%s", sequences['seqDescription'][i,1], names(dnaset))
    write(c("*", seqDnaString, length(pqs)), file= paste("../results/", id, ".fa", sep=""),  append=TRUE)
    writeXStringSet(dnaset, file=paste("../results/", id, ".fa", sep=""), format="fasta", append=TRUE)
  }
  for(seqDesc in sequences['seqDescription'][,1]) {
    print(seqDesc);
  }
  id <- id + 1
  write.table(id, file = "../pqsfinder-backend/config.txt", sep = "\t",
              row.names = FALSE)
  return(id - 1)
}

