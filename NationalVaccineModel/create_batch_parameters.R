library(rjson)

args <- commandArgs(trailingOnly = T)

fname <- paste0(args[1])
parfile <- fromJSON(file=args[1])

max_iter <- as.numeric(args[2])

for (i in 1:max_iter) {
  parfile$start_sims <- i
  
  js <- rjson::toJSON(parfile)
  
  write(jsonlite::prettify(js), file=paste0(gsub(".json", "", fname), "_", i, ".json"))
}


