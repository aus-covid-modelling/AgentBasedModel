library(dplyr)
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

print(max_iter, " copies of ", fname, " created.")
print("Here is what one of them looks like: ")
print(jsonlite::prettify(rjson::toJSON(parfile)))

subJobs <- stringr::str_sub(fname, end=-6L) %>% stringr::str_split("/") 

jobName <- subJobs %>% sapply(., tail, 1)
folderPaths <- subJobs[1] %>% unlist() %>% head(., length(.)-1) %>% paste0(., collapse="/")

jobFname <- paste0(jobName, "job.script")
jobFile <- file(jobFname)

jobScript <- c(
  "#!/bin/bash",
  paste0("#SBATCH --job-name=CommonwealthModelling", jobName),
  "#SBATCH --account=cm37",
  "#SBATCH --time=00:15:00",
  "#SBATCH --ntasks=1",
  "#SBATCH --mem-per-cpu=20480",
  "#SBATCH --cpus-per-task=1",
  "#SBATCH --qos=shortq",
  "#SBATCH --array=1-20",
  "",
  paste0("./Run ", folderPaths, "/", jobName, "_${SLURM_ARRAY_TASK_ID}.json")
)

writeLines(jobScript, con=jobFile)

close(jobFile)
print("Written job file to ", jobFname)


combinerFname <- paste0(jobName, ".combine.job.script")
combinerFile <- file(combinerFname)

jobScript <- c(
  "#!/bin/bash",
  paste0("#SBATCH --job-name=CommonwealthModellingCombiner", jobName),
  "#SBATCH --account=cm37",
  "#SBATCH --time=00:30:00",
  "#SBATCH --ntasks=1",
  "#SBATCH --mem-per-cpu=40960",
  "#SBATCH --cpus-per-task=1",
  "#SBATCH --qos=shortq",
  "",
  paste0("Rscript combine_outputs.R ", folderPaths, "/", jobName)
)

writeLines(jobScript, con=combinerFile)

close(combinerFile)
print("Written combiner file to ", combinerFname)