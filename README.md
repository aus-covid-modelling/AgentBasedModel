# Agent Based Model framework utilised for National Cabinet

# Compiling and running
This code relies on C++14. A `makefile` is included for compilation using `gnu make`.

To run the code, you will need to create a JSON parameters file. `NationalVaccineModel/parameters/dev/` has an example of the file and includes all relevant inputs.

Once compiled, the code can be run with `./Run parameters/dev/parameters-dev.json`, or with whatever JSON file you have created. Line-listed outputs for infected cases will appear in the `output_directory`, with 1 file per simulation.

## Combining outputs
The R script `combine_outputs.R` will take the output files and combine them into a single output file for further processing. It takes one argument, the parameter file JSON, and will calculate everything else required. The combined output will appear in the higher level folder.

## Generating HPC parallel jobs
The R script `create_batch_parameters.R` will generate parameter JSONs, and slurm job scripts for submitting to the SLURM HPC queueing system. It takes two arguments: the path to the 'master' JSON file and the desired number of iterations. The slurm scripts will appear in the current working directory.