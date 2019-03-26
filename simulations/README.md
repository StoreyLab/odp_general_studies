Steps to reproduce results
===

The following document outlines the steps to reproduce the simulations results in the paper. This document is designed for specific clusters used at Princeton. See specific analysis files for replicating/integrating into another cluster.

## Analyze studies

The first step is to analyze the experiments using all five testing procedures. The analysis steps are shown in `00-analyze_experiments`. To run the analysis on `gen-comp1` (the cluster at Princeton is called `gen-comp1`, adjust appropriately), 

```
sbatch ./job_analyze.sh kidney
sbatch ./job_analyze.sh endotoxin
sbatch ./job_analyze.sh dose
sbatch ./job_analyze.sh smoker
```

This will save all the `edge` objects to the directory `../data/raw_data/`. Next we need to combine all the q-value objects into a single data frame. This will reduce the size of the objects and make it easier to perform calculations:

```
sbatch --mem=18000 R -f 01-combine_objects.R
```

The following files are created:

- `../data/endotoxin_comp.rds`
- `../data/kidney_comp.rds`
- `../data/dose_comp.rds`
- `../data/smoker_comp.rds`

## Simulations

Now that we have significance results for the studies, we can run the simulations. The location of the file is `02-main_simulations.R`, where we use `02-run_simulations.R` as a wrapper to run the functions in the former file. Running these simulations are slightly more complicated than the previous step: Here we open up a new R session and run the following command:

```
source("02-main_simulations.R")
#timecourse_run(replicates = 1:50, seed = 1)
#timecourse_run(replicates = 51:100, seed = 2)
#timecourse_run(replicates = 101:150, seed = 3)
#timecourse_run(replicates = 151:200, seed = 4)
#timecourse_run(replicates = 201:250, seed = 5)
#timecourse_run(replicates = 251:300, seed = 6)
#timecourse_run(replicates = 301:350, seed = 7)
#timecourse_run(replicates = 351:400, seed = 8)
#timecourse_run(replicates = 401:450, seed = 9)
timecourse_run(replicates = 451:500, seed = 10)
```

This code will run 50 replicates under each simulation setting (limit of 1000 total jobs). To reach an adequate number of simulations per setting, the replicates argument increases along with the seed, e.g., `replicates = 51:100` and `seed = 2`.

Once these simulations finish, the raw files will be located in `/data/raw_data/rdafiles/`. Now the significance results need to be merged for further analysis. The file `03-summarize_simulations.R` contains the code to process these files. To execute this code,

```
sbatch ./job_summarize.sh
```

and a new file, called `../data/simulation_results.rda` will be created. Once you have the simulation files and the analysis files locally, proceed to the final step.

## Creating main figures

The code to create all the main figures in the paper is located at `../analysis/main_figures.Rmd`.
