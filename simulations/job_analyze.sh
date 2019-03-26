#!/bin/bash

#SBATCH -J 00-exp_run # a single job name for the array
#SBATCH --qos=1day # which queue
#SBATCH --mem=18000 # mem needed in MB
#SBATCH -o ../data/simulation/00-output%A%a.o # Standard output
#SBATCH -e ../data/simulation/00-output%A%a.e # Standard error
~/tmp/R-3.4.1/builddir/bin/R -f 00-analyze_experiments.R --args $@
  