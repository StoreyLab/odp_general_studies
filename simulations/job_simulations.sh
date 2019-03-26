#!/bin/bash

#SBATCH -J 00-tc_sim # a single job name for the array
#SBATCH --qos=1day # which queue
#SBATCH --mem=18000 # mem needed in MB
#SBATCH -o ../data/simulation/00-output%A%a.o # Standard output
#SBATCH -e ../data/simulation/00-output%A%a.e # Standard error
~/tmp/R-3.4.1/builddir/bin/R -f 02-run_simulations.R --args $@
