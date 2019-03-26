Experiment data
----

The experiment data for the endotoxin and kidney studies are location in `raw_data/`. The file `filter_probes.rds` contains the probes used in the paper (see `../analysis/filter_probes.R`). 

Analyzed experiments
-----

These files are the combined results using F-test, moderated F-test, bootstrap F-test, bootstrap moderated F-test, and the optimal discovery procedure:

- `dose_comp.rds`: Dose study
- `kidney_comp.rds`: Kidney study
- `endotoxin_comp.rds`: Endotoxin study
- `smoker_comp.rds`: Smoker study

Simulation study
-----

The simulation results are contained in `simulution_results.rds`.

Intermediate files
-----

If the simulation code is run then there should be a directory `raw_data/rdafiles/` and other objects created in `raw_data/`. These files are excluded due to space constraints. However, running the code in `../simulation/` will produce these files.