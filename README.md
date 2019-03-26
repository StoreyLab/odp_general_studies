This contains the code used in [the pre-print](https://doi.org/10.1101/571992).

Overview of directories:

- `simulations/`:
    - `00-analyze-experiments.R`: analysis of the experiments in the paper.
    - `01-combine_objects.R`: merges the datasets from `00-analyze-experiments.R`.
    - `02-main_simulations.R`: contains the steps code used to generate the simulated data in Figure 5. Note that `02-run_simulations.R` is code used to execute the code on the cluster.
    - `03-summarize_simulations.R`: summarizes the simulation results.

- `data/`:
    - Data output from `simulations/`. See subdirectory for descriptions.
    
- `analysis/`:
    - `main_figures.Rmd`: R markdown file to reproduce the figures in the manuscript.
    - `kidney_filter_probes.R`: code used to filter kidney probes.
    - `figures/`: figure used in manuscript (generated from `main_figures.Rmd`).
    
- `manuscript/`: contains manuscript.
    