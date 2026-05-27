## Pathway diversity and ADRIA

Scripts for pathway diversity analysis using ADRIA-CoralBlox.

## Getting started
Create a `config.toml` at the project root:
```
[operation]
num_cores = 1     # No. of cores to use. Values <= 0 will use all available cores.
threshold = 1e-8  # Result values below this will be set to 0.0 (to save disk space)
debug = false     # Disable multi-processing to allow error messages to be shown

[results]
output_dir = "./Outputs"  # Change this to point to where you want to store simulation results

[PathwayDiversityAnalysis]
domain_path      = "<path to domain>"                  # The default analysis uses the GBR_MCB_GBR_2026-03-30_v080 version
coral_param_path = "<path to coral calibrated params>" # Calibrated coral parameters (version from 2026_03)
plot_output_path = "<path to plot outputs>"            # Path where plots and analysis will be saved
```
