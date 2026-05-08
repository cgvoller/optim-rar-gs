## Repository Map

Below is a file tree, showcasing the structure of this repository.

```text
.
├── main.R
├── config.R
├── priors.R
├── boundaries.R
├── functions/
│   ├── general_fun.R
│   ├── loss_function.R
│   └── ...
├── figures/
│   └── stopping_heatmap.R
├── jennison_benchmarks/
│   ├── 01_tab_no_early_stopping.R
│   ├── 02_tab_no_early_stopping_delay.R
│   ├── 03_tab_early_stopping.R
│   └── ...
└── prior_predictive_loss/
    ├── early_stopping/
    │   ├── grouped/
    │   └── pooled/
    └── no_early_stopping/
        ├── grouped/
        └── pooled/
```

## How to navigate the repository

| Path | Description |
|---|---|
| `main.R` | Main entry point for running the project. |
| `config.R` | Loads packages and custom functions. |
| `priors.R` | Defines priors used in the Bayesian scripts. |
| `boundaries.R` | Defines group sequential boundaries. |
| `functions/` | Contains reusable helper functions used across the simulations. |
| `jennison_benchmarks/` | Contains benchmark scripts based on Jennison’s paper: [“Comment: Group Sequential Designs with Response-Adaptive Randomisation”](https://projecteuclid.org/journals/statistical-science/volume-38/issue-2/Comment-Group-Sequential-Designs-with-Response-Adaptive-Randomisation/10.1214/23-STS865D.pdf). |
| `prior_predictive_loss/` | Contains the main scenario scripts where data is generated under a particular prior. |
| `prior_predictive_loss/early_stopping/` | Home for scenarios with early stopping (Pooled and Grouped). |
| `prior_predictive_loss/no_early_stopping/` | Home for scenarios without early stopping (Pooled and Grouped). |
| `figures/` | Contains scripts for generating plots. |


> Reference: Jennison, C.  
> [“Comment: Group Sequential Designs with Response-Adaptive Randomisation”](https://projecteuclid.org/journals/statistical-science/volume-38/issue-2/Comment-Group-Sequential-Designs-with-Response-Adaptive-Randomisation/10.1214/23-STS865D.pdf),  
> *Statistical Science*, 38(2).

