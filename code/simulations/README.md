The script `simulations.R` simulates the community data sets and fits the models for the three scenarios:

* Deterministic community with traits strongly correlated with species
  abundances
* Deterministic community with traits poorly correlated with species
  abundances
* Stochastic community with traits strongly correlated with species
  abundances

A list with the simulated coomunities are stored in `simulated_data.RData`.
For each data set in each the competing models are fit and stored in a list, respectively:

* det_rt.RData
* det_wt2.RData
* sto_rt.RData

The fits to data sets in the last two scenarios had much more convergence errors when ran from the single script (`simulations.R`) than when ran in separate scripts (`simulations_det_wt.R`, `simulations_sto_rt.R`). That is why we have separated scripts for these last two scenarios. Therefore, the simulations.R was used only to simulate the communities and fit the models to the 1st scenario.
