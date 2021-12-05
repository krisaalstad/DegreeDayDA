# DegreeDayDA
Data assimilation in a degree-day snowmelt model. Currently this is just a synthetic experiment for two grid cells near Finse in the 2019 water year.

A degree day model takes daily temperature (in degrees Celsius) and precipitation (in mm/day) as inputs to simulate the evolution of the snow water equivalent (SWE). Snow accumulation is diagnosed using a fixed (assumed known) threshold temperature for snowfall and by applying a snowfall multiplier (uncertain parameter) to the precipitation to account for potential biases and unresolved proceses like snowdrift. Snow ablation is diagnosed for positive degree days (days with a temperature > 0 degrees Celsius) by scaling positive temperatures with a degree-day factor (uncertain parameter) which parametrizes all the unresolved processes related to snowmelt and sublimation. 

The idea in a synthetic experiment is to:
1. Perform a so-called "truth run" by setting some true value for the parameters in the respective grid cells. For the rest of the experiment (other than the validation), we pretend that we don't know the true values of the two parameters (degree day factor, snowfall multiplier) and states (SWE) other than in the final validation. 
2. Generate synthetic observations by perturbing the true SWE for a subset of days in the truth runs to mimic reality with sparse and imperfect observations.
3. Set a prior distribution on the uncertain parameters (remember, we wouldn't know the truth in practice) and run an ensemble of degree-day models for the entire water year with different values for the two parameters. 
4. Assimialte the synthetic observations at the end of the water year to update the parameters. Then rerun the water year with the posterior (updated) parameters to get the posterior SWE.
5. Perform some validation, i.e. compare the prior and posterior to the truth, to see if the assimilation is working as intended. Potentially compare the effect of using different: (i) Assimilation schemes, (ii) Observation errors, (iii) Observation density, (iv) Observation types, (v) Snow models, (vi) Localization routines.

## Contents
- main.m : Runs the synthetic experiment.
- DDM.m : The degree-day snowmelt model.
- EnKA.m : Analysis (update) step for the ensemble Kalman-based data assimilation schemes. 
- TopoSCALE.mat : Downscaled meteorological forcing for the two grid cells. 


## Currently implemented DA schemes
- Ensemble Kalman-based schemes
	- [x] Ensemble smoother
	- [x] Iterative ensemble smoother (ES-MDA)

## TO DO

- Add particle-based schemes and MCMC.
- Investigate different approaches to propagate information in space so as to update snowpack states in unobserved locations.

