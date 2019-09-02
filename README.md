# non-euclidean-scr

Some simulation experiments investigating the behaviour of non-Euclidean distance in SCR models.

One set of simulations is for the paper in prep "Understanding snow leopard populations and their spatial ecology through Spatial Capture Recapture Analysis across three sites in South Gobi, Mongolia". To rerun these run *tost_simulation.R* and *tost_sim_analysis.R*. These simulations simulate capture histories assuming the same covariate affects both activity centre density and conductance/movement, and see whether the true parameter estimates can be reliably recovered. 

A second simulation in *does_getting_D_and_euc_right_matter.R* tests the effect of model misspecification. It is less complete 

* Line 10 loads the Tost data- Up to line 64 generates the ACs, either uniformly or as a function of ruggedness (true parameter = 1.5, line 45). 
* Up to line 115 generates capture histories again either uniform or ~ ruggedness (true parameter = 0.3, line 76).
* On line 125 you specify which simulated history you want (currently set so both density and distance are function of ruggedness)
* Up to line 167 fits 4 different models, corresponding to getting the D part right/wrong, and the distance part right/wrong. 

For this capture history model 4 is the correct one, and it recovers the true parameter values well.
