Scripts and data sets for Gerhardt et al 2021

fig1ddata 

> experimental data: inducer concentration, mean, CV

directMethod.m
firstReactionMethod.m

> files for running Gillespie SSA (add source- mathworks file exchange)
> Used only for LuxR-AHL circuit stochastic model. TetR-aTc modeled in SimBiology

>> LuxR-AHL (high noise) circuit
luxR_simplesteadystate.m, luxsimpleode.m --> Deterministic model 

luxrpropensities --> Stochastic model; propensity function for Gillespie SSA
luxrsim_exn --> Gillespie simulation script with extrinsic noise

