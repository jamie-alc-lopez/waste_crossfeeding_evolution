%This script generates the default parameter structure used

clear;clc

par.Keq = [1,1]; %Keq values
par.K = [1,1]; %Monod half-saturation constants if used
par.m = 0; %Number of microbes
par.n_rxn = 2; %Number of reactions
par.E = []; %Metabolic enzymes
par.T = []; %Transport enzymes
par.beta = 2e2; %Transport coefficient
par.alpha = [1.8e-1,1.8e-1]; %Reaction yield values
par.nu = 1; %Metabolic enzyme coefficient
par.supply = 1; %Nutrient supply rate
par.D = 0.1; %Dilution rate
par.max_internal_time = inf; %Max amount of simulation time permitted
par.max_runtime = 5; %Max amount of real time permitted 
par.exceeded_time = 0; %Tracking variable for whether time is exceeded
par.VcVC = 2e-15; %Ratio of cell to chemostat volume
par.min_growth = 5e-4; %Minimum growth rate
par.tol = 1e-13; %Integrator tolerance
par.pseudo_tol = 1e-12; %dy tolerance for stopping pseudo SS integrations
par.growth_tol = 1e-3; %Relative tolerance for invader growth improvement
par.extol = 1e2; %Threshold for determining species extinction
par.neglect_pseudo_growth_dilution = 1; %Whether to neglect growth dilution in finding optima
par.plt = 0; %Whether to plot individual epoch
par.plt_invasion = 0; %Whether to plot entire evo trajectory
par.save = 0; %whether to save data from individual epoch
par.save_invasion = 0; %Whether to save data from entire evo trajectory
par.invasion_biomass = 1e11; %Invasion biomass for new strain classes
par.n_rounds = 100; %Max number of evolution epochs
par.max_pseudo_rounds = 100; %Max rounds of simulation for calculation invader pseudo-equi
par.osm_coeff = 0; %Osmotic toxicity coefficient
par.no_saturation = true; %Whether to remove saturation
par.no_thermo = false; %Whether to remove thermodynamic toxicity
par.fast_pseudo = 0; %Whether to use fast analytical pseudo-equi approximation
par.no_dilution = 0; %Whether to exclude dilution
par.only_P1P2 = 0; %Similation only P1P2 consortia 
par.existing_environment = 0; %Whether to start from existing environment
par.reporting = 'off'; %Reporting of fmincon
par.plt_opt_results = 0; %Whether to plot optimization outcome

save('default_parameter_struct.mat');