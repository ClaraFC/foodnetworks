clc; clear all; close all; 
SETUP_PATH = 'setup.m';
run(SETUP_PATH)

%% Perform time integration to observe biomass dynamics of individual food web
parameters.solver = 'ode23s';
parameters.t_max = 1e2; 
parameters.t_max_plot = parameters.t_max;
parameters.slope_ub = 0.01;
parameters.linearY = true;
parameters.sd = 0.1; parameters.sr = 0;
parameters.array_sd = 0:0.1:1;
parameters.array_sr = 0:0.1:1;
parameters.binitopt.integ.mode = 'lognormal';
parameters.binitopt.integ.cv = 0.01;
parameters.event_termination = true;
disp(parameters)

parameters.exclude_networks_list=false;

parameters.event_termination = true;
disp(parameters)
fname = sprintf('out/', ...
    parameters.solver, parameters.sd, parameters.sr, parameters.binitopt.integ.cv);
process = 'simulate biomass dynamics';
target_name = NaN;

foodwebs_dir = '/Users/cseunitclara/Documents/OIST/dir/'

names=dir('/Users/cseunitclara/Documents/OIST/dir/*.m')

for i=1:length(names);
    c{i}=[ foodwebs_dir,names(i).name];
end

foodweblist = c; compID = NaN; parmode = true; parprocesses = 7;

outputfname1 = dynamic_model_triple(parameters, process, plot_settings, fname, target_name, foodweblist, compID, parmode, parprocesses);


