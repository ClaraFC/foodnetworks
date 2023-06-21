clc; clear all; close all; 
SETUP_PATH = 'setup.m';
run(SETUP_PATH)

%% Perform time integration to observe biomass dynamics of individual food web
parameters.solver = 'ode23s';
parameters.t_max = 1e2; 
parameters.t_max_plot = parameters.t_max;
parameters.slope_ub = 0.01;
parameters.linearY = true;
parameters.sd = NaN; 
parameters.sr = NaN;
parameters.array_sd = 0:0.1:1;
parameters.binitopt.integ.mode = 'lognormal';
parameters.binitopt.integ.cv = 0.01;
parameters.binitopt.jacobian.mode = 'empirical';
parameters.event_termination = true;
disp(parameters)

fname = sprintf('time_integ/indv/timeInteg_%s_sd%.2f_sr%.2f_percent%.3f', ...
    parameters.solver, parameters.sd, parameters.sr, parameters.binitopt.integ.cv);
fname = 'jacobian_20230609'
process = 'analyze fixed points' %'simulate biomass dynamics';
target_name = 'Bothnian_Bay_Sandberg_2000';
foodweblist = []; compID = NaN; parmode = false; parprocesses = 0;
outputfname1 = dynamic_model_triple(parameters, process, plot_settings, fname, target_name, foodweblist, compID, parmode, parprocesses);

% Run eigenmode expansion using initial biomass values from time
% integration
plot_settings.trace_eigenmodes = false;
if plot_settings.trace_eigenmodes
    process = 'eigenmode expansion';
    plotdata = load([outputfname1,'.mat']);
    x0 = plotdata.time_integ_struc.b_init - plotdata.time_integ_struc.b_steady;
    outputfname2 = dynamic_model_triple(parameters, process, plot_settings, target_name, x0);
    if parameters.logY
        outputfname1 = [outputfname1,'_log'];
        outputfname2 = [outputfname2,'_log'];
    end
    p1 = openfig(strcat(outputfname1,".fig"));
    p2 = openfig(strcat(outputfname2, ".fig"));
    data3 = copyobj(get(gca(p2),'Children'), gca(p1));
    y_lim = get(gca(p2), 'ylim');
    set(gca(p1), 'ylim', y_lim);
    ax_children = get(gca(p1), 'children');
    leg_strings = get(ax_children, 'displayname');
    nComp = length(plotdata.time_integ_struc.b_init);
    legend(data3(nComp:end), leg_strings(nComp:end));
    saveas(gca(p1), [outputfname2, '_overlayed.fig']);
    saveas(gca(p1), [outputfname2, '_overlayed.png']);
end

%% To create triangular plot
base = BaseClass(parameters);
inputfname = outputfname1;
PlotClass(base,'control space stability', plot_settings, inputfname, true);