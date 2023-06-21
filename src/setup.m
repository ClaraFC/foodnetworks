parameters = {};
parameters.exclude_networks_list = []; % IDs of networks to exclude from analysis
parameters.array_sd = NaN;  % range of donor-control strengths. 
                            % Log scale, 10.^[-6:0.1:0] or regular scale 0:0.01:1
parameters.array_sr = NaN;  % range of recipient-control strengths 
                            % (by default will be set to array_sd values)
parameters.sd = NaN;        % Donor control strength. Any value from [0,1] 
parameters.sr = NaN;        % Recipient control strength. Any value from [0,1], where sd + sr <= 1
parameters.binitopt.integ.mode = '';     % Specify initial biomass values for time integration
                                         % Options: 'normrand', randomly sample from normal distribution. Adjust variance in BaseClass.m -> get_initbiomass function
                                         %          'percentrand', randomly sample some +/- percentage of empirical biomass, b0
                                         %          'lognormal', lognormal distribution with mean at the empirical biomass, b0, and coefficient of variation cv  
parameters.binitopt.integ.percent = NaN; % If using 'percentrand' for mode, set percent between [0, 100] eg. 0.001
parameters.bibitopt.integ.cv = NaN;      % If using 'lognormal' for mode, select the coefficient of variation
parameters.binitopt.jacobian.mode = '';  % Specify about which biomass to calculate Jacobian
parameters.binitopt.jacobian.percent = NaN;
parameters.event_termination = true;  % Set to true if integration should 
                                     % terminate when negative biomass 
                                     % values occur. Otherwise set to false.
parameters.solver = '';      % name of ode solver (eg. ode45, ode23s) 
                             % to be used in analyses/time_integ.m ->
                             % solve_ODE function
parameters.t_max = NaN;      % maximum timepoint to run integration up to
parameters.slope_ub = NaN;   % Maximum slope magnitude of trajectory which qualifies as a straight line
parameters.t_max_plot = NaN; % maximum timepoint to plot time integrated 
                             % biomass trajectories in PlotClass.m 
                             % -> plot_biomass_timecourse function 
parameters.logX = false;     % true or false to plot x-axis on log scale
parameters.logY = false;
parameters.linearX = false;  % true or false to plot x-axis on linear scale
parameters.linearY = false;
parameters.local = true;
parameters.cluster = false;  % true or false if running script on cluster
parameters.compute_locally_ko = false; % sets local paths

plot_settings = {};
plot_settings.fname = '';     % filename to save plots with
                              % default: generic filename according to plot type
plot_settings.colormap = NaN; % color scheme for triangle heatmap plots
                              % default: jet
plot_settings.title = '';     % specify plot title
                              % default: generic title according to plot type
plot_settings.drawplot = true;
plot_settings.newfigure = false;
plot_settings.exclude_networks = false;
plot_settings.trace_eigenmodes = false;  % when running individual time integration,
                                         % set to true to overlay eigenmode expansion trajectories

if parameters.local
    SRC_DIR = '/Users/cseunitclara/Documents/OIST/foodweb-dynamics-rotationproj-master/src/';
elseif parameters.cluster
    SRC_DIR = getenv('SRC_DIR');
else
    SRC_DIR = '/bucket/DieckmannU/Clara/foodweb-dynamics-rotationproj/src/';
end

addpath(genpath(SRC_DIR));