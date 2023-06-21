function varargout = dynamic_model_triple(parameters, process, plot_settings, varargin)
    % ~~ Input arguments ~~
    % parameters: defined in sample_main_script.m
    % process: analyses to carry out. Options: 'analyze fixed points' or 'simulate biomass dynamics'
    % 
    % ~~ Additional input arguments ~~
    % when process is 'simulate biomass dynamics'
    % outputfname: main file name to store 3 outputs 
    %              .mat file with time and biomass values of trajectories,
    %              .png files with log or linear axes as specified
    % target_name: name of food web to plot timecourse for eg. 'n007-BlackSea_433'
    % foodweblist: [only for batch integration] list of food web IDs
    % compID:      ID of a specific compartment to plot the timecourse of
    % parmode:     [only for batch integration] 
    %              true of false if using parallel processing in matlab
    % parprocesses:[only for batch integration] number of workers to run
    %               if parmode is true

    %warning on backtrace
    base = BaseClass(parameters);
    names = base.get_names(base);
    filepaths = get_input_filepaths(base.input_dir);
    file_number = length(filepaths);
    
    switch process
        case 'analyze fixed points'
            outputfname = varargin{:};
            FPstability_obj = FPstability(base, file_number, outputfname);
        case 'simulate biomass dynamics'
            [outputfname, target_name, foodweblist, compID, parmode, parprocesses] = varargin{:};
        case 'eigenmode expansion'
            [target_name, x0] = varargin{:};
    end
    
    switch process
        case 'analyze fixed points'
                for i = 1:file_number
                    filepath = filepaths{i};

                    [fname, crc, n, l, b0, p, q, r, F] = read_input_file(filepath);
                    netInfo = pack_netInfo(i, names{i}{:}, b0, F, p, q, r, n, l, NaN, NaN, NaN);

                    disp([num2str(i), '/', num2str(length(filepaths)), ' (', ...
                       num2str(round(100*i/file_number)), '%) ', filepath ', ', ...
                       fname, ', n = ', num2str(n)])

                    b_init = base.get_initbiomass(netInfo, base.binitopt.jacobian);
                    FPstability_obj = FPstability_obj.run_analysis(netInfo, b_init, b_init);
                end
                varargout{1} = FPstability_obj.save();

                if plot_settings.drawplot
                    plot = PlotClass(base,'control space stability', plot_settings, outputfname, true, plot_settings);
                end
        case 'simulate biomass dynamics'
            if ~isnan(target_name)
                i = find(strcmp([names{:}], target_name));
                filepath = filepaths{i};
            
                [fname, crc, n, l, b0, p, q, r, F] = read_input_file(filepath);
                netInfo = pack_netInfo(i, target_name, b0, F, p, q, r, n, l, base.sd, base.sr, base.sl);

                disp([num2str(i), '/', num2str(length(filepaths)), ' (', ...
                num2str(round(100*i/file_number)), '%) ', filepath ', ', ...
                fname, ', n = ', num2str(n)])

                time_integ_obj = time_integ(base, outputfname, i, target_name);
                b_init = base.get_initbiomass(netInfo, base.binitopt.integ);
                fprintf('\nEmpirical biomass values: %s', mat2str(b0'));
                fprintf('\nInitial biomass values: %s', mat2str(b_init'));
                time_integ_obj = time_integ_obj.simulate_dynamics(netInfo, b_init);
                time_integ_obj = time_integ_obj.get_stability_heuristic;
                time_integ_obj = time_integ_obj.save;
                time_integ_obj.analyse_ODE_solution('Stable state simulation', target_name, i, compID, ...
                plot_settings, base.plotting, base.printing);
                varargout{1} = time_integ_obj.outputfname;
            else
                disp(foodweblist)
                if parmode
                    poolobj = parpool(parprocesses);
                    parfor w = 1:length(foodweblist)
                        i = foodweblist(w);
                        run_batch_time_integration(w, base, filepaths, names, outputfname);
                    end
                    delete(poolobj)
                else
                    for w = 1:length(foodweblist)
                        %i = foodweblist(w);
                        run_batch_time_integration(w, base, filepaths, names, outputfname);
                    end
                end
                varargout={NaN};
            end
        case 'profile biological system'
            for i = 1:file_number
                filepath = filepaths{i};
        
                [fname, crc, n, l, b0, p, q, r, F] = read_input_file(filepath);
                netInfo = pack_netInfo(i, names{i}{:}, b0, F, p, q, r, n, l, base.sd, base.sr, base.sl);
        
                disp([num2str(i), '/', num2str(length(filepaths)), ' (', ...
                   num2str(round(100*i/file_number)), '%) ', filepath ', ', ...
                   fname, ', n = ', num2str(n)])

                bfixedpoint = netInfo.b0;
                J = BaseClass.compute_jacobian_triple(netInfo, bfixedpoint, netInfo.b0);
                D = eig(J);

                % Statistic 1: Fraction of positive real parts
                posCountFrac(i) = sum(real(D) > 0) / netInfo.n;

                % Statistic 2: Relative magnitude of positive real parts
                posMagFrac(i) = sum(D(real(D) > 0)) / sum(abs(D));
            end

            % Plot histograms
            subplot(1,2,1);
            histogram(posCountFrac);
            titlestr = sprintf('Fraction of eigenvalues with positive real parts, \ns_d = %.2f and s_r = %.2f', netInfo.sd, netInfo.sr);
            title(titlestr);

            subplot(1,2,2);
            histogram(posMagFrac);
            titlestr = sprintf('Fraction of positive real part magnitude, \ns_d = %.2f and s_r = %.2f', netInfo.sd, netInfo.sr);
            title(titlestr);
        case 'eigenmode expansion'
            netID = find(strcmp([names{:}], target_name));
            filepath = filepaths{netID};
        
            [fname, crc, n, l, b0, p, q, r, F] = read_input_file(filepath);
            netInfo = pack_netInfo(netID, names{netID}{:}, b0, F, p, q, r, n, l, base.sd, base.sr, base.sl);

            disp([num2str(netID), '/', num2str(length(filepaths)), ' (', ...
            num2str(round(100*netID/file_number)), '%) ', filepath ', ', ...
            fname, ', n = ', num2str(n)])

            % Run eigenmode expansion
            time_integ_obj = time_integ(base, NaN, netInfo.netID, netInfo.name);
            time_integ_obj = time_integ_obj.eigenmode_expansion(netInfo, x0);
            time_integ_obj = time_integ_obj.get_stability_heuristic();
            outputfolder = [base.output_dir, 'time_integ/indv/eigenmode_expansion/',target_name,'/'];
            mkdir(outputfolder);
            outputfname = [outputfolder,netInfo.name];
            time_integ_struc.t_out = time_integ_obj.t_out;
            time_integ_struc.b_out = time_integ_obj.b_out;
            time_integ_struc.b_init = time_integ_obj.b_init;
            time_integ_struc.b_steady = time_integ_obj.b_steady;
            time_integ_struc.rawcaseType = time_integ_obj.rawcaseType;
            save(outputfname, 'time_integ_struc');
            compID = NaN;
            plotobj = PlotClass(base, 'biomass timecourse', plot_settings, outputfname, true, netInfo.name, netInfo.netID, compID);
            varargout{1} = outputfname;
    end    
end

function run_batch_time_integration(i, base, filepaths, names, outputfname)
    file_number = length(filepaths);
    filepath = filepaths{i};

    [name, crc, n, l, b0, p, q, r, F] = read_input_file(filepath);

    disp([num2str(i), '/', num2str(length(filepaths)), ' (', ...
        num2str(round(100*i/file_number)), '%) ', filepath ', ', ...
        name, ', n = ', num2str(n)])

    time_integ_obj = time_integ(base, outputfname, i, names{i}{:});
    for dID = 1:base.N_sd
        sd = round(base.array_sd(dID), base.precision);
        for rID = 1:base.N_sr
            sr = round(base.array_sr(rID), base.precision);
            if round(sd+sr, base.precision) <= 1
                sl = round(1 - sd - sr, base.precision);
                netInfo = pack_netInfo(i, names{i}{:}, b0, F, p, q, r, n, l, sd, sr, sl);
                b_init = base.get_initbiomass(netInfo, base.binitopt.integ);
                bfixedpoint = base.get_initbiomass(netInfo, base.binitopt.integ);
                fprintf('\nFood web: %s, ID: %d, Control space: sd = %.3f, sr = %.3f\n', name, i, sd, sr)
                fprintf('\nEmpirical biomass values: %s', mat2str(b0'));
                fprintf('Initial biomass values: %s\n', mat2str(b_init))
                time_integ_obj = time_integ_obj.simulate_dynamics(netInfo, b_init);
                time_integ_obj = time_integ_obj.update_summary(netInfo, bfixedpoint);
                time_integ_obj = time_integ_obj.check_consistency();
            end
        end
    end
    time_integ_obj = time_integ_obj.save;
end

function filepaths = get_input_filepaths(input_dir)
    files = dir(input_dir);
    filepaths = {};
    for i = 1:length(files)
        if endsWith(files(i).name, '.m') && ~startsWith(files(i).name, '._')
            filepaths{length(filepaths)+1} = fullfile(input_dir, files(i).name);
        end
    end
end

function [name, crc, n, l, b0, p, q, r, F] = read_input_file(filepath)
    run(filepath);
    if ~exist('title', 'var')
        title = '';
    end
    if ~exist('crc', 'var')
        crc = NaN;
    end
    name = food_web_filename;
end

function netInfo = pack_netInfo(netID, name, b0, F, p, q, r, n, l, sd, sr, sl)
    netInfo.netID = netID;
    netInfo.name = name;
    netInfo.b0 = b0;
    netInfo.F = F;
    netInfo.p = p;
    netInfo.q = q;
    netInfo.r = r;
    netInfo.n = n;
    netInfo.l = l;
    netInfo.sd = sd;
    netInfo.sr = sr;
    netInfo.sl = sl;
end