classdef time_integ < BaseClass
    properties (Constant)
        t_interval = 1e3;
        thr = 0.01;     % maximum allowed steady-state biomass deviation 
                        % to be within thr * empirical biomass
    end
    
    properties
        base;
        eventfnc;
        options2;
        count = 0;
        outputfname;
        summary = struct([]);
        t_out;
        b_out;
        b_steady; 
        b_init;
        terminal_event = '';
        neval_points = 0;
        batchmode = false;
        singular = false;
        extinction_flags;
        rawcaseType;
        caseType;
        stableres;
    end
    
    methods
        function obj = time_integ(base, outputfname, i, name)            
            obj.base = base;
            str = sprintf('_%03d', i);
            obj.outputfname = strcat(outputfname, str, '_', name);
        end
        
        function obj = save(obj)
            folder = split(obj.outputfname, '/');
            sub_output_dir = strcat(strjoin(folder(1:end-1),'/'), '/');
            dir_name = strcat(obj.base.output_dir, sub_output_dir);
            if not(isfolder(dir_name))
                mkdir(dir_name);
            end
            
            if ~isempty(obj.summary)
                filename = strcat(obj.base.output_dir, obj.outputfname, '_summary');
                T = struct2table(obj.summary);
                writetable(T, [filename, '.csv']);
            else
                time_integ_struc.base = obj.base;
                time_integ_struc.t_out = obj.t_out;
                time_integ_struc.b_out = obj.b_out;
                time_integ_struc.b_steady = obj.b_steady;
                time_integ_struc.b_init = obj.b_init;
                time_integ_struc.rawcaseType = obj.rawcaseType;
                filename = strcat(obj.base.output_dir, obj.outputfname, '_', datestr(now,'yyyymmdd'));
                obj.outputfname = filename;
                save(strcat(filename, '.mat'), 'time_integ_struc', '-v7.3')
            end

            settings_fpath = strcat(obj.base.output_dir, sub_output_dir, 'settings.mat');
            if ~exist(settings_fpath)
                settings = {};
                settings.thr = time_integ.thr;
                settings.parameters = obj.base;
                save(settings_fpath, 'settings', '-v7.3');
            end
        end
        
        function obj = update_summary(obj, netInfo, bfixedpoint)
            idx = length(obj.summary)+1;
            obj.summary(idx).netID = netInfo.netID;
            obj.summary(idx).netName = netInfo.name;
            obj.summary(idx).sd = netInfo.sd;
            obj.summary(idx).sr = netInfo.sr;
            J = BaseClass.compute_jacobian_triple(netInfo, bfixedpoint, netInfo.b0);
            [domeigval,~] = BaseClass.get_domeigcomp(J);
            clear J
            
            obj.summary(idx).domeigreal = real(domeigval);
            obj.summary(idx).domeigimag = imag(domeigval);
            db = BaseClass.compute_db(netInfo, netInfo.b0);
            obj.summary(idx).normdb0 = norm(db);
            obj.summary(idx).terminal_event = obj.terminal_event;
            obj.summary(idx).neval_points = obj.neval_points;
        end
        
        function obj = get_stability_heuristic(obj)
            nComp = size(obj.b_out,2);
            unstable_comp = [];
            obj.stableres = NaN(1,nComp);
            obj.rawcaseType = zeros(1,3);
            idx = find(obj.t_out >= (obj.base.t_max - time_integ.t_interval));
            if isempty(idx)
                start_idx = 1;
            else
                start_idx = idx(1);
            end
            
            % Loop through each compartment and classify steady-state
            % behavior
            for c=1:nComp
                b_rel = obj.b_out(start_idx:end,c) ./ obj.b_steady(c);
                res = abs(1 - b_rel);

                % Fit straight line to end of cycle and extract slope
                coef = polyfit(obj.t_out(start_idx:end),obj.b_out(start_idx:end, c),1);
                slope = coef(1); 
                b_out_avg = mean(obj.b_out(start_idx:end, c));
                isconverging = (b_out_avg >= obj.b_steady(c) && slope < 0) || (b_out_avg <= obj.b_steady(c) && slope > 0);
                fprintf('\nCompartment: %i, slope is: %.6f, converging: %i', c, slope, isconverging);
   
                if abs(slope) <= obj.base.slope_ub || isconverging
                    chk_bout_emp = res <= time_integ.thr;
                    if all(chk_bout_emp)
                        obj.rawcaseType(1) = obj.rawcaseType(1) + 1;
                        obj.stableres(c) = res(end);
                    else
                        obj.rawcaseType(2) = obj.rawcaseType(2) + 1;
                        obj.stableres(c) = res(end);
                    end
                else
                    obj.rawcaseType(3) = obj.rawcaseType(3) + 1;
                    unstable_comp = [unstable_comp c];
                end
            end

            if any(obj.stableres > 1)
                fprintf('\nCheck large residuals: %s', mat2str(obj.stableres))
                fprintf('\nCorresponding empirical biomass: %s', mat2str(obj.b_steady(~isnan(obj.stableres))'))
                fprintf('\nCorresponding time integrated biomass: %s', mat2str(obj.b_out(end, ~isnan(obj.stableres))'))
            end
            fprintf('\nUnstable compartments: %s', mat2str(unstable_comp));
            obj.caseType = round(obj.rawcaseType / nComp, obj.base.precision);
        end

        function obj = check_consistency(obj)
            errmsg = '';
            obj = obj.get_stability_heuristic;

            % Update summary
            idx = length(obj.summary);
            obj.summary(idx).maxres = nanmax(obj.stableres);
            obj.summary(idx).minres = nanmin(obj.stableres);
            obj.summary(idx).moderes = mode(obj.stableres);
            obj.summary(idx).case1 = obj.caseType(1);
            obj.summary(idx).case2 = obj.caseType(2);
            obj.summary(idx).case3 = obj.caseType(3);
            M = max(obj.caseType);
            nMax = length(find(M==obj.caseType));
            
            if nMax > 1
                I = -1;
            elseif obj.caseType(1) == 1
                I = 1;
            elseif obj.caseType(1) + obj.caseType(2) == 1
                I = 2;
            else
                I = 3;
            end
            
            switch I
                case 1
                    if ~(obj.summary(idx).domeigreal < 0)
                        errmsg = 'Failed convergence to empirical biomass';
                    end
                case 2
                    if ~(obj.summary(idx).domeigreal < 0)
                        errmsg = 'Failed convergence to a stable biomass';
                    end
                case 3
                    if obj.summary(idx).domeigreal < 0
                        errmsg = 'Failed convergence, unstable dynamics';
                    end
                case -1
                    errmsg = 'Unclear network stability';
            end
            
            obj.summary(idx).error = errmsg;
            obj.summary(idx).singular_err = obj.singular;
        end
                        
        function obj = simulate_dynamics(obj, netInfo, b_init)            
            obj = obj.solve_ODE(netInfo, b_init);
            obj.b_init = b_init;
            obj.b_steady = netInfo.b0;
            if any(obj.b_steady < 0)
                warning('Negative steady-state solution');
                obj.b_steady
                netInfo.b0
                netInfo.p
            end
        end

        function obj = solve_ODE(obj, netInfo, b_init)  
            global init_flag
            init_flag = 1;
            obj.base.set_random_seed;
            t_span = [0 obj.base.t_max];
            obj.terminal_event = '';
            lastwarn('', '');
            
            if isnan(obj.base.event_termination) || obj.base.event_termination
                obj.eventfnc = @(t, b) obj.OOMEvent(t, b, netInfo.b0);
                obj.options2 = odeset(obj.base.options, 'Events', obj.eventfnc);
                if strcmp(obj.base.solver,'ode45') 
                    [t_out, b_out, t_out_e, b_out_e, i_e] = ode45(@(t,b) obj.differential(t, b, netInfo) , t_span, b_init, obj.options2);
                elseif strcmp(obj.base.solver, 'ode23s')
                    [t_out, b_out, t_out_e, b_out_e, i_e] = ode23s(@(t,b) obj.differential(t, b, netInfo) , t_span, b_init, obj.options2);
                end
            elseif ~obj.base.event_termination
                if strcmp(obj.base.solver,'ode45') 
                    [t_out, b_out] = ode45(@(t,b) obj.differential(t, b, netInfo) , t_span, b_init, obj.base.options);
                elseif strcmp(obj.base.solver, 'ode23s')
                    [t_out, b_out] = ode23s(@(t,b) obj.differential(t, b, netInfo) , t_span, b_init, obj.base.options);
                end
            end

            [warnMsg, warnId] = lastwarn();
            
            if(~isempty(warnId))
                disp([warnMsg, warnId]);
                obj.singular = true;
            end

            obj.b_out = b_out; obj.t_out = t_out;
            obj.neval_points = size(obj.b_out,1);
            b = obj.b_out(end,:)';
            obj.extinction_flags = (b <= netInfo.b0 * obj.base.Delta);
            fprintf('\nPost-termination check...')
            %disp(extinction_flags)
            if any(obj.extinction_flags)
                fprintf('\nCurrent b_out values: %s\n', mat2str(b))
                terminal_node = find(obj.extinction_flags); living = terminal_node <= netInfo.l;
                if strcmp(obj.terminal_event,'')
                    obj.terminal_event = sprintf('extinction at node %i of %i', terminal_node, netInfo.n);
                end

                if ~living
                    obj.terminal_event = sprintf('depletion at node %i of %i where l = %i', terminal_node, netInfo.n, netInfo.l);
                end
            end
            
            if ~strcmp(obj.terminal_event, '')
                fprintf('Terminating integration because: %s\n', obj.terminal_event);
            end
        end
        
        function db = differential(obj, t, b, netInfo)
            db = BaseClass.compute_db(netInfo, b);            
        end
        
        function obj = eigenmode_expansion(obj, netInfo, x0)
            T = 0:0.01:obj.base.t_max;
            xt = NaN(netInfo.n,length(T));
            J = BaseClass.compute_jacobian_triple(netInfo, netInfo.b0, netInfo.b0);
            [V,D] = eig(J);
            eigenvalues = diag(D);
            cnst = V \ x0;

            for t=1:length(T)
                xt(:,t) = V * (exp(eigenvalues * T(t)) .* cnst);
                [Mb, Ib] = max(xt(:,t));
                if Mb > 1e4
                    break;
                end
            end

            % Remove timepoints with NaN
            tend = find(isnan(xt(1,:)));
            if ~isempty(tend)
                xt = xt(:,1:tend-1);
                T = T(1:tend-1);
            end
            
            obj.t_out = T';
            obj.b_out = xt' + repmat(netInfo.b0', length(obj.t_out), 1);
            obj.b_init = x0 + netInfo.b0;
            obj.b_steady = netInfo.b0;
        end 

        function [value,isterminal,direction] = OOMEvent(obj, t, b, b0)            
            isterminal = 1;              % Halt integration 
            direction = 0;               % The zero can be approached from either direction
            
            extinction_flags = (b <= b0 * obj.base.Delta);
            value1 = int8(any(extinction_flags));  % integration terminates when at least 1 extinction has occurred
            value = 1-value1;  % event is terminal when value is 0
        end
        
        function nres = analyse_ODE_solution(obj, text, name, foodwebID, compID, plot_settings, plotting, printing)
            if plotting
                fprintf('\nBegin plotting...');
                plotobj = PlotClass(obj.base, 'biomass timecourse', plot_settings, obj.outputfname, true, name, foodwebID, compID);
            end
            if printing
                fprintf('\n%s', text);
                fprintf('\nb_steady:')
                fprintf('\n%s', mat2str(obj.b_steady.'));
                fprintf('\nb_out final:')
            end
            
            b_out_final = obj.b_out(end, :).';
            res = b_out_final - max(obj.b_steady, 0);
            nres = norm(res);
            if printing
                if nres <= 1e-2
                    disp('the same');
                else
                    disp(b_out_final.');
                    disp('res:');
                    disp(res.');
                    warning(['DIFFERENT --- DID NOT CONVERGE TO STEADY STATE, norm = ' num2str(norm(res))]);
                end
            end
        end        
    end
end