classdef PlotClass < BaseClass
    properties (Constant)
        FontSizeX = 12;
    end
    
    properties
        base
        plotname
        plotdata
        plot_settings
        inputdata
        ofolder
        mask
    end
        
    methods (Static)
        function obj = get_mask(obj)
            sdvals = obj.inputdata.stability_struc.sd;
            srvals = obj.inputdata.stability_struc.sr;
            obj.mask = NaN(length(sdvals), length(srvals));
            for sd=1:length(sdvals)
                for sr=1:length(srvals)
                    if round(sdvals(sd)+srvals(sr), BaseClass.precision) <= 1
                        obj.mask(sd,sr) = 1;
                    end
                end
            end
        end
        
        function stability_FP = plot_controlspace_stability(obj)
            stability_FP_struc = obj.inputdata.stability_struc;
            stability_FP = stability_FP_struc.data;
            dims = size(stability_FP);
            stability_FP_avg = round(nanmean(stability_FP, 3), BaseClass.precision);
            stability_FP_avg_min = min(stability_FP_avg,[],'all');
            stability_FP_avg_max = max(stability_FP_avg,[],'all');
            obj = PlotClass.get_mask(obj);
            stability_FP_avg = stability_FP_avg .* obj.mask;
            hm = heatmap(stability_FP_avg);
            % formatting
            hm.NodeChildren(3).YDir='normal';
            bin = 10;
            value = stability_FP_struc.sd;
            pos = 1:1:length(value);
            CustomXLabels = compose(strcat("%.", num2str(BaseClass.precision),"f"), value);
            CustomXLabels(mod(pos-1,bin) ~= 0) = " ";
            hm.XDisplayData = pos;
            hm.XDisplayLabels = CustomXLabels;
            hm.YDisplayData = pos;
            hm.YDisplayLabels = CustomXLabels;
            if isnan(obj.plot_settings.colormap)
                hm.Colormap = jet;
            else
                hm.Colormap = obj.plot_settings.colormap;
            end
            hm.ColorLimits = [min(0, stability_FP_avg_min),max(stability_FP_avg_max,1)];
            xlabel('s_r');
            ylabel('s_d');

%             figure;
%             ternpcolor(sd', sr', stability_FP_avg); ternlabel('DONOR', 'RECIPIENT', 'LOTKA-VOLTERRA');
%             shading interp
            if strcmp(obj.plot_settings.title, '')
                titlestr = {'Stability of biomass dynamics at empirically observed biomasses','shown in control space'};
            else
                titlestr = obj.plot_settings.title;
            end
            
            if isfield(obj.inputdata.stability_struc, 'acc')
                titlestr{end+1} = [num2str(obj.inputdata.stability_struc.acc), '% consistency'];
            end
            
            title(titlestr);
            fprintf('\nSaving figure to: %s', obj.ofolder);
            if strcmp(obj.plot_settings.fname, '')
                plotfname = 'controlspace_stability';
            else
                plotfname = obj.plot_settings.fname;
            end
            
            saveas(hm, [obj.ofolder, plotfname, '.fig']);
            saveas(hm, [obj.ofolder, plotfname, '.png']);
        end 
        
        function plot_consistency_check(obj)
            % Load and extract stability information from directory with
            % consistency check results
            stability_mat = NaN(obj.base.N_sd, obj.base.N_sr, obj.base.ALL_networks);
            grid_locations = obj.base.get_grid(obj.base);
            totalpoints = 0; totalconsistency = 0;

            csvFileList = dir(strcat(obj.ofolder, 'timeInteg*.csv'));
            if obj.plot_settings.exclude_networks
                updated_csvFileList = obj.base.filter_networks(csvFileList);
            else
                updated_csvFileList = csvFileList;
            end
            
            for f=1:length(updated_csvFileList)
                file = updated_csvFileList(f);
                data = readtable(strcat(file.folder, '/', file.name));
                nRows = size(data,1);
                consistency_score = zeros(nRows,1);
                if strcmp(obj.plot_settings.stability_criteria, 'case1')
                    stable_score = data.('case1') == 1;
                    unstable_score = data.('case1') == 0;
                    consistency_score = (stable_score & data.('domeigreal')<0) | (unstable_score & data.('domeigreal')>0); 
                elseif strcmp(obj.plot_settings.stability_criteria, 'case1&2')
                    stablefrac = data.('case1') + data.('case2');
                    stable_score = stablefrac == 1;
                    unstable_score = stablefrac ~= 1;
                    consistency_score = (stable_score & data.('domeigreal')<0) | (unstable_score & data.('domeigreal')>=0);
                end
                totalconsistency = totalconsistency + sum(consistency_score);
                totalpoints = totalpoints + nRows;

                for row=1:nRows
                    netID = data.('netID')(row); sd = data.('sd')(row); sr = data.('sr')(row);
                    [dID, rID] = find(grid_locations(:,:,1) == sd & grid_locations(:,:,2) == sr);
                    if strcmp(obj.plot_settings.stability_criteria, 'case1')
                        stability_mat(dID, rID, netID) = data.('case1')(row) == 1;
                    elseif strcmp(obj.plot_settings.stability_criteria, 'case1&2')
                        stablefrac = data.('case1')(row) + data.('case2')(row);
                        stability_mat(dID, rID, netID) = stablefrac == 1;
                    end
                end
            end
            save(strcat(obj.ofolder, 'timeInteg_stability.mat'), 'stability_mat');
            
            % Plot results
            if obj.base.array_sr == 0
                % Draw linear trend at sr = 0
                sr0_stability = stability_mat(:,1,:);
                sr0_stability = reshape(sr0_stability, length(obj.base.array_sd),obj.base.ALL_networks);

                if obj.plot_settings.newfigure
                    figure();
                end
                plot(obj.base.array_sd, nanmean(sr0_stability, 2), '-o');
            else
                % Draw triangular heat map across simplex space
                tmp = load([obj.ofolder, 'timeInteg_stability.mat']);
                obj.inputdata.stability_struc.data = tmp.stability_mat;
                obj.inputdata.stability_struc.acc = totalconsistency*100 / totalpoints;
                obj.inputdata.stability_struc.sd = obj.base.array_sd;
                obj.inputdata.stability_struc.sr = obj.base.array_sr;
                PlotClass.plot_controlspace_stability(obj);
            end
        end
        
        function plot_sr0_stability_curve(obj)
            % Stability curve at s_r = 0
            if obj.plot_settings.newfigure
                figure();
            end
            xvals = obj.inputdata.stability_struc.sd;
            stability_FP_stableonly = obj.inputdata.stability_struc.data;
            stability_FP_stableonly(isnan(stability_FP_stableonly)) = 0;
            stability_FP_avg_stableonly = round(mean(stability_FP_stableonly, 3), BaseClass.precision);
            yvals = stability_FP_avg_stableonly(:,1)';
            p = plot(xvals, yvals, '.-');
            if ~isnan(obj.plot_settings.color)
                p.Color = obj.plot_settings.color;
            end
            xlabel('s_d')
            ylabel('Local Stability Fraction')
            if ~strcmp(obj.plot_settings.title,'')
                title(obj.plot_settings.title)
            else
                title('Local Stability for s_r = 0')
            end
            set(gca,'DefaultTextFontSize',18);
            ylim([0 1])

            if obj.base.linearX
                xlabel('s_d', 'FontSize', PlotClass.FontSizeX);
                fname = strcat(obj.ofolder, 'local_stability_sr0.png');
                saveas(gca, fname);
            end
            
            if obj.base.logX
                set(gca, 'XScale', 'log');
                xlabel('log s_d', 'FontSize', PlotClass.FontSizeX);
                fname = strcat(obj.ofolder, 'local_stability_sr0_logX.png');
                saveas(gca, fname);
            end

%             ranges = [1, 61, 122, 183, 244];
%             for p=1:4
%                 figure('Position', get(0, 'Screensize'));
%                 xvals = obj.inputdata.stability_struc.sd;
%                 subarray = stability_FP(:,:,ranges(p):ranges(p+1));
%                 dims = size(subarray);
%                 yvals = reshape(subarray(:,1,:), dims(1), dims(3));
%                 for w=1:dims(3)
%                     data_stable = yvals(:,w);
%                     data_unstable = data_stable;
%                     data_stable(data_stable == 0) = NaN;
%                     data_unstable(data_stable == 1) = NaN;
%                     data_unstable(data_unstable == 0) = 1;
%                     scatter(xvals, w .* data_stable', 's', 'MarkerEdgeColor', '#a4c3f5', 'MarkerFaceColor', '#a4c3f5');
%                     hold on;
%                     scatter(xvals, w .* data_unstable', 'sr', 'MarkerFaceColor', 'r');
%                 end
%                 hold off;
%                 legend('stable', 'unstable');
%                 yticks(1:dims(3));
%                 yticklabels(string(names(ranges(p):ranges(p+1))));
%                 ax=gca;
%                 ax.FontSize = 7;
%                 set(gca,'TickDir','out');
%                 set(gca,'TickLabelInterpreter', 'none');                
%                 if obj.base.linearX
%                     xlabel('s_d', 'FontSize', PlotClass.FontSizeX);
%                     fname = strcat(BaseClass.output_dir, 'raster_plot', string(p), '.png');
%                     saveas(ax, fname);
%                 end
%                 
%                 if obj.base.logX
%                     set(gca, 'XScale', 'log');
%                     xlabel('log s_d', 'FontSize', PlotClass.FontSizeX);
%                     fname = strcat(BaseClass.output_dir, 'raster_plot', string(p), '_logX.png');
%                     saveas(ax, fname);
%                 end
%             end
        end
        
        function plot_biomass_timecourse(obj, name, foodwebID, compID, outputfile)
            t_out = obj.inputdata.time_integ_struc.t_out;
            b_out = obj.inputdata.time_integ_struc.b_out;
            b_steady = obj.inputdata.time_integ_struc.b_steady;
            nComp = size(b_out,2);

            if ~isnan(compID)
                b_out = b_out(:,compID);
                b_steady = b_steady(compID);
            end

            if obj.plot_settings.newfigure
                figure();
            end
            tc = obj.base.t_max_plot; % Plots are generated from t=0 to t=tc.
            ind = t_out <= tc;
            ind(1) = true;

            if obj.plot_settings.trace_eigenmodes
                h1 = plot(t_out(ind), b_out(ind, :), '--');
            else
                h1 = plot(t_out(ind), b_out(ind, :), '.-');
            end
            y_lim = get(gca, 'ylim');
            set(gca, 'ylim', y_lim);

            for i = 1:size(b_steady, 1)
                if ~isnan(compID)
                    bsteadyval = b_steady;
                else
                    bsteadyval = b_steady(i);
                end
                hl(i) = line([0; max([t_out(ind).'])], [bsteadyval, bsteadyval]);
                set(hl(i), 'LineStyle', ':', 'Color', 'k');
            end

            if nComp < 16
                legendstr = strsplit([num2str(1:nComp), ' Steady-state']);
                legend([h1', hl(1)], legendstr, 'Location', 'best');
            else
                legend([h1(1), hl(1)], 'Integration', 'Steady-state', 'Location', 'best');
            end
            xlabel('t');
            ylabel('b(t)');

            title({['Biomass dynamics of ', name], ...
                ['at sd = ', num2str(obj.base.sd), ' and sr = ', num2str(obj.base.sr)], ...
                'Stability heuristic [empirical biomass, other biomass, unstable]:', ...
                [mat2str(obj.inputdata.time_integ_struc.rawcaseType)]}, ...
                'interpreter', 'none');
            
            if obj.base.linearY
                drawnow;
                saveas(gca, [outputfile, '.png']);
                saveas(gca, [outputfile, '.fig']);
            end
            
            if obj.base.logY
                set(gca, 'YScale', 'log');
                ylabel('log b(t)');
                y_lim = get(gca, 'ylim');
                set(gca, 'ylim', y_lim);
                drawnow;
                fname = strcat(outputfile, '_log');
                saveas(gca, [fname, '.png']);
                saveas(gca, [fname, '.fig']);
            end
        end
    end
    
    methods
        function obj = PlotClass(base, plotname, plot_settings, inputfname, loadfile, varargin)
            % ~~ Input Arguments ~~
            % base: instance of base class with parameters detailed in
            %       sample_main_script.m
            % plotname: type of plot to make. 
            %           Options: 'control space stability', triangular plot
            %                    for any 3D matrix with 1s and 0s
            %                    'consistency check', triangular plot for
            %                    visualizing time integration-based stability and computing consistency score against Jacobian-based stability
            %                    'sr0 stability curve', linear trend at recipient control strength = 0, 
            %                    for different values of donor and LV control to visualize nonmonotonicity
            %                    'biomass timecourse', timecourse of all
            %                    compartments in a single food web at a specific point in control space
            % inputfname: file path of data to be plotted in a '.mat' file,
            %             leave out '.mat' extension
            % loadfile:   true or false on whether to load file specified
            %             in inputfname
            % 
            % ~~ Additional Input Arguments ~~
            % plot_settings: is detailed in sample_main_script.m 

            obj.base = base;
            obj.plotname = plotname;
            obj.plot_settings = plot_settings;
            inputfnamefull = [ inputfname, '.mat'];
            temp = split(inputfname,'/');
            obj.ofolder = strcat(strjoin(temp(1:end-1), '/'),'/');
            if loadfile
                outputfile = inputfname;
                fprintf('\nLoading input data from, %s', inputfnamefull)
                obj.inputdata = load(inputfnamefull);
            end

            switch obj.plotname
                case 'control space stability'
                    obj.plotdata = PlotClass.plot_controlspace_stability(obj);
                case 'consistency check'
                    PlotClass.plot_consistency_check(obj);
                case 'sr0 stability curve'
                    PlotClass.plot_sr0_stability_curve(obj);
                case 'biomass timecourse'
                    if length(varargin) == 3
                        [name, foodwebID, compID] = varargin{:};
                        PlotClass.plot_biomass_timecourse(obj, name, foodwebID, compID, outputfile);
                    else
                        disp('Error!')
                    end
            end
        end
    end
end

