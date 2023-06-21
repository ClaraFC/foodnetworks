classdef FPstability < BaseClass
    properties
        base
        stability_mat
        domeig_mat
        eig_mat
        valid_fixedpoint_mat
        outputfname
    end
    
    methods (Static)
        function [stability, domeigval] = check_stability(J)
            [domeigval,~] = BaseClass.get_domeigcomp(J);

            % checking stability criteria for continuous time
            if real(domeigval) < 0
                stability = 1;
            else
                stability = 0;
            end
        end
    end

    methods
        function obj = FPstability(base, file_number, outputfname)
            % Initializing class instance
            obj.base = base;
            obj.stability_mat = NaN(base.N_sd, base.N_sr, file_number);
            obj.domeig_mat = NaN(base.N_sd, base.N_sr, file_number);
            obj.eig_mat = NaN(base.N_sd, base.N_sr, file_number, base.nCompartmentsMAX);
            obj.outputfname = outputfname;
        end
       
        function grid_locations = get_stability_grid(obj)
            grid_locations = NaN(obj.base.N_sd, obj.base.N_sr, 2);
            for dID = 1:obj.base.N_sd
                sd = round(obj.base.array_sd(dID), obj.base.precision);
                for rID = 1:obj.base.N_sr
                    sr = round(obj.base.array_sr(rID), obj.base.precision);
                    if round(sd+sr, obj.base.precision) <= 1
                        grid_locations(dID, rID, 1) = sd;
                        grid_locations(dID, rID, 2) = sr;
                    end
                end
            end
        end

        function obj = run_analysis(obj, netInfo, b_fixedpoint, b_balancing)
            if ~isnan(netInfo.sd) & ~isnan(netInfo.sr)             
                J = BaseClass.compute_jacobian_triple(netInfo, b_fixedpoint, b_balancing);
                [stability, domeigval] = FPstability.check_stability(J);
                fprintf('netID: %i, re. J: %.3f\n', netInfo.netID, real(domeigval));
                if real(domeigval) >= 0 
                    disp('Check!')
                end
            else
                % Here we envisage two for loops over a simplex of (sd, sr, sl)
                for dID = 1:obj.base.N_sd
                    netInfo.sd = round(obj.base.array_sd(dID), obj.base.precision);
                    for rID = 1:obj.base.N_sr
                        netInfo.sr = round(obj.base.array_sr(rID), obj.base.precision);
                        if round(netInfo.sd + netInfo.sr, obj.base.precision) <= 1
                            netInfo.sl = round(1 - netInfo.sd - netInfo.sr, obj.base.precision);
                            J = BaseClass.compute_jacobian_triple(netInfo, b_fixedpoint, b_balancing);
                            [obj.stability_mat(dID,rID,netInfo.netID), obj.domeig_mat(dID,rID,netInfo.netID)] = ...
                                FPstability.check_stability(J);
                            obj.eig_mat(dID, rID, netInfo.netID, 1:netInfo.n) = eig(J);
                        end
                    end
                end
            end
        end
                
        function filepath = save(obj)
            folder = split(obj.outputfname, '/');
            dir_new = strcat(obj.base.output_dir, strjoin(folder(1:end-1),'/'));
            if ~exist(dir_new, 'dir')
                mkdir(dir_new);
            end
            
            stability_struc.data = obj.stability_mat;
            stability_struc.domeig = obj.domeig_mat;
            stability_struc.eig = obj.eig_mat;
            stability_struc.sd = obj.base.array_sd;
            stability_struc.sr = obj.base.array_sr;
            filepath = [obj.base.output_dir, '/' obj.outputfname, '.mat'];
            save(filepath, 'stability_struc');
            disp(['Saved results as ' filepath]);
        end
    end
end