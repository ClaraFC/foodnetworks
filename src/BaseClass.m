

classdef BaseClass
    properties (Constant = true)
        plotting = true;
        printing = true;
        precision = 6;
        Delta = 0; % if node's biomass falls below Delta * b0 we treat it as extinct
        options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'InitialStep',0.001, 'Refine', 200);
        nCompartmentsMAX = 125;
        ALL_networks = 245;
    end
    
    properties (SetAccess = private)
        exclude_networks_list = [];
        N_networks = NaN;
        array_sd = NaN;  % range of donor-control strengths
        array_sr = NaN;  % range of recipient-control strengths
        sd = NaN; sr = NaN; sl = NaN;
        N_sd; N_sr; sd_vals; sr_vals;
        binitopt = NaN;
        solver = '';    % name of ode solver (eg. ode45, ode23s)
        slope_ub = NaN;
        t_max = NaN;
        t_max_plot = NaN;
        linearX = false;   % true or false for adjusting figure scales
        logX = false;
        linearY = false;
        logY = false;
        cluster = NaN;
        input_dir = '';
        output_dir = '';
        event_termination = NaN;
    end
    
    methods (Static)
        % Getting names of food webs
        function nameList = get_names(obj)
            fileList = dir(strcat(obj.input_dir,'*.m'));
            folderNames = {fileList.name}';
            folderNames = folderNames(~startsWith(folderNames ,'.'));
            nameList = cell(1,length(folderNames));
            for f=1:length(folderNames)
                fn = folderNames(f);
                fn = fn{1,1}; % converting cell array element to char
                run(strcat(obj.input_dir,fn));
                name = split(food_web_filename, '.');
                nameList{1,f} = name(1);
            end
        end
        
        function grid = get_grid(obj)
            grid = NaN(obj.N_sd, obj.N_sr, 2);
            for dID = 1:obj.N_sd
                s_d = round(obj.array_sd(dID), obj.precision);
                for rID = 1:obj.N_sr
                    s_r = round(obj.array_sr(rID), obj.precision);
                    if round(s_d+s_r, obj.precision) <= 1
                        grid(dID, rID, 1) = s_d;
                        grid(dID, rID, 2) = s_r;
                    end
                end
            end
        end

        function set_random_seed
            rng('default');
        end

        function [domeigval,domeigvec] = get_domeigcomp(J)
            [V,D] = eig(J);
            [sorted_eigval,idx] = sort(diag(D),'descend','ComparisonMethod','real');
            domeigval = sorted_eigval(1);
            domeigidx = idx(1);
            domeigvec = V(:,domeigidx);
        end

        function [Cd, Cr, Cl] = get_consumption_intensity(F, b0, n)
            Cd = F ./ repmat(b0.', n, 1);
            Cr = F ./ repmat(b0, 1, n);
            Cl = F ./ (b0 * b0.');
        end
        
        function J = compute_jacobian_triple(netInfo, b, b_balancing)
            if ~(netInfo.sd>=0 && netInfo.sd<=1 && netInfo.sr>=0 && netInfo.sr<=1 && netInfo.sl>=0 && netInfo.sl<=1)
                disp([netInfo.sd, netInfo.sr, netInfo.sl])
                error('Incorrect value of parameter sd, sr or, sl')
            end
            
            [Cd, Cr, Cl] = BaseClass.get_consumption_intensity(netInfo.F, b_balancing, netInfo.n);
            J_offdiag = netInfo.sd * Cd - netInfo.sr * Cr' ...
                + netInfo.sl * (Cl .* repmat(b, 1, netInfo.n)) ...
                - netInfo.sl * (Cl' .* repmat(b, 1, netInfo.n));
            J_diag = netInfo.sd * (diag(Cd) - diag(netInfo.sd * sum(Cd,1)))...
                + netInfo.sr * (diag(sum(Cr,2)) - diag(Cr))...
                + netInfo.sl * (diag(sum(repmat(b',netInfo.n,1) .* Cl, 2)) - netInfo.sl * diag(sum(repmat(b,1,netInfo.n) .* Cl, 1)))...
                - diag(netInfo.q./b_balancing + netInfo.r./b_balancing);
            J = J_offdiag - diag(J_offdiag) + J_diag;
        end
        
        function db = compute_db(netInfo, b)

            [Cd, Cr, Cl] = BaseClass.get_consumption_intensity(netInfo.F, netInfo.b0, netInfo.n);
            db = netInfo.sd * (sum(Cd .* repmat(b', netInfo.n, 1), 2) - sum(Cd .* repmat(b', netInfo.n, 1), 1)') ...
                + netInfo.sr * (sum(Cr .* repmat(b,1,netInfo.n),2) - sum(Cr .* repmat(b,1,netInfo.n),1)') ...
                + netInfo.sl * (sum(Cl .* repmat(b,1,netInfo.n) .* repmat(b',netInfo.n,1),2) - sum(Cl .* repmat(b',netInfo.n,1) .* repmat(b,1,netInfo.n),1)') ...
                - (netInfo.q./netInfo.b0 + netInfo.r./netInfo.b0).*b + netInfo.p;

            verifyVectorization = false;
            if verifyVectorization
                db2 = zeros(netInfo.n, 1);
                for i = 1:netInfo.n
                    donor = 0;
                    receipment = 0;
                    lv = 0;
                    rest = - (netInfo.q(i)/netInfo.b0(i) + netInfo.r(i)/netInfo.b0(i))*b(i) + netInfo.p(i);
                    for k = 1:netInfo.n
                        donor = donor + netInfo.sd * Cd(i,k) * b(k) - netInfo.sd * Cd(k,i) * b(i);
                        receipment = receipment + netInfo.sr * Cr(i,k) * b(i) - netInfo.sr * Cr(k,i) * b(k);
                        lv = lv + netInfo.sl * Cl(i,k) * b(i) *b(k) - netInfo.sl * Cl(k,i) * b(i) * b(k);
                    end
                    db2(i) = donor + receipment + lv + rest;
                end
                nd = norm(db-db2);
                if nd > 1e-3
                    warning(['Norm of difference ' num2str(nd)]);
                end
            end
            
            %[Mb, Ib] = max(b);
            %if Mb > (100 * max(netInfo.b0))
            %    disp(['Suspicious growth for ', netInfo.name, ' with ID: ', num2str(netInfo.netID), ', Comp: ', num2str(Ib), ' at sd = ', num2str(netInfo.sd), ' and sr = ', num2str(netInfo.sr)])
            %end
        end
        
        function totSize = GetSize(this, varargin) 
           props = properties(this); 
           totSize = 0; 

           for ii=1:length(props) 
              currentProperty = getfield(this, char(props(ii))); 
              s = whos('currentProperty');
              if length(varargin) == 1 && strcmp(props{ii}, varargin{1})
                  totSize = s.bytes / (1024^2); % in Mb
                  break
              else
                totSize = totSize + (s.bytes / (1024^2)); 
              end
           end       
        end
    end
    
    methods
        % Initializing class instance
        function obj = BaseClass(parameters)
            if exist('parameters', 'var')
                % unpack input parameters
                obj.exclude_networks_list = parameters.exclude_networks_list;
                obj.N_networks = BaseClass.ALL_networks - length(obj.exclude_networks_list);
                obj.array_sd = parameters.array_sd;
                if isnan(parameters.array_sr)
                    obj.array_sr = obj.array_sd;
                else
                    obj.array_sr = parameters.array_sr;
                end
                obj.sd = parameters.sd;
                obj.sr = parameters.sr;
                obj.sl = 1 - parameters.sd - parameters.sr;
                obj.N_sd = length(obj.array_sd);  % number of donor-control points to analyze
                obj.N_sr = length(obj.array_sr);  % number of recipient-control points to analyze
                obj.sd_vals = round(obj.array_sd, obj.precision); 
                obj.sr_vals = round(obj.array_sr, obj.precision);

                obj.event_termination = parameters.event_termination;
                obj.binitopt = parameters.binitopt;
                obj.solver = parameters.solver;
                obj.slope_ub = parameters.slope_ub;
                obj.t_max = parameters.t_max;
                obj.t_max_plot = parameters.t_max_plot;
                obj.linearX = parameters.linearX;
                obj.logX = parameters.logX;
                obj.linearY = parameters.linearY;
                obj.logY = parameters.logY;

                if parameters.cluster
                    obj.input_dir = '/bucket/DieckmannU/Clara/dir/';
                    obj.output_dir = '/bucket/DieckmannU/Clara/foodweb-dynamics-rotationproj-master/src/out/';
                else
                    obj.input_dir = '/Users/cseunitclara/Documents/OIST/dir/';
                    obj.output_dir = '/Users/cseunitclara/Documents/OIST/foodweb-dynamics-rotationproj-master/src/out';
                end
                if parameters.compute_locally_ko
                    obj.input_dir ='/Users/cseunitclara/Documents/OIST/dir/';
                    obj.output_dir = '/Users/cseunitclara/Documents/OIST/foodweb-dynamics-rotationproj-master/src/out';
                end
                
                if ~exist(obj.output_dir, 'dir')
                    disp(['Creating output directory: ', obj.output_dir]);
                    mkdir(obj.output_dir)
                    
                end
            end
        end
        
        function new_network_FileList = filter_networks(obj,old_network_FileList)
            include_networks = 1:obj.ALL_networks;
            include_networks(obj.exclude_networks_list) = [];

            filenames = {old_network_FileList(:).name};
            func = @(name) str2num(extractBefore(extractAfter(name,'timeInteg_ode23s_'),'_'));
            current_netIDs = cellfun(func, filenames, 'UniformOutput', true);

            for n=1:length(include_networks)
                ID = include_networks(n);
                if ~any(current_netIDs == ID)
                    fprintf('\nError! Missing network ID %i data', ID);
                else
                    pos = find(current_netIDs == ID);
                    new_network_FileList(n) = old_network_FileList(pos);
                end
            end
            new_network_FileList = new_network_FileList';
        end
      
        function b_init = get_initbiomass(obj, netInfo, options)
            BaseClass.set_random_seed;
            switch options.mode
                case 'empirical'
                    b_init = netInfo.b0;
                case 'normrand'
                    b_init = netInfo.b0 .* max(1e-8, normrnd(1,0.3,netInfo.n,1)); % The use of max function prevents negative biomass
                case 'percentrand'
                    rand_factor = [1-options.percent, 1+options.percent]; rand_int = randi([1, 2], netInfo.n, 1);
                    rand_int(rand_int == 1) = rand_factor(1); rand_int(rand_int == 2) = rand_factor(2);
                    b_init = netInfo.b0 .* rand_int;
                case 'lognormal'
                    mu = log(netInfo.b0/sqrt(1 + options.cv^2));
                    sigma = sqrt(log(1 + options.cv^2));
                    b_init = lognrnd(mu,repmat(sigma, netInfo.n, 1));
                case 'eigendirection'
                    % Evaluate Jacobian at the empirically-observed
                    % biomass
                    J = BaseClass.compute_jacobian_triple(netInfo, netInfo.b0, netInfo.b0);
                    [domeigval,domeigvec] = BaseClass.get_domeigcomp(J);

                    % Displace empirically-observed biomass until total
                    % biomass changes by 1%
                    percent_change = 0; b_init = netInfo.b0; step_size = 1e-5;
                    while percent_change < options.percent
                        % Displace compartments
                        new_b_init = b_init + step_size.*real(domeigvec);
                        if any(new_b_init < 0)
                            break;
                        else
                            b_init = new_b_init;
                            percent_change = abs(sum(b_init) - sum(netInfo.b0)) / sum(netInfo.b0);
                        end
                    end
            end

            if any(b_init < 0)
                fprintf('Initial biomass values: %s', mat2str(b_init'))
                error('Error! Found negative initial biomass values! Please check again...');
            end
        end
    end
end