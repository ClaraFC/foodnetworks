function dynamic_model()
close all
plotting = true;
printing = true;
input_dir =  '../dat';
output_dir = '../out';

if ~exist(output_dir, 'dir')
    mkdir(output_dir)
    disp(['Created output directory: ', output_dir]);
end

filepaths = get_input_filepaths(input_dir);

res = {};
Delta = 0; %0.01; % if the stock falls below a product of initial stock and Delta we treat it as extinct

t_max = 1e4;
if (t_max < 1e4)
    warning('t_max < 1e4')
end
tc = 30; % Plots are generated from t=0 to t=tc.

file_number = length(filepaths);

for i = 118 %1:file_number
    filepath = filepaths{i};
    
    [name, crc, n, l, b0, p, q, r, F] = read_input_file(filepath);
    
    disp([num2str(i), '/', num2str(length(filepaths)), ' (', ...
        num2str(round(100*i/file_number)), '%) ', filepath ', ', ...
        name, ', n = ', num2str(n)])
    
    % Here we envisage two for loops over a grid of (s_d, s_r)
    s_d = 1;
    s_r = 0;
    db_dt = check_fixedpoint(s_d, b0, b0, F, q, r, p, n);
    disp(db_dt)
    disp('Empirical biomass')
    disp(b0)
    J = compute_jacobian(s_d, s_r, b0, F, q, r, n);
    
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8,'InitialStep',0.001, 'Refine', 200);
    
    b_init = b0 .* normrnd(1,0.3,n,1); % pertubed initial biomass
    
    [t_out, b_out, b_steady] = simulate_dynamics(b0, b_init, J, p, q, r, F, t_max, n, l, options, Delta);
    analyse_ODE_solution(t_out, b_out, tc, b_steady, 'Stable state simulation', ...
        plotting, printing);
    
    res{i} = struct('crc', crc,...
        'name', name, ...
        'min_b_steady', min(b_steady) ... % or whatever we want to report
        );
    
    % For plotting, ternary diagrams can be useful:
    % https://www.mathworks.com/matlabcentral/fileexchange/2299-alchemyst-ternplot
    
    %     if (i >= 3)
    %         break
    %     end
    
end

save([output_dir, '/res' '_' datestr(now,'yyyymmdd') '.mat']);
save_as_csv(res, output_dir, t_max);

end

function J = compute_jacobian(s_d, s_r, b0, F, q, r, n)
if ~(s_d>=0 && s_d<= 1 && s_r >=0 && s_r <= 1)
    error('Incorrect value of parameter s_d or s_r')
end

if (s_d == 1) % donor-control model
    Cd = F ./ repmat(b0.', n, 1);
    J = Cd - diag(sum(Cd).' + q./b0 + r./b0);
end
% Eventually, we would like to have the general case implemented for any
% feasible values of s_d and s_r
end

function nres = analyse_ODE_solution(t_out, b_out, tc, b_steady, text, plotting, printing)
global extinction_times;
if plotting
    plot_solution(t_out, b_out, b_steady, tc);
end
if printing
    disp(text);
    disp('Extinction times:');
    disp(extinction_times);
    disp('b_steady:')
    disp(b_steady.');
    disp('b_out final:')
end
b_out_final = b_out(end, :).';

res = b_out_final - max(b_steady, 0);
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

function [t_out, b_out, b_steady, ext_times] = simulate_dynamics(b0, b_init, J, p, q, r, F, t_max, n, l, options, Delta)
global extinction_times;
extinction_times = Inf(1, n);
set_extinction_flags(b_init == 0, 0);

[t_out, b_out] = solve_ODE(b0, b_init, J, p, q, r, F, n, t_max, l, options, Delta);

b_steady = b0; %asymptotic_steady_state(J, p, t_max, false, false);

if any(b_steady < 0)
    warning('Negative steady-state solution');
    b_steady
    b_init
    J
    p
end

ext_times = extinction_times;

ext_indicator = t_out >= repmat(ext_times, size(t_out, 1), 1);
b_out(ext_indicator) = 0;

end

function [t_out, b_out] = solve_ODE(b0,b_init, J, p, q, r, F, n, t_max, l, options, Delta)
t_span = [0 t_max];
[t_out, b_out] = ode23s(@(t,b) differential(t,b,J,p,q,r,l,b0,b_init,F,n,Delta) , t_span, b_init, options);
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

function [times, stocks, steady_state] = read_ecomodels_trajectories(filepath)
run(filepath);
end

function plot_solution(t_out, b_out, b_steady, tc)
% i = 1;
% b_out = b_out(:,i); b_steady = b_steady(:,i);
ind = t_out <= tc;
ind(1) = true;
h1 = plot(t_out(ind), b_out(ind, :), '.-');
xlabel('t');
ylabel('b(t)');
y_lim = get(gca, 'ylim');
set(gca, 'YScale', 'log');
%y_lim(1) = 0;
set(gca, 'ylim', y_lim);
for (i = 1:size(b_steady, 1))
    hl(i) = line([0; max([t_out(ind).'])], [b_steady(i), b_steady(i)]);
    set(hl(i), 'LineStyle', ':');
end
legend([h1(1), hl(1)], 'Integration', 'Steady state');
drawnow;
end

function set_extinction_flags(extinction, time)
global extinction_flags;
global extinction_times;
extinction_flags = extinction;
ext_time = min(time, extinction_times);
extinction_times(extinction) = ext_time(extinction);
end

function extinction = get_extinction_flags(time)
global extinction_flags;
global extinction_times;
extinction = extinction_flags;
extinction(extinction_times > time) = false;
end

function b_steady = asymptotic_steady_state(J, p, t, dont_remove_extinct_nodes, fix_negative_b0)
if nargin < 4
    fix_negative_b0 = false;
    extinction = get_extinction_flags(t);
    [J, p, ~] = remove_extinct_compartments(J, p, NaN('like', p), extinction);
end
b_steady = -pinv(J) * p;
if(fix_negative_b0)
    b_steady(b_steady <= 0) = 0;
end
end


function db = differential(t, b, J, p, q, r, l, b0, b_init, F, n, Delta)

update_extinction_flags(b, l, t, b_init, Delta);
extinction = get_extinction_flags(t);
[J, p, b] = remove_extinct_compartments(J, p, b, extinction);

s_d = 1;
db = check_fixedpoint(s_d, b0, b, F, q, r, p, n); 
%db = J * b + p;
depleted = b <= 0 & (1:length(b) > l).';
db(depleted) = - b(depleted); % max(0, dx(depleted));
end

function [J, p, b] = remove_extinct_compartments(J, p, b, extinction)
J_non_diag =  J .* (1 - eye(size(J)));
if (sum(extinction) > 0)
    
    if (sum(extinction) == 1)
        J_diag = diag(J) + J_non_diag(extinction,:).';
    elseif (sum(extinction) > 1)
        J_diag = diag(J) + sum(J_non_diag(extinction,:)).';
    end
    
    J = J_non_diag + diag(J_diag);
    
    J(extinction, :) = 0;
    J(:, extinction) = 0;
    p(extinction) = 0;
    b(extinction) = 0;
end
end

function update_extinction_flags(b, l, t, b_init, Delta)
extinction = get_extinction_flags(t);

extinction = extinction | (b < b_init * Delta & ((1:length(b)).' <= l));
set_extinction_flags(extinction, t);

end



function filepaths = get_input_filepaths(input_dir)
files = dir(input_dir);
filepaths = {};
for i = 1:length(files)
    if endsWith(files(i).name, '.m')
        filepaths{length(filepaths)+1} = fullfile(input_dir, files(i).name);
    end
end
end

function [str] = save_as_csv(res, output_dir, t_max)
fns = fieldnames(res{length(res)});
str = strjoin(fns,';');
for i = 1:length(res)
    str = [str '\n'];
    fns = fieldnames(res{i});
    for j = 1:length(fns)
        fn = fns{j};
        % disp(fn);
        str = [str num2str(getfield(res{i},fn))];
        if j < length(fns)
            str = [str ';'];
        end
    end
end
formatOut = 'yyyymmdd';
filepath = [output_dir, '/res_tmax_' num2str(t_max) '_' datestr(now,formatOut) '.csv'];
fileID = fopen(filepath,'w');
fprintf(fileID, str);
fclose(fileID);
disp(['Results saved in file ' filepath]);
end

function db_dt = check_fixedpoint(s_d, b0, b, F, q, r, p, n)
    Cd = F ./ repmat(b0.', n, 1);
    F_th_ij = s_d * (Cd .* repmat(b', n, 1));
    F_th_ji = s_d * (Cd .* repmat(b', n, 1));
    db_dt = sum(F_th_ij, 2) - sum(F_th_ji, 1)' - (q./b0 + r./b0).*b + p;
%     db_dt = sum(F_th_ij, 2) - sum(F_th_ji, 1)' - (q + r) + p;
end