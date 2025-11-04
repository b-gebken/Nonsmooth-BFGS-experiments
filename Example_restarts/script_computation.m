% A script that performs the numerical experiment shown in Example 3 and
% saves the result.
%
% Computation with vpa was not feasible due to runtime, so the default
% accuracy is used for this experiment.

symbolic_flag = 0; % Flag for symbolic or non-symbolic computation
save_flag = 1; % Flag for saving the results 

%% Example parameters
rng("default")

N_runs = 10;
result_cell = cell(1,N_runs);

% Load functions for Problem 8.4 in [Lewis, Wylie (2019)]
addpath('../LW2019_84');

n = 100;
m = 80;

%% Set algorithm options

algo_options.N_iter = 1440;
algo_options.c1 = 0.5;
algo_options.c2 = 0.75;

algo_options.step_threshold = 0;
algo_options.descent_threshold = 0;
algo_options.disp_flag = 1;
algo_options.H0 = eye(n);

%% Generate random problem data

problem_data.n = n;
problem_data.x0 = 5*(2*rand(n,1) - 1);

M_cell = cell(1,m);
for j = 1:m
    tmp = 2*rand(n)-1;
    M_cell{j} = tmp * tmp';
end
d_arr = rand(1,m);

if(symbolic_flag)
    % Generate g symbolically
    digits(500);
    lambda = vpa(rand(m,1)); lambda = lambda/sum(lambda);
    tmp = 2*vpa(rand(n,m))-1; g_arr = tmp - tmp*lambda;
    for j = 1:m
        M_cell{j} = vpa(M_cell{j});
    end
    d_arr = vpa(d_arr);
    [problem_data,algo_options] = vpa_wrapper(problem_data,algo_options,500);
else
    % Generate g non-symbolically (sufficient for default accuracy, but
    % yields negative f values in vpa)
    lambda = rand(m,1); lambda = lambda/sum(lambda);
    tmp = 2*rand(n,m)-1; g_arr = tmp - tmp*lambda;
end

problem_data.f = @(x) lw2019_84_f(x,g_arr,M_cell,d_arr);
problem_data.grad_f = @(x) lw2019_84_grad(x,g_arr,M_cell,d_arr);

%% Run BFGS for different reset patterns
addpath('../BFGS_method');

algo_options.reset_period = m;
result_m = BFGS_method(problem_data,algo_options);

algo_options.reset_period = m-1;
result_m_minus_1 = BFGS_method(problem_data,algo_options);

%% Save in a file
if(save_flag)
    cDT = datetime('now');
    cDT.Format = 'dd-MMM-yyyy_HH_mm_ss';
    filename = string(cDT);

    save(filename,'problem_data','algo_options','result_m','result_m_minus_1');
end

%% Paper plots

script_plot