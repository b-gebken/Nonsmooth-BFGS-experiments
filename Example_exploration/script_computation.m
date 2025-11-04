% A script that performs the numerical experiment shown in Example 2 and
% saves the result.
%
% The runtime is relatively long due to the use of vpa. For much faster but
% less accurate computation, change symbolic_flag from 1 to 0.

symbolic_flag = 1; % Flag for symbolic or non-symbolic computation
save_flag = 1; % Flag for saving the results 

%% Example parameters
rng("default")

N_runs = 100;
result_cell = cell(1,N_runs);

% Add path to functions for Problem 8.4 in [Lewis, Wylie (2019)]
addpath('../LW2019_84');

n = 10;
m = 6;

%% Set algorithm options

algo_options.N_iter = m-1;
algo_options.c1 = 0.0001;
algo_options.c2 = 0.5;
algo_options.reset_period = Inf;
algo_options.step_threshold = 0;
algo_options.descent_threshold = 0;
algo_options.disp_flag = 0;
algo_options.H0 = eye(n);

%% Generate random problem data

problem_data.n = n;

% Array of radii for the initial points
radii_arr = linspace(2,-30,N_runs);

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
else
    % Generate g non-symbolically (sufficient for default accuracy, but
    % yields negative f values in vpa)
    lambda = rand(m,1); lambda = lambda/sum(lambda);
    tmp = 2*rand(n,m)-1; g_arr = tmp - tmp*lambda;
end

problem_data.f = @(x) lw2019_84_f(x,g_arr,M_cell,d_arr);
problem_data.grad_f = @(x) lw2019_84_grad(x,g_arr,M_cell,d_arr);

%% Apply BFGS method
addpath('../BFGS_method');

for i = 1:N_runs
    fprintf([repmat('    ',1,i == 1),repmat('\b',1,6*(i > 1)),'%.4f',repmat('\n',1,i == N_runs)], i/N_runs);
    
    % Generate random initial point with given distance to minimum
    x0 = randn(n,1); x0 = x0/norm(x0,2); % https://en.wikipedia.org/wiki/N-sphere#Uniformly_at_random_on_the_(n_%E2%88%92_1)-sphere
    problem_data.x0 = 10^radii_arr(i) * x0;

    % vpa wrapper
    if(symbolic_flag)
        [problem_data,algo_options] = vpa_wrapper(problem_data,algo_options,500);
    end

    % Run BFGS method
    result_cell{i} = BFGS_method(problem_data,algo_options);

    result_cell{i}.problem_data = problem_data;
    result_cell{i}.algo_options = algo_options;
end

%% Additional computations for plots

disp('Computing eigenvalues of H...')
for i = 1:N_runs
    fprintf([repmat('    ',1,i == 1),repmat('\b',1,6*(i > 1)),'%.4f',repmat('\n',1,i == N_runs)], i/N_runs);

    eig_mat = real(cell2mat(cellfun(@(in) double(sort(eig(in))),result_cell{i}.H_cell,'UniformOutput',false)));
    result_cell{i}.eig_mat = eig_mat;
end

%% Save the results
if(save_flag)
    cDT = datetime('now');
    cDT.Format = 'dd-MMM-yyyy_HH_mm_ss';
    filename = string(cDT);

    save(filename,'problem_data','algo_options','result_cell');
end

%% Paper plots

script_plot