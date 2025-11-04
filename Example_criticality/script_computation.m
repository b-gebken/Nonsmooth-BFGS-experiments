% A script that performs the numerical experiment shown in Example 1 and
% saves the result.
%
% The runtime is relatively long due to the use of vpa. For much faster but
% less accurate computation, change symbolic_flag from 1 to 0.

symbolic_flag = 1; % Flag for symbolic or non-symbolic computation
save_flag = 1; % Flag for saving the results 

%% Example parameters
rng("default")

N_runs = 10;
result_cell = cell(1,N_runs);

% Add path to functions for Problem 8.4 in [Lewis, Wylie (2019)]
addpath('../LW2019_84');

n = 10;
m = 6;

%% Set algorithm options

algo_options.N_iter = 1000;
algo_options.c1 = 0.0001;
algo_options.c2 = 0.5;
algo_options.reset_period = Inf;
algo_options.step_threshold = 0;
algo_options.descent_threshold = 0;
algo_options.disp_flag = 1;

H0 = 2*rand(n) - 1; H0 = H0'*H0;
algo_options.H0 = H0;


%% Apply BFGS method
addpath('../BFGS_method');

problem_data.n = n;

for i = 1:N_runs

    % Generate random problem data
    problem_data.x0 = 2*rand(n,1) - 1;

    M_cell = cell(1,m);
    for j = 1:m
        tmp = 2*rand(n)-1;
        M_cell{j} = tmp * tmp';
    end
    d_arr = rand(1,m);

    if(symbolic_flag)
        % Generate g symbolically and convert the structs
        digits(500);
        lambda = vpa(rand(m,1)); lambda = lambda/sum(lambda);
        tmp = 2*vpa(rand(n,m))-1; g_arr = tmp - tmp*lambda;
        for j = 1:m
            M_cell{j} = vpa(M_cell{j});
        end
        [problem_data,algo_options] = vpa_wrapper(problem_data,algo_options,500);
    else
        % Generate g non-symbolically (sufficient for default accuracy, but
        % yields negative f values in vpa) 
        lambda = rand(m,1); lambda = lambda/sum(lambda);
        tmp = 2*rand(n,m)-1; g_arr = tmp - tmp*lambda;
    end

    problem_data.f = @(x) lw2019_84_f(x,g_arr,M_cell,d_arr);
    problem_data.grad_f = @(x) lw2019_84_grad(x,g_arr,M_cell,d_arr);
    problem_data.gradselfuns = @(x) lw2019_84_gradselfuns(x,g_arr,M_cell,d_arr);

    % Run BFGS method
    result_cell{i} = BFGS_method(problem_data,algo_options);

    result_cell{i}.problem_data = problem_data;
    result_cell{i}.algo_options = algo_options;

end

%% Additional computations for plots

x_min = zeros(n,1);

disp('Computing eigenvalues of H...')
ind = min(m-1,n);
for i = 1:N_runs
    fprintf([repmat('\b',1,6*(i > 1)),'%.4f',repmat('\n',1,i == N_runs)], i/N_runs);

    eig_mat = real(cell2mat(cellfun(@(in) double(sort(eig(in))),result_cell{i}.H_cell,'UniformOutput',false)));
    result_cell{i}.eig_mat = eig_mat;
end

disp('Computing gradients of all selection functions in all x...')
for i = 1:N_runs
    fprintf([repmat('\b',1,6*(i > 1)),'%.4f',repmat('\n',1,i == N_runs)], i/N_runs);
    N_iter_i = size(result_cell{i}.x_arr,2);
    secant_mem_arr = zeros(1,N_iter_i);
    grad_normal_arr = zeros(1,N_iter_i);

    grads_opt = result_cell{i}.problem_data.gradselfuns(x_min);

    fprintf(' ');
    for j = 1:N_iter_i
        fprintf([repmat('\b',1,6*(j > 1)),'%.4f',repmat('\n',1,j == N_iter_i)], j/N_iter_i);
        secant_mem_arr(j) = max(vecnorm(double( result_cell{i}.H_cell{j} * (grads_opt(:,2:end) - grads_opt(:,1)) ),inf,1));
        grad_normal_arr(j) = max(vecnorm(double( result_cell{i}.H_cell{j} * grads_opt ),inf,1));
    end
    fprintf(repmat('\b',1,6+2));

    result_cell{i}.secant_mem_arr = secant_mem_arr;
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