% A simple example with a nonsmooth function from [Guo, Lewis (2018)].

%% Define objective function and initial point
n = 2;
f = @(x) x(1,:).^2 + abs(x(2,:));
grad_f = @(x) [2*x(1,:); sign(x(2,:))];

problem_data.f = f;
problem_data.grad_f = grad_f;
problem_data.n = n;
problem_data.x0 = [1;2];

%% Set options for BFGS method

algo_options.N_iter = 50;
algo_options.H0 = eye(n);
algo_options.c1 = 0.0001;
algo_options.c2 = 0.5;
algo_options.reset_period = Inf;
algo_options.step_threshold = 10^-15;
algo_options.descent_threshold = 10^-15;
algo_options.disp_flag = 1;

%% Run BFGS method

% To compute with variable-precision arithmetic (vpa), the input structs
% have to be converted.
% [problem_data,algo_options] = vpa_wrapper(problem_data,algo_options,500);

result = BFGS_method(problem_data,algo_options);

x_arr = result.x_arr;
f_arr = result.f_arr;
grad_arr = result.grad_arr;
H_cell = result.H_cell;
t_arr = result.t_arr;
y_arr = result.y_arr;
s_arr = result.s_arr;
p_arr = result.p_arr;

%% Plots

x_min = zeros(n,1);
f_min = f(x_min);

figure;
lw = 1.5;
ms = 12;

subplot(1,2,1) % ----------------------------------------------------------
h1 = plot(0:size(x_arr,2)-1,log10(vecnorm(double(x_arr - x_min),2,1)),'.-','LineWidth',lw,'MarkerSize',ms);
hold on
h2 = plot(0:size(x_arr,2)-1,log10(f_arr - f_min),'.-','LineWidth',lw,'MarkerSize',ms);

grid on
title('Optimality')
legend([h1,h2],{'$|| x^k - x^* ||$','$f(x^k) - f(x^*)$'},'Location','northeast','Interpreter','latex','FontSize',12);
xlabel('k')

subplot(1,2,2) % ----------------------------------------------------------
hold on;

eig_mat = real(cell2mat(cellfun(@(in) sort(double(eig(in))),H_cell,'UniformOutput',false)));
for i = 1:n
    plot(0:size(x_arr,2)-1,log10(abs(eig_mat(i,:))),'.-','LineWidth',lw,'MarkerSize',ms);
end

I = min(eig_mat,[],1) < 0;
if(any(I))
    xline(find(I),'r--');
end

I = min(eig_mat,[],1) == 0;
if(any(I))
    xline(find(I),'k--');
end

grid on
title('Eigenvalues of H_k')
xlabel('k')
