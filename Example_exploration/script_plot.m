% A script that produces the plots for Example 2.
%
% For saving the figure, the toolbox export_fig
% (https://github.com/altmany/export_fig) is used.

save_flag = 0; % Flag for saving the results 

% Load results
filename = '';
if(~isempty(filename))
    disp('Loading file...')
    tic;
    load(filename)
    t = toc;
    disp(['    ... done! (',num2str(t),'s)'])
end

N_runs = numel(result_cell);
f_min = 0;

n = 10;
m = 6;

% Offset for arrays starting at 1 instead of 0. For easier comparison to
% paper notation.
nul = 1; 

%% (a)
figure;

lw = 2.2;
ms = 10;

for i = 1:N_runs
    eig_mat = result_cell{i}.eig_mat;

    eig_k_plus_1 = zeros(1,algo_options.N_iter);
    for k = 0:m-2
        eig_k_plus_1(nul+k) = eig_mat(k+1,nul+k);
    end

    rad = log10(norm(result_cell{i}.problem_data.x0));

    h1 = plot(rad,min(log10(eig_k_plus_1)),'k.','MarkerSize',ms);
    hold on
    h2 = plot(rad,max(log10(eig_mat(n,nul+(0:m-2)))),'ko','MarkerSize',5,'LineWidth',1.2);
end

legend([h1,h2],{'$\min_{k \in \{0,\dots,m-2\}} \lambda_{k+1}^k$','$\max_{k \in \{0,\dots,m-2\}} \lambda_n^k$'},...
    'Interpreter','latex','FontSize',15,'Location','southeast');

xlim([-30-0.5,2+0.5])
ylim([-0.6,1.6])
xlabel('$\| x^0 - x^* \|$','Interpreter','latex')
grid on
axis square

set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

% Log tick labeling (x)
old_ticks = xticks;
new_ticks_cell = cell(numel(old_ticks),1);
for i = 1:numel(old_ticks)
    new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
end
xticklabels(new_ticks_cell)

% Log tick labeling (y)
old_ticks = yticks;
new_ticks_cell = cell(numel(old_ticks),1);
for i = 1:numel(old_ticks)
    new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
end
yticklabels(new_ticks_cell)

if(save_flag)
    export_fig 'plot1' '-png' '-r500' '-transparent'
end

%% (b)
figure;

for i = 1:N_runs
    H_cell = result_cell{i}.H_cell;
    y_arr = result_cell{i}.y_arr;
    s_arr = result_cell{i}.s_arr;

    theta_mat = NaN(m-1,m-2);
    for k_ = 1:m-2
        for k = 0:k_-1
            theta_mat(nul+k,nul+k_) = norm(H_cell{nul+k_} * y_arr(:,nul+k),2);
        end
    end

    rad = log10(norm(result_cell{i}.problem_data.x0));
    h1 = plot(rad,log10(max(max(theta_mat))),'k.','MarkerSize',ms);
    hold on
end

legend(h1,{'Value (9)'},...
    'Interpreter','latex','FontSize',15,'Location','northwest');

xlim([-30-0.5,2+0.5])
ylim([-32,6])
xlabel('$\| x^0 - x^* \|$','Interpreter','latex')
grid on
axis square

set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

% Log tick labeling (x)
old_ticks = xticks;
new_ticks_cell = cell(numel(old_ticks),1);
for i = 1:numel(old_ticks)
    new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
end
xticklabels(new_ticks_cell)

% Log tick labeling (y)
old_ticks = yticks;
new_ticks_cell = cell(numel(old_ticks),1);
for i = 1:numel(old_ticks)
    new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
end
yticklabels(new_ticks_cell)

if(save_flag)
    export_fig 'plot2' '-png' '-r500' '-transparent'
end