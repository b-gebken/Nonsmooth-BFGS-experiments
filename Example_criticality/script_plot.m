% A script that produces the plots for Example 1.
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

%% (a)
figure;

lw = 2.2;
ms = 10;

colororder("gem12");

for i = 1:N_runs
    plot(0:size(result_cell{i}.x_arr,2)-1,log10(result_cell{i}.f_arr - f_min),'-','LineWidth',lw);
    hold on;
end

h1 = plot(-100,100,'k-','LineWidth',lw); % Dummy plot for legend

legend(h1,{'$f(x^k) - f(x^*)$'},...
    'Interpreter','latex','FontSize',15,'Location','southwest');

xlim([0,algo_options.N_iter])
ylim([-90,10])
xlabel('$k$','Interpreter','latex')
xticks(0:200:1000)
grid on
axis square

set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

% Log tick labeling
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

ind = min(m-1,n);
for i = 1:N_runs
    eig_mat = result_cell{i}.eig_mat;
    h1 = plot(0:size(result_cell{i}.x_arr,2)-1,log10(eig_mat(ind,:)),'r-','LineWidth',lw);
    hold on

    if(ind+1 <= n)
        h2 = plot(0:size(result_cell{i}.x_arr,2)-1,log10(eig_mat(ind+1,:)),'b-','LineWidth',lw);
    end

    h3 = plot(0:size(result_cell{i}.x_arr,2)-1,log10(eig_mat(n,:)),'g-','LineWidth',lw);
end

legend([h1,h2,h3],{['$\lambda^k_',num2str(ind),'$'],['$\lambda^k_',num2str(ind+1),'$'],['$\lambda^k_{',num2str(n),'}$']},...
    'Interpreter','latex','FontSize',15,'Location','southwest');

xlim([0,algo_options.N_iter])
ylim([-85,10])
xlabel('$k$','Interpreter','latex')
xticks(0:200:1000)
grid on
axis square

set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

% Log tick labeling
old_ticks = yticks;
new_ticks_cell = cell(numel(old_ticks),1);
for i = 1:numel(old_ticks)
    new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
end
yticklabels(new_ticks_cell)

if(save_flag)
    export_fig 'plot2' '-png' '-r500' '-transparent'
end

%% (c)
figure;

colororder("gem12");

for i = 1:N_runs
    plot(0:size(result_cell{i}.x_arr,2)-1,log10(result_cell{i}.secant_mem_arr),'-','LineWidth',lw);
    hold on
end

% dummy plot
h1 = plot(-100,100,'k-','LineWidth',lw);

legend(h1,{'Value (5)'},...
    'Interpreter','latex','FontSize',15,'Location','southwest');

xlim([0,algo_options.N_iter])
ylim([-45,5])
xlabel('$k$','Interpreter','latex')
xticks(0:200:1000)
grid on
axis square

set(gca,'linewidth',1.1);
set(gca,'fontsize',15);

% Log tick labeling
old_ticks = yticks;
new_ticks_cell = cell(numel(old_ticks),1);
for i = 1:numel(old_ticks)
    new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
end
yticklabels(new_ticks_cell)

if(save_flag)
    export_fig 'plot3' '-png' '-r500' '-transparent'
end