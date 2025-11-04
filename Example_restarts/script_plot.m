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

addpath('../LW2019_84');

f_min = 0;

n = 100;
m = 80;

reset_period = m;
x_arr = result_m.x_arr;
f_arr = result_m.f_arr;
t_arr = result_m.t_arr;
f = problem_data.f;

lw = 1.5;
ms = 10;

reset_inds = 1:reset_period:size(x_arr,2);

%% (a)
figure;

plot(0:size(x_arr,2)-1,log10(result_m_minus_1.f_arr - f_min),'r.-','MarkerSize',ms,'LineWidth',lw);
hold on
plot(0:size(x_arr,2)-1,log10(f_arr - f_min),'k.-','MarkerSize',ms,'LineWidth',lw);
xline(reset_inds,'k--','LineWidth',1.25)

h1 = plot(-100,-100,'k.-','MarkerSize',13,'LineWidth',1.15); % Dummy plot for legend
h2 = plot(-100,-100,'r.-','MarkerSize',13,'LineWidth',1.15); % Dummy plot for legend
legend([h1,h2],{'$m$','$m-1$'},...
    'Interpreter','latex','FontSize',15,'Location','southwest');

xlim([0,1440])
ylim([-12.5,8])
xlabel('$k$','Interpreter','latex')
xticks(0:250:1500)
axis square

set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

% Log tick labeling (y)
old_ticks = yticks;
new_ticks_cell = cell(numel(old_ticks),1);
for i = 1:numel(old_ticks)
    new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
end
yticklabels(new_ticks_cell)

    % Subfigure -----------------------------------------------------------
    
    bk_ax = axes('position',[.45 .46 .35 .45]);
    set(gca,'XColor','none','YColor','none')
    
    ax = axes('position',[.45 .515 .4 .4]);
    
    min_ind = 801-20;
    max_ind = 961+20;
    zoom_inds = min_ind:max_ind;
    
    plot(zoom_inds-1,log10(f_arr(zoom_inds) - f_min),'k.-','MarkerSize',ms,'LineWidth',1);
    hold on
    xline(reset_inds(reset_inds >= min_ind & reset_inds <= max_ind)-1,'k--','LineWidth',1.25);
    
    xlim([min_ind,max_ind])
    ylim([-8.2,-6.3])
    axis square
    ax.ClippingStyle = "rectangle";
    
    xticks(800:100:1000)
    yticks([-9,-8,-7])
    
    % Log tick labeling (y)
    old_ticks = yticks;
    new_ticks_cell = cell(numel(old_ticks),1);
    for i = 1:numel(old_ticks)
        new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
    end
    yticklabels(new_ticks_cell)
    
    set(gca,'linewidth',1.1)
    set(gca,'fontsize',10)

if(save_flag)
    export_fig 'plot1' '-png' '-r500' '-transparent'
end

%% (b)
figure;

I_arr = zeros(1,size(x_arr,2));
for k = 1:size(x_arr,2)
    [~,tmp_I] = f(x_arr(:,k));
    I_arr(k) = tmp_I;
end

I_arr_cut = I_arr(1:floor(size(x_arr,2)/reset_period) * reset_period);
I_arr_mat = reshape(I_arr_cut,[reset_period,floor(size(x_arr,2)/reset_period)]);

I_unique = NaN(1,floor(size(x_arr,2)/reset_period));
for i = 1:floor(size(x_arr,2)/reset_period)
    I_unique(i) = numel(unique(I_arr_mat(:,i)));
end

plot(1:floor(size(x_arr,2)/reset_period),I_unique,'k.-','MarkerSize',12,'LineWidth',lw);
hold on
yline(80,'k--','LineWidth',1.25)

grid on

xlim([1,18])
ylim([29,83])
xlabel('Number of resets','Interpreter','latex');

axis square

h1 = plot(-100,-100,'k.-','MarkerSize',13,'LineWidth',1.15); % dummy
legend(h1,{'Unique sel. funs.'},...
    'Interpreter','latex','FontSize',15,'Location','southeast');

set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

if(save_flag)
    export_fig 'plot2' '-png' '-r500' '-transparent'
end

%% (c)
figure;

h1 = plot(0:size(x_arr,2)-2,log10(t_arr),'k.','MarkerSize',11,'LineWidth',lw);
hold on
xline(reset_inds,'k--','LineWidth',1.25)

legend(h1,{'$t_k$'},...
    'Interpreter','latex','FontSize',15,'Location','southwest');

set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

xlim([0,1441])
ylim([-18,0])
xlabel('$k$','Interpreter','latex');
xticks(0:250:1500)

axis square

set(gca,'linewidth',1.1)
set(gca,'fontsize',15)

% Log tick labeling (y)
old_ticks = yticks;
new_ticks_cell = cell(numel(old_ticks),1);
for i = 1:numel(old_ticks)
    new_ticks_cell{i} = ['10^{',num2str(old_ticks(i)),'}'];
end
yticklabels(new_ticks_cell)

if(save_flag)
    export_fig 'plot3' '-png' '-r500' '-transparent'
end