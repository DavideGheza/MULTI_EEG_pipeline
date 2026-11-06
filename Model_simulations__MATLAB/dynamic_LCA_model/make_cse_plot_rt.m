function make_cse_plot_rt(output,lag)
% Plot simulated RT by current congruency and lagged congruency.
% Author: Davide Gheza and Wouter Kool
figure
c = cool;

if nargin<2 || lag == 1
    means = (output.prev_sum_rt(:,:));
elseif lag == 2
    means = (output.prev_2_sum_rt(:,:));
elseif lag == 3
    means = (output.prev_3_sum_rt(:,:));
end

c4 = [1,86,171,256]; % default colors if 4 dimensions
tasks = length(unique(output.task));
if tasks>4
    c4 = sort(round(rand(tasks,1)*256));
end

for i = 1:tasks
    plot(0:tasks-1,means(i,:),'LineWidth',1.5,'MarkerSize',15,'Color',c(c4(i),:));
    hold on
end
hold off

xlabel('# congruent dimensions')
ylabel('Reaction time')

set(gca,'FontName','Arial','fontsize',14,'Xtick',0:1:tasks-1)
xlim([-0.5 tasks-0.5])

legend(arrayfun(@num2str, 0:tasks, 'UniformOutput', false))



end
