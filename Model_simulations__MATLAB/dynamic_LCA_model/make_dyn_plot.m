function make_dyn_plot(output)
% Plot within-trial energy trajectories by dimension and congruency.
% Author: Davide Gheza

c = spring;

figure

c4 = [1,86,171,256]; % default colors if 4 dimensions
tasks = length(unique(output.task));
if tasks>4
    c4 = sort(round(rand(tasks,1)*256));
end

max_rt = max(output.rt);

for i = 1:tasks
    plot(1:max_rt,squeeze(output.path_en_con(1:max_rt,i,1)),'LineWidth',2.5,'MarkerSize',15,'Color',c(c4(i),:),'LineStyle',":");
    hold on
    plot(1:max_rt,squeeze(output.path_en_con(1:max_rt,i,2)),'LineWidth',2.5,'MarkerSize',15,'Color',c(c4(i),:),'LineStyle',"-");
    % plot(1:max_rt,squeeze(output.path_en_con(1:max_rt,i,1)-output.path_en_con(1:max_rt,i,2)),'LineWidth',2.5,'MarkerSize',15,'Color',c(c4(i),:),'LineStyle','--');
end
hold off

xlabel('Time (cycles)')
ylabel('Energy')

cut_cycles = 0.3; % percentage of cycles of max_rt to keep (long RTs are rare and average is noisy)
set(gca,'FontName','Arial','fontsize',14)
xlim([0 max_rt*cut_cycles])

taskLegends = arrayfun(@(x) sprintf('Non-cued dim. %d', x), 1:tasks, 'UniformOutput', false);
for i = 1:length(taskLegends)
    Legends{(i*(i==1))+(i*2-1)*(i>1)} = [taskLegends{i} ' - Incongruent'];
    Legends{(i*(i==1))+(i*2-1)*(i>1)+1} = [taskLegends{i} ' - Congruent'];
end
legend(Legends, 'Location', 'best')


end
