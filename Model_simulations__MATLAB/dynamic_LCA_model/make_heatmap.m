function make_heatmap(output)
% Plot simulated within- and across-dimension adaptation heatmaps.
% Author: Davide Gheza

% Accuracy effects.
figure(101)

d = output.prevcon_task_int_acc-output.previnc_task_int_acc;
d = round(d,3);
heatmap({'C','S','E','P'},{'C','S','E','P'},mean(d,3))
set(gca,'fontsize',14)
colormap turbo
xlabel('Previous non-cued dimension')
ylabel('Current non-cued dimension')
colorbar('off');

% Custom colormap: red to white
nColors = 256;
customMap = [linspace(1,1,nColors)', linspace(0,1,nColors)', linspace(0,1,nColors)']; % R from 1->1, G&B from 1->0
colormap(customMap)

% % Save out accuracy panel
% aspectRatio = 1;
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% newWidth = screenposition(4)/1.5 * aspectRatio;  % width and height based on aspect ratio
% newHeight = screenposition(4)/1.5;  % Keep the height the same
% set(gcf, ...
%     'PaperPosition', [0 0 newWidth newHeight], ...
%     'PaperSize', [newWidth newHeight]);
% 
% filename = 'acc_sims_heatmap_accuracy';
% print('-dpdf', '-vector', filename);




% RTs
figure(100)

d = output.prevcon_task_int_rt-output.previnc_task_int_rt;
heatmap({'C','S','B','P'},{'C','S','B','P'},mean(d,3))
set(gca,'fontsize',14)
colormap turbo
xlabel('Previous non-cued dimension')
ylabel('Current non-cued dimension')

% Custom colormap: red to white
nColors = 256;
customMap = [linspace(1,1,nColors)', linspace(0,1,nColors)', linspace(0,1,nColors)']; % R from 1->1, G&B from 1->0
colormap(customMap)

end

