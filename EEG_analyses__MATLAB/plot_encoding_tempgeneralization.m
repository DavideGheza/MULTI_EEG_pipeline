%% Plot temporal generalization of full-model encoding axes
%% Used for quick inspection only. No data or plots are saved out.
% Loads repeated-cross-validation encoding outputs and plots train-time by
% test-time projections for a selected axis and proportion-congruency condition. Available
% axes (`model`): SR, SI, dimension-general or dimension-specific conflict,
% previous conflict, and cue.
% Author: Davide Gheza

%% 1. Load settings

clear

% choose PC
PC_idx = [1,4];

model='DSC'; % 'SR' 'SI' 'DGC' 'DSC' 'PrevC' 'Cue'
if strcmp('SR',model)
    ax_name = 'encoded_dim_spec_SR';
elseif strcmp('SI',model)
    ax_name = 'encoded_dim_gen_SI';
elseif strcmp('DGC',model)
    ax_name = 'encoded_dim_gen_conf';
elseif strcmp('DSC',model)
    ax_name = 'encoded_dim_spec_conf';
elseif strcmp('PrevC',model)
    ax_name = 'encoded_prev_conf';
elseif strcmp('Cue',model)
    ax_name = 'encoded_cue';
end

noncued_contrast_all = nan(4,230,230,4,26);
cued_proj_all = nan(4,230,230,26); % store trial average per sbj to inpsect axes
noncued_proj_all = nan(4,230,230,4,26);



%Select onset or response events
onset = true;
if onset
    suffix = '';
else
    suffix = 'response_';   % response-locked
end

%Struct 'sets' is loaded from the 'scripts_and_settings' folder:
if ispc
    load('C:\Users\gdavide\OneDrive\Documenti\MATLAB\dg_workspace22\dimensions\analyses\scripts_and_settings\preprocessing_settings.mat');
else
    load('/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/scripts_and_settings/preprocessing_settings_mac.mat');
end


%Add path to sub-folder of 'main_dir' to store datasets:
sets.rsa_data_dir=[sets.main_saving_dir 'rsa_data/'];

% sbjs
sbjIDX = [22,24:48];
% tf times
times = -800:10:1498;


%% loop through subjects
% parsing cued/noncued and trial averaging step was done within encoding script, to reduce dimension of file

for sn = 1:length(sbjIDX)
    
% load single trial projections (~ RSA fits)
load(fullfile(sets.rsa_data_dir, ['RSA_pc/multiclass/singletrial/encoding_fullmodel_proj_wtempgen/encoding_fullmodel_' suffix 'sbj_' num2str(sbjIDX(sn)) '.mat']));
    cued_proj_all(:,:,:,sn) = cued.(ax_name);          
    noncued_proj_all(:,:,:,:,sn) = noncued.(ax_name);
    noncued_contrast_all(:,:,:,:,sn) = prevc_by_axis.(ax_name).contrast_prevI_minus_prevC;% note: included only in latest version of fullmodel encoding 





    fprintf('\tLoaded sbj %g\n', sn) 
end


   


%% plot simple projections timeseries, i.e., on diagonal of t train and t test 
% (across PC, for cued / each non cued dimension)
% note: for conflict models, cued is expected to be zero (no congruency for cued dim)


fontsize = 18;
colors = parula(4);

% extract simple time-course of train and test on same time: reshape dims 2&3, take diagonal
diag_lin_idx = sub2ind([230 230], 1:230, 1:230);
resh = reshape(cued_proj_all, 4, 230*230, 26); % collapse dims 2&3 - Matlab is column-major
cued_diag_vals = resh(:, diag_lin_idx, :);
% same for non-cued
resh = reshape(noncued_proj_all, 4, 230*230, 4, 26); % collapse dims 2&3 - Matlab is column-major
non_cued_diag_vals = resh(:, diag_lin_idx, :, :);
% average over PCs
cued_diag_vals = squeeze(mean(cued_diag_vals(PC_idx,:,:),1, 'omitnan'));      
non_cued_diag_vals = squeeze(mean(non_cued_diag_vals(PC_idx,:,:,:),1, 'omitnan')); % average over PCs
mean_vals = squeeze(mean(cued_diag_vals, 2)); 
mean_vals_2 = squeeze(mean(non_cued_diag_vals, 3)); 
sem_vals = squeeze(std(cued_diag_vals, 0, 2)) / sqrt(size(cued_diag_vals, 2));
sem_vals_2 = squeeze(std(non_cued_diag_vals, 0, 3)) / sqrt(size(non_cued_diag_vals, 3));

% Create the figure and plot
figure

hold on
% Plot mean and SEM as shaded areas
shaded_area_x = [times, fliplr(times)];
shaded_area_y = [mean_vals' - sem_vals', fliplr(mean_vals' + sem_vals')];
plot(times, mean_vals, 'LineWidth', 2, 'LineStyle', '-', 'Color', [0,0,0]);
fill(shaded_area_x, shaded_area_y, [0.5,0.5,0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

for i=1:4
        shaded_area_x = [times, fliplr(times)];
        shaded_area_y = [squeeze(mean_vals_2(:,i))' - squeeze(sem_vals_2(:,i))', ...
                        fliplr(squeeze(mean_vals_2(:,i))' + squeeze(sem_vals_2(:,i))')];
        fill(shaded_area_x, shaded_area_y, ...
                 colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(times, mean_vals_2(:,i), 'LineWidth', 2, 'LineStyle', '--', 'Color', colors(i,:));
end

        
        plot(times, zeros(length(times), 1), 'Color', 'black', 'LineStyle', ':', 'LineWidth', 2)        % Plot the zero line
        xlim([times(1) times(end)])
xlabel('Time [ms]', 'FontSize', fontsize)
ylabel('Projection', 'FontSize', fontsize)
set(gca, 'FontSize', fontsize)
legend('cued', '', 'non cued (c)','', 'non cued (s)','', 'non cued (e)','', 'non cued (p)', 'location', 'southeast')        % Add legend
hold off
ylim([-0.03 0.04])









 %% temporal generalization 
        

colors = parula(4);
newtitles = {'C', 'S', 'E', 'P'};  

fontsize = 18;

% extract slice
cued_slice_vals = squeeze(mean(cued_proj_all(PC_idx,:,:,:),[1,4]));      % mean over PC and sbjs: 230 train x 230 test
% same for non-cued
noncued_slice_vals = squeeze(mean(noncued_proj_all(PC_idx,:,:,:,:),[1,5]));      % mean over PC and sbjs: 230  x 230 x 4 distractors
        
        


figure
sgtitle('Temporal generalization: train time (x) by test time (y)', 'FontSize', 20)

tick_ms = -800:300:1400;  

for dist = 1:4
    subplot(2,2,dist)
    hold on
    imagesc(times, times, squeeze(noncued_slice_vals(:,:,dist)))
    axis xy
    xlim([times(1) times(end)])
    ylim([times(1) times(end)])
    axis square
    colormap(parula)
    colorbar
    clim([-0.02 0.02])
    
    % ---- time ticks ----
    xticks(tick_ms)
    yticks(tick_ms)
    xtickangle(0)
    
    xlabel('Train time (ms)')
    ylabel('Test time (ms)')
    
    % ---- reference lines ----
    xline(0,'k:','LineWidth',1.5)
    yline(0,'k:','LineWidth',1.5)
    plot([times(1) times(end)], [times(1) times(end)], ...
         'k--', 'LineWidth', 1)   % diagonal
    
    set(gca,'FontSize', 16)
    title(newtitles{dist}, 'FontSize', 18)
end



% mean across distractors
figure
sgtitle('Temporal generalization: distractor average', 'FontSize', 20)

hold on
imagesc(times, times, squeeze(mean(noncued_slice_vals,3)))
axis xy
xlim([times(1) times(end)])
ylim([times(1) times(end)])
axis square
colormap(parula)
colorbar
% clim([-0.005 0.005])

% ---- time ticks ----
xticks(tick_ms)
yticks(tick_ms)
xtickangle(0)

xlabel('Train time (ms)')
ylabel('Test time (ms)')

% ---- reference lines ----
xline(0,'k:','LineWidth',1.5)
yline(0,'k:','LineWidth',1.5)
plot([times(1) times(end)], [times(1) times(end)], ...
     'k--', 'LineWidth', 1)   % diagonal

set(gca,'FontSize', 16)





% cued projection 
figure
sgtitle('Temporal generalization: cued dimension', 'FontSize', 20)

hold on
imagesc(times, times, cued_slice_vals)
axis xy
xlim([times(1) times(end)])
ylim([times(1) times(end)])
axis square
colormap(parula)
colorbar
% clim([-0.005 0.005])

% ---- time ticks ----
xticks(tick_ms)
yticks(tick_ms)
xtickangle(0)

xlabel('Train time (ms)')
ylabel('Test time (ms)')

% ---- reference lines ----
xline(0,'k:','LineWidth',1.5)
yline(0,'k:','LineWidth',1.5)
plot([times(1) times(end)], [times(1) times(end)], ...
     'k--', 'LineWidth', 1)   % diagonal

set(gca,'FontSize', 16)






%% ── prevI minus prevC contrast: timeseries ────────────────────────────────
% Test whether previous trial state modulates current conflict-axis magnitude
% stratify current projections by previous congruency on each distractor
%  noncued_contrast_all: [4_pc × 230_train × 230_test × 4_dist × 26_subj]


% extract diagonal (train time == test time) and average over PCs
resh_c = reshape(noncued_contrast_all, 4, 230*230, 4, 26);  % [pc × t² × dist × subj]
contrast_diag_vals = resh_c(:, diag_lin_idx, :, :);         % [pc × 230 × dist × subj]
contrast_diag_vals = squeeze(mean(contrast_diag_vals(PC_idx,:,:,:), 1, 'omitnan')); % [230 × dist × subj]

mean_contrast = squeeze(mean(contrast_diag_vals, 3));        % [230 × 4]
sem_contrast  = squeeze(std(contrast_diag_vals,  0, 3)) / sqrt(size(contrast_diag_vals, 3)); % [230 × 4]

% dimension-average per subject first, then group mean/SEM
contrast_diag_avg_sn = squeeze(mean(contrast_diag_vals, 2)); % [230 × nSubj]
mean_contrast_avg    = mean(contrast_diag_avg_sn, 2);        % [230 × 1]
sem_contrast_avg     = std(contrast_diag_avg_sn, 0, 2) / sqrt(size(contrast_diag_avg_sn, 2));

shaded_area_x = [times, fliplr(times)];
dimLabels = {'C','S','E','P'};

% Figure 1: all 4 dimensions overlaid
figure('Color','w')
sgtitle('prevI – prevC contrast: timeseries (diagonal)', 'FontSize', 20)
hold on
for i = 1:4
    shaded_area_y = [mean_contrast(:,i)' - sem_contrast(:,i)', ...
                     fliplr(mean_contrast(:,i)' + sem_contrast(:,i)')];
    fill(shaded_area_x, shaded_area_y, colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(times, mean_contrast(:,i), 'LineWidth', 2, 'Color', colors(i,:));
end
plot(times, zeros(length(times),1), 'k:', 'LineWidth', 1.5)
xlim([times(1) times(end)])
xlabel('Time [ms]', 'FontSize', fontsize)
ylabel('prevI – prevC projection', 'FontSize', fontsize)
legend({'','C','','S','','E','','P'}, 'Location','southeast')
set(gca, 'FontSize', fontsize)
box off
hold off

% Figure 2: average only
figure('Color','w')
hold on
shaded_area_y = [mean_contrast_avg' - sem_contrast_avg', fliplr(mean_contrast_avg' + sem_contrast_avg')];
fill(shaded_area_x, shaded_area_y, [0.3 0.3 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(times, mean_contrast_avg, 'k-', 'LineWidth', 2.5);
plot(times, zeros(length(times),1), 'k:', 'LineWidth', 1.5)
xlim([times(1) times(end)])
xlabel('Time [ms]', 'FontSize', fontsize)
ylabel('prevI – prevC projection (avg)', 'FontSize', fontsize)
title('prevI – prevC contrast: dimension average', 'FontSize', fontsize)
set(gca, 'FontSize', fontsize)
box off
hold off





%% ── prevI minus prevC contrast: temporal generalization ───────────────────


% mean over PC_idx and subjects → [230_train × 230_test × 4_dist]
contrast_slice_vals = squeeze(mean(noncued_contrast_all(PC_idx,:,:,:,:), [1,5], 'omitnan'));

figure
sgtitle('Temporal generalization: prevI – prevC contrast (per distractor)', 'FontSize', 20)

for dist = 1:4
    subplot(2,2,dist)
    hold on
    imagesc(times, times, squeeze(contrast_slice_vals(:,:,dist)))
    axis xy
    xlim([times(1) times(end)])
    ylim([times(1) times(end)])
    axis square
    colormap(parula)
    colorbar
    % clim([-0.002 0.002])   % uncomment and tune if needed

    xticks(tick_ms); yticks(tick_ms); xtickangle(0)
    xlabel('Train time (ms)')
    ylabel('Test time (ms)')
    xline(0,'k:','LineWidth',1.5)
    yline(0,'k:','LineWidth',1.5)
    plot([times(1) times(end)], [times(1) times(end)], 'k--', 'LineWidth', 1)
    set(gca,'FontSize', 16)
    title(newtitles{dist}, 'FontSize', 18)
    hold off
end

% dimension-averaged contrast TG
figure
sgtitle('Temporal generalization: prevI – prevC contrast (distractor average)', 'FontSize', 20)
hold on
imagesc(times, times, squeeze(mean(contrast_slice_vals, 3, 'omitnan')))
axis xy
xlim([times(1) times(end)])
ylim([times(1) times(end)])
axis square
colormap(parula)
colorbar
% clim([-0.002 0.002])   % uncomment and tune if needed

xticks(tick_ms); yticks(tick_ms); xtickangle(0)
xlabel('Train time (ms)')
ylabel('Test time (ms)')
xline(0,'k:','LineWidth',1.5)
yline(0,'k:','LineWidth',1.5)
plot([times(1) times(end)], [times(1) times(end)], 'k--', 'LineWidth', 1)
set(gca,'FontSize', 16)
hold off


