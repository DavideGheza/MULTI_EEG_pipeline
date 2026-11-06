%%  Dimensions EEG decoding SIDE and task-ID - V3 
% Load classification based distances - plot and export accuracy timeseries 
% (03/06/2025)
% based on script for diagnostics on variance structure (06/11,14/2024)


%% 1. Load settings

clear

% Time or tf domain
tf = true;

% choose plot type: 1=accuracy over time, 2=prob. over time, 3=prob.structure
pt = 1;


%Struct 'sets' is loaded from the 'scripts_and_settings' folder:
if ispc
    load('C:\Users\gdavide\OneDrive\Documenti\MATLAB\dg_workspace22\dimensions\analyses\scripts_and_settings\preprocessing_settings.mat');
else
    load('/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/scripts_and_settings/preprocessing_settings_mac.mat');
end


%Add path to sub-folder:
sets.rsa_data_dir=[sets.main_saving_dir 'rsa_data/RSA_pc/multiclass/singletrial/'];



%%
sbjIDX = [22,24:48];


for sn = 1:length(sbjIDX)
    
rng(10);


map = [1,5,2,6,3,7,4,8]; % order of LDA results with respect to RSA models


if tf
    times = [-800:10:1498];
else
    times = [-1000:2:1999];
end

% pre-load distances
load(fullfile(sets.rsa_data_dir, ['distances_mltclass_repCV_effHPT_bal_sbj_' num2str(sbjIDX(sn)) '.mat'])) % local renamed copy (see README)

if pt==1
    figure;
    % Look at the results
    hold on
    plot(times, 100*squeeze(nanmean(result.lda, 1))+50, 'linewidth', 2)
    plot([times(1) times(end)], [12.5 12.5], 'k-')
    xlim([times(1) times(end)])
    xlabel('Time [ms]')
    ylabel('Classification accuracy')
    legend('LDA')
    hold off

    figure;
    % Look at the results by session
    hold on
    plot(times, 100*squeeze(result.lda(1,:,:,:))+50, 'linewidth', 2)
    plot(times, 100*squeeze(result.lda(2,:,:,:))+50, 'linewidth', 2)
    plot(times, 100*squeeze(result.lda(3,:,:,:))+50, 'linewidth', 2)
    plot(times, 100*squeeze(result.lda(4,:,:,:))+50, 'linewidth', 2)
    plot([times(1) times(end)], [12.5 12.5], 'k-')
    xlim([times(1) times(end)])
    xlabel('Time [ms]')
    ylabel('Classification accuracy (LDA) - over sessions')
    legend('ses1', 'ses2', 'ses3', 'ses4')
    hold off

elseif pt==2    
    % plot probability over diagonal (for correct class) - mean over classes and PC/CNG
    figure
    hold on
    [n_ses,n_cond,n_cond2,n_times] = size(result.lda_prob_full);
    tmp = reshape(result.lda_prob_full,[n_ses,n_cond*n_cond2,n_times]);
    linIndsDiag = sub2ind([n_cond,n_cond],1:n_cond,1:8);
    diag = tmp(:,linIndsDiag,:);
    plot(times, 100*squeeze(nanmean(nanmean(diag, 1), 2)), 'linewidth', 2)
    plot([times(1) times(end)], [12.5 12.5], 'k-')
    xlim([times(1) times(end)])
    xlabel('Time [ms]')
    ylabel('Classification accuracy')
    legend('LDA_probability_for_correct')
    hold off




elseif pt==3
    figure;
    % Look at the variance structure in classification probability vectors
    timewin = [450,550];
    timewin = [100,700];
    timesidxs = [dsearchn(times',timewin(1)), dsearchn(times',timewin(2))];
    for s=1:4
        subplot(2,2,s)
        full_prob_matrix = squeeze(mean(result.lda_prob_full(s,map,map,(timesidxs(1):timesidxs(2))),4, 'omitnan')); % mean over time window
        group_full_prob_matrix(sn,s,:,:) = full_prob_matrix;
        h = heatmap(full_prob_matrix);
        set(gca,'fontsize',14)
    %     colormap turbo
        ylabel('true class') % (average over test trials)
        xlabel('tested class')
        title({'LDA prob - sbj ' num2str(sbjIDX(sn)) })
        ax = struct(gca);
        ylabel(ax.Colorbar, 'mean p');
    end
    % "take CRL prob, average them by class. Each trial has vector of 8, average trials based on class. Fill out 8x8 matrix with these results. "

elseif isnan(pt)
end


keep.group_lda(sn,:,:) = result.lda;
keep.group_lda_prob_full(sn,:,:,:,:) = result.lda_prob_full;
keep.group_lda_cond_acc(sn,:,:,:) = result.lda_cond_acc;





end



% export to R 
save Fig4_classification_accuracy.mat keep


    
% plot group-average accuracy by condition
figure;
colors = winter(4);
% Look at the results by session
hold on
for co = 1:size(keep.group_lda_cond_acc,3)/2
    plot(times, squeeze(mean(100*keep.group_lda_cond_acc(:,:,co,:),(1:2))), 'linewidth', 2, 'Color', colors(co,:))
end
for co = 5:size(keep.group_lda_cond_acc,3)
    plot(times, squeeze(mean(100*keep.group_lda_cond_acc(:,:,co,:),(1:2))), 'linewidth', 2, 'LineStyle','--', 'Color', colors(co-4,:))
end
plot([times(1) times(end)], [12.5 12.5], 'k-')
xlim([times(1) times(end)])
xlabel('Time [ms]')
ylabel('accuracy (%)')
title('Group-average accuracy (LDA) - by class')
legend('Left Color', 'Left Shape', 'Left Edge', 'Left Pattern', 'Right Color', 'Right Shape', 'Right Edge', 'Right Pattern')
hold off
ylim([5 40])


if pt==3
    % plot group-average variance structure in classification probability vectors
    figure
    for s=1:4
        subplot(2,2,s)
        matrix = mean(group_full_prob_matrix(:,s,:,:),1); % group mean, by ses
        h = heatmap(squeeze(matrix));
        set(gca,'fontsize',14)
    %     colormap turbo
        ylabel('true class') % (average over test trials)
        xlabel('tested class')
        title({'mean LDA prob'})
        ax = struct(gca);
        ylabel(ax.Colorbar, 'p');
    end
end








