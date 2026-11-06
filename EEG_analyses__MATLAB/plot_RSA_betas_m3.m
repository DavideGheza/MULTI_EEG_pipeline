%% Prepare and plot model-3 RSA conjunction time courses
% Loads participant fits for feature-specific conjunction coding in cued and
% non-cued dimensions. Produces overall, dimension-resolved, and
% proportion-congruency exports used in Figures 5 and 6.
% Author: Davide Gheza


%% 1. Load settings

clear

fits_sr = [];

fits_sr_dstr = [];

fits_sr_dstr_con = [];
fits_sr_dstr_inc = [];

fits_sr_singledstr = [];
fits_sr_singledstr_con = [];
fits_sr_singledstr_inc = [];

fits_sr_dstr_bytask = [];
fits_sr_bytask = [];

CLR = 1; % was classification profile transformed from [0 1] to centered log-ratio? 
% Stratify results by proportion-congruency condition.
stratify='PC';
% bootstrap (e.g., for controlling cng within each PC, match cng trials within PC, oversampling low cng trials)
bst = 1;
bootstrap_n = 10000;
% Time or tf domain
tf = true;

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
%If the folder doesn't exist yet, create one:
if ~exist(sets.rsa_data_dir,'dir')
    mkdir(sets.rsa_data_dir);
end


%%
sbjIDX = [22,24:48];
for sn = 1:length(sbjIDX)
    



if tf
    times = [-800:10:1498];
else
    times = [-1000:2:1999];
end




% pre-load single trial fits
load(fullfile(sets.rsa_data_dir, ['RSA_pc/multiclass/singletrial/fits/fits_singletrials_m3_sbj_' num2str(sbjIDX(sn)) '.mat'])) 
keep(sn).fits = fits; % all single trial data
% load behavioral data
behav = load([sets.sessions_data_dir 'V3_behav_seqinfo_trimmed_sbj_' num2str(sbjIDX(sn))]); % behavioral file with addidional info on subsequent and previous trials (10/01/2024)
behav = behav.behav_trim;

    if ~bst
    fits_sr = cat(3,fits_sr, squeeze(mean(fits.slope_sr,1,'omitnan')));% 1. compute mean over trials (trials not part of a session are NAN for that session) 2. catenate subjects over 3rd dim
    
    fits_sr_dstr = cat(3,fits_sr_dstr, squeeze(mean(fits.slope_sr_dstr,1,'omitnan')));
    
    fits_sr_dstr_con = cat(3,fits_sr_dstr_con, squeeze(mean(fits.slope_sr_dstr_con,1,'omitnan')));
    fits_sr_dstr_inc = cat(3,fits_sr_dstr_inc, squeeze(mean(fits.slope_sr_dstr_inc,1,'omitnan')));

    fits_sr_singledstr = cat(4,fits_sr_singledstr, squeeze(mean(fits.slope_sr_singledstr,1,'omitnan')));
    
    fits_sr_dstr_bytask_sn = nan(4,230,4,4); % preallocate sessions, times, non-cued dim, and cued task
    fits_sr_bytask_sn = nan(4,230,4); % preallocate sessions, times, and cued task
    for task=1:4
        fits_sr_dstr_bytask_sn(:,:,:,task) = squeeze(mean(fits.slope_sr_singledstr(behav.task==task,:,:,:),1,'omitnan'));
        fits_sr_bytask_sn(:,:,task) = squeeze(mean(fits.slope_sr(behav.task==task,:,:),1,'omitnan'));
    end
    fits_sr_dstr_bytask = cat(5,fits_sr_dstr_bytask, fits_sr_dstr_bytask_sn);
    fits_sr_bytask = cat(4,fits_sr_bytask, fits_sr_bytask_sn);
    
    else
    
    % bootstrap results to match congruency within PC (decimate and draw with replacement) # - define trial indices (03/28/2023-10/06/2023)
    % idxs = cong_control_ST_function(behav,bootstrap_n);
    % bootstrap results to match congruency within each cued task and PC (07/01/2024)
    idxs_singlecued = cong_control_ST_singlecued_function(behav,bootstrap_n);
    % bootstrap results to match congruency levels within distractor and PC (independently on each non-cued dimension - already quite balanced by default) (07/06/2024)
    [idxs, idxs_singledstr] = cong_control_ST_singledstr_function(behav,bootstrap_n);
    % preallocate sessions and times
    fits_sr_sn = nan(4,230);
    
    fits_sr_dstr_sn = nan(4,230);
    
    fits_sr_dstr_con_sn = nan(4,230);
    fits_sr_dstr_inc_sn = nan(4,230);

    fits_sr_singledstr_sn = nan(4,230,4);
    fits_sr_singledstr_con_sn = nan(4,230,4);
    fits_sr_singledstr_inc_sn = nan(4,230,4);
    fits_sr_dstr_bytask_sn = nan(4,230,4,4); % preallocate sessions, times, non-cued dim, and cued tasks ---> this is further splitting by cued tasks (4).
    fits_sr_bytask_sn = nan(4,230,4); % preallocate sessions, times, and cued task

        for s=1:4 % loop over session (PC) 
            fits_sr_sn(s,:) = squeeze(mean(fits.slope_sr(idxs(s).bst,s,:),1,'omitnan')); 
            
            fits_sr_dstr_sn(s,:) = squeeze(mean(fits.slope_sr_dstr(idxs(s).bst,s,:),1,'omitnan'));
            
            fits_sr_dstr_con_sn(s,:) = squeeze(mean(fits.slope_sr_dstr_con(idxs(s).bst,s,:),1,'omitnan'));
            fits_sr_dstr_inc_sn(s,:) = squeeze(mean(fits.slope_sr_dstr_inc(idxs(s).bst,s,:),1,'omitnan'));
            
            for dist=1:4
                fits_sr_singledstr_sn(s,:,dist) = squeeze(mean(fits.slope_sr_singledstr(idxs_singledstr(s,dist).bst,s,:,dist),1,'omitnan'));
                fits_sr_singledstr_con_sn(s,:,dist) = squeeze(mean(fits.slope_sr_singledstr(idxs_singledstr(s,dist).dist_con,s,:,dist),1,'omitnan'));
                fits_sr_singledstr_inc_sn(s,:,dist) = squeeze(mean(fits.slope_sr_singledstr(idxs_singledstr(s,dist).dist_inc,s,:,dist),1,'omitnan'));
            end
            for task=1:4
                fits_sr_dstr_bytask_sn(s,:,:,task) = squeeze(mean(fits.slope_sr_singledstr(idxs_singlecued(s,task).bst,s,:,:),1,'omitnan'));
                fits_sr_bytask_sn(s,:,task) = squeeze(mean(fits.slope_sr(idxs_singlecued(s,task).bst,s,:),1,'omitnan'));
            end
        end
        fits_sr = cat(3,fits_sr, fits_sr_sn);
        
        fits_sr_dstr = cat(3,fits_sr_dstr, fits_sr_dstr_sn);
        
        fits_sr_dstr_con = cat(3,fits_sr_dstr_con, fits_sr_dstr_con_sn);
        fits_sr_dstr_inc = cat(3,fits_sr_dstr_inc, fits_sr_dstr_inc_sn);

        fits_sr_singledstr = cat(4,fits_sr_singledstr, fits_sr_singledstr_sn);
        fits_sr_singledstr_con = cat(4,fits_sr_singledstr_con, fits_sr_singledstr_con_sn);
        fits_sr_singledstr_inc = cat(4,fits_sr_singledstr_inc, fits_sr_singledstr_inc_sn);

        fits_sr_dstr_bytask = cat(5,fits_sr_dstr_bytask, fits_sr_dstr_bytask_sn);
        fits_sr_bytask = cat(4,fits_sr_bytask, fits_sr_bytask_sn);
    
    end



fprintf('\tLoaded sbj %g\n', sn) 
end






% plots for singletrials

% permute 
fits_sr = permute(fits_sr,[3,1,2]);

fits_sr_dstr = permute(fits_sr_dstr,[3,1,2]);

fits_sr_dstr_con = permute(fits_sr_dstr_con,[3,1,2]);
fits_sr_dstr_inc = permute(fits_sr_dstr_inc,[3,1,2]);

fits_sr_singledstr = permute(fits_sr_singledstr,[4,1,2,3]);
fits_sr_singledstr_con = permute(fits_sr_singledstr_con,[4,1,2,3]);
fits_sr_singledstr_inc = permute(fits_sr_singledstr_inc,[4,1,2,3]);

fits_sr_dstr_bytask = permute(fits_sr_dstr_bytask,[5,1,2,3,4]);
fits_sr_bytask = permute(fits_sr_bytask,[4,1,2,3]);


titles={'Neu1', 'MI', 'MC', 'Neu2', 'MI-MC', 'MI-NeuAVG', 'MC-NeuAVG'};
contrasts = {2,3;2,[1,4];3,[1,4]};
sesn = 4;


    
% reported 2025                                                 (05/22/25)
% store and export to R, for plotting timeseries - assumes stratify='PC'; 
% single-trial fits are balanced with cong_control function (balance cng within PC), and then averaged within PC 
% fair to average timeseries across PCs
avg.fits_sr = fits_sr;
avg.fits_sr_dstr = fits_sr_dstr;
avg.fits_sr_dstr_con = fits_sr_dstr_con; % added yet not used/reported
avg.fits_sr_dstr_inc = fits_sr_dstr_inc; % added yet not used/reported

% save Fig5_Fig6_m3_RSA_betas.mat avg                % run with bst=1

% single-trial fits stratified by dimension (cued, or distractor) 
% since intent of this plot is to show RSA fit by dimension across PCs, bootstrapping is unnecessary and not reccomended.
avg2.fits_sr_bytask = fits_sr_bytask;                    % 4d: sbjs, ses, time, task
avg2.fits_sr_singledstr = fits_sr_singledstr;            % 4d: sbjs, ses, time, dist
avg2.fits_sr_dstr_bytask = fits_sr_dstr_bytask;          % 5d: sbjs, ses, time, dist, task

% save Fig5_m3_RSA_betas_bydim_bst0.mat avg2         % run with bst=0

% Same as before, but WITH bootstrapping for PC effects (single-trial fits stratified by dimension (cued, or distractor) 

% save Fig6_m3_RSA_betas_bydim.mat avg2         % run with bst=1


        %% cued dim and distractors' conjunction

        colors = num2cell([0 0 1; 0 0 0],2);
        
        if CLR
            ylims = [-0.5 0.8]; 
        else
            ylims = [-0.03 0.03]; 
        end
        
        figure
        sgtitle('RSA models fits for cued dimensions', fontsize=20)
        for ses = 1:sesn
        subplot(2,4,ses)
        h = plot(times,squeeze([mean(fits_sr(:,ses,:),1); mean(fits_sr_dstr(:,ses,:),1) ])', 'linewidth', 2);
        set(h, {'color'}, colors);
        hold on
        plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
        ylim(ylims)
        xlim([times(1) times(end)])
        xlabel('Time [ms]', 'FontSize', 18)
        ylabel('RSA beta', 'FontSize', 18)
        set(gca,'FontSize', 18)
        legend('cued SR', 'non-cued SR', 'location', 'south'), 
        title(titles{ses}, 'FontSize', 18);
        hold off
        end
            subplot(2,4,5)
            h = plot(times,squeeze([(mean(fits_sr(:,contrasts{1,1},:),1)-mean(fits_sr(:,contrasts{1,2},:),1)), (mean(fits_sr_dstr(:,contrasts{1,1},:),1)-mean(fits_sr_dstr(:,contrasts{1,2},:),1))  ]), 'linewidth', 2); 
            set(h, {'color'}, colors);
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]')
            ylabel('RSA beta - MI-MC')
            legend('cued SR', 'non-cued SR', 'location', 'south'), 
            title(titles(5), 'FontSize', 18);
                subplot(2,4,6)
            h = plot(times,squeeze([(mean(fits_sr(:,contrasts{2,1},:),1)-mean(fits_sr(:,contrasts{2,2},:),1:2)), (mean(fits_sr_dstr(:,contrasts{2,1},:),1)-mean(fits_sr_dstr(:,contrasts{2,2},:),1:2))  ]), 'linewidth', 2); 
            set(h, {'color'}, colors);
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MI-NeuAVG', 'FontSize', 18)
            legend('cued SR', 'non-cued SR', 'location', 'south'), 
            title(titles(6), 'FontSize', 18);
                subplot(2,4,7)
            h = plot(times,squeeze([(mean(fits_sr(:,contrasts{3,1},:),1)-mean(fits_sr(:,contrasts{3,2},:),1:2)),  (mean(fits_sr_dstr(:,contrasts{3,1},:),1)-mean(fits_sr_dstr(:,contrasts{3,2},:),1:2))  ]), 'linewidth', 2);
            set(h, {'color'}, colors);
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MC-NeuAVG', 'FontSize', 18)
            legend('cued SR', 'non-cued SR', 'location', 'south'), 
            title(titles(7), 'FontSize', 18);
        
    

        RSA_permutest(fits_sr_dstr,times)    

        RSA_permutest(fits_sr,times)    


        %% single cued dimensions' conjunctions (symmetrical to single distractors)
        % Note: balanced cng by cued taks and PC. If bst==1, some missing conditions caues some cng levels by ses by task have 0 trial count (see below)

        colors = num2cell([autumn(4)],2);
        ylims = [-1 1.3]; 
        
        figure
        sgtitle('RSA models fits for cued dimensions by cued task', fontsize=20)
        for ses = 1:sesn
        subplot(2,4,ses)
        h = plot(times,squeeze([mean(fits_sr_bytask(:,ses,:,1),1); mean(fits_sr_bytask(:,ses,:,2),1); mean(fits_sr_bytask(:,ses,:,3),1);  mean(fits_sr_bytask(:,ses,:,4),1)])', 'linewidth', 2);
        set(h, {'color'}, colors);
        hold on
        plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
        ylim(ylims)
        xlim([times(1) times(end)])
        xlabel('Time [ms]', 'FontSize', 18)
        ylabel('RSA beta', 'FontSize', 18)
        set(gca,'FontSize', 18)
        legend('C SR', 'S SR', 'E SR', 'P SR', 'location', 'southwest', 'FontSize', 18), 
        title(titles{ses}, 'FontSize', 18);
        hold off
        end
        
            
            
            subplot(2,4,5)
            h = plot(times,squeeze([(mean(fits_sr_bytask(:,contrasts{1,1},:,1),1)-mean(fits_sr_bytask(:,contrasts{1,2},:,1),1)), (mean(fits_sr_bytask(:,contrasts{1,1},:,2),1)-mean(fits_sr_bytask(:,contrasts{1,2},:,2),1)), (mean(fits_sr_bytask(:,contrasts{1,1},:,3),1)-mean(fits_sr_bytask(:,contrasts{1,2},:,3),1)), (mean(fits_sr_bytask(:,contrasts{1,1},:,4),1)-mean(fits_sr_bytask(:,contrasts{1,2},:,4),1))   ]), 'linewidth', 2);
            set(h, {'color'}, colors);
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MI-MC', 'FontSize', 18)
            set(gca,'FontSize', 18)
            title(titles(5), 'FontSize', 18);
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            hold off
            legend('C SR', 'S SR', 'E SR', 'P SR', 'location', 'southwest', 'FontSize', 18), 
                subplot(2,4,6)
            h = plot(times,squeeze([(mean(fits_sr_bytask(:,contrasts{2,1},:,1),1)-mean(fits_sr_bytask(:,contrasts{2,2},:,1),1:2)), (mean(fits_sr_bytask(:,contrasts{2,1},:,2),1)-mean(fits_sr_bytask(:,contrasts{2,2},:,2),1:2)), (mean(fits_sr_bytask(:,contrasts{2,1},:,3),1)-mean(fits_sr_bytask(:,contrasts{2,2},:,3),1:2)), (mean(fits_sr_bytask(:,contrasts{2,1},:,4),1)-mean(fits_sr_bytask(:,contrasts{2,2},:,4),1:2))   ]), 'linewidth', 2); 
            set(h, {'color'}, colors);
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MI-NeuAVG', 'FontSize', 18)
            set(gca,'FontSize', 18)
            title(titles(6), 'FontSize', 18);
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            hold off
            legend('C SR', 'S SR', 'E SR', 'P SR', 'location', 'southwest', 'FontSize', 18), 
                subplot(2,4,7)
            h = plot(times,squeeze([(mean(fits_sr_bytask(:,contrasts{3,1},:,1),1)-mean(fits_sr_bytask(:,contrasts{3,2},:,1),1:2)), (mean(fits_sr_bytask(:,contrasts{3,1},:,2),1)-mean(fits_sr_bytask(:,contrasts{3,2},:,2),1:2)), (mean(fits_sr_bytask(:,contrasts{3,1},:,3),1)-mean(fits_sr_bytask(:,contrasts{3,2},:,3),1:2)), (mean(fits_sr_bytask(:,contrasts{3,1},:,4),1)-mean(fits_sr_bytask(:,contrasts{3,2},:,4),1:2))   ]), 'linewidth', 2);
            set(h, {'color'}, colors);
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MC-NeuAVG', 'FontSize', 18)
            set(gca,'FontSize', 18)
            title(titles(7), 'FontSize', 18);
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            hold off
            legend('C SR', 'S SR', 'E SR', 'P SR', 'location', 'southwest', 'FontSize', 18), 







        %% single distractors - Note: cng is balanced for each distractor (but already ~50% by design)
        
        colors = num2cell([autumn(4)],2);
        
        if CLR
            ylims = [-0.3 0.3]; 
        else
            ylims = [-0.03 0.03]; 
        end
        
        figure
        sgtitle('RSA models fits for non-cued dimensions', fontsize=20)
        for ses = 1:sesn
        subplot(2,4,ses)
        h = plot(times,squeeze([mean(fits_sr_singledstr(:,ses,:,1),1); mean(fits_sr_singledstr(:,ses,:,2),1); mean(fits_sr_singledstr(:,ses,:,3),1);  mean(fits_sr_singledstr(:,ses,:,4),1)])', 'linewidth', 2);
        set(h, {'color'}, colors);
        hold on
        plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
        ylim(ylims)
        xlim([times(1) times(end)])
        xlabel('Time [ms]', 'FontSize', 18)
        ylabel('RSA beta', 'FontSize', 18)
        set(gca,'FontSize', 18)
        % legend('C non-cued conjunction', 'S non-cued conjunction', 'E non-cued conjunction', 'P non-cued conjunction', 'location', 'north', 'FontSize', 18), 
        title(titles{ses}, 'FontSize', 18);
        hold off
        end
        
            
            
            subplot(2,4,5)
            h = plot(times,squeeze([(mean(fits_sr_singledstr(:,contrasts{1,1},:,1),1)-mean(fits_sr_singledstr(:,contrasts{1,2},:,1),1)), (mean(fits_sr_singledstr(:,contrasts{1,1},:,2),1)-mean(fits_sr_singledstr(:,contrasts{1,2},:,2),1)), (mean(fits_sr_singledstr(:,contrasts{1,1},:,3),1)-mean(fits_sr_singledstr(:,contrasts{1,2},:,3),1)), (mean(fits_sr_singledstr(:,contrasts{1,1},:,4),1)-mean(fits_sr_singledstr(:,contrasts{1,2},:,4),1))   ]), 'linewidth', 2);
            set(h, {'color'}, colors);
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MI-MC', 'FontSize', 18)
            set(gca,'FontSize', 18)
            title(titles(5), 'FontSize', 18);
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            hold off
            % legend('C non-cued conjunction', 'S non-cued conjunction', 'E non-cued conjunction', 'P non-cued conjunction', 'location', 'southwest', 'FontSize', 18), 
                subplot(2,4,6)
            h = plot(times,squeeze([(mean(fits_sr_singledstr(:,contrasts{2,1},:,1),1)-mean(fits_sr_singledstr(:,contrasts{2,2},:,1),1:2)), (mean(fits_sr_singledstr(:,contrasts{2,1},:,2),1)-mean(fits_sr_singledstr(:,contrasts{2,2},:,2),1:2)), (mean(fits_sr_singledstr(:,contrasts{2,1},:,3),1)-mean(fits_sr_singledstr(:,contrasts{2,2},:,3),1:2)), (mean(fits_sr_singledstr(:,contrasts{2,1},:,4),1)-mean(fits_sr_singledstr(:,contrasts{2,2},:,4),1:2))   ]), 'linewidth', 2); 
            set(h, {'color'}, colors);
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MI-NeuAVG', 'FontSize', 18)
            set(gca,'FontSize', 18)
            title(titles(6), 'FontSize', 18);
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            hold off
            % legend('C non-cued conjunction', 'S non-cued conjunction', 'E non-cued conjunction', 'P non-cued conjunction', 'location', 'southwest', 'FontSize', 18), 
                subplot(2,4,7)
            h = plot(times,squeeze([(mean(fits_sr_singledstr(:,contrasts{3,1},:,1),1)-mean(fits_sr_singledstr(:,contrasts{3,2},:,1),1:2)), (mean(fits_sr_singledstr(:,contrasts{3,1},:,2),1)-mean(fits_sr_singledstr(:,contrasts{3,2},:,2),1:2)), (mean(fits_sr_singledstr(:,contrasts{3,1},:,3),1)-mean(fits_sr_singledstr(:,contrasts{3,2},:,3),1:2)), (mean(fits_sr_singledstr(:,contrasts{3,1},:,4),1)-mean(fits_sr_singledstr(:,contrasts{3,2},:,4),1:2))   ]), 'linewidth', 2);
            set(h, {'color'}, colors);
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MC-NeuAVG', 'FontSize', 18)
            set(gca,'FontSize', 18)
            title(titles(7), 'FontSize', 18);
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            hold off
            legend('C non-cued conjunction', 'S non-cued conjunction', 'E non-cued conjunction', 'P non-cued conjunction', 'location', 'southwest', 'FontSize', 18), 
        
    
        RSA_permutest(mean(fits_sr_singledstr,4),times,stratify)   
        RSA_permutest(fits_sr_singledstr(:,:,:,1),times,stratify)    
        RSA_permutest(fits_sr_singledstr(:,:,:,2),times,stratify)  
        RSA_permutest(fits_sr_singledstr(:,:,:,3),times,stratify)   
        RSA_permutest(fits_sr_singledstr(:,:,:,4),times,stratify)    


        %% distractor congruent vs incongruent

        colors = num2cell([autumn(2)],2);
        
        if CLR
            ylims = [-1 1]; 
        else
            ylims = [-0.03 0.03]; 
        end
        
        figure
        sgtitle('RSA models fits for non-cued dimensions - when congruent vs incongruent', fontsize=20)
        for ses = 1:sesn
        subplot(2,4,ses)
        h = plot(times,squeeze([mean(fits_sr_dstr_con(:,ses,:),1); mean(fits_sr_dstr_inc(:,ses,:),1)])', 'linewidth', 2);
        set(h, {'color'}, colors);
        hold on
        plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
        ylim(ylims)
        xlim([times(1) times(end)])
        xlabel('Time [ms]', 'FontSize', 18)
        ylabel('RSA beta', 'FontSize', 18)
        set(gca,'FontSize', 18)
        legend('congruent non-cued conjunctions', 'incongruent non-cued conjunctions', 'location', 'southwest', 'FontSize', 18), 
        title(titles{ses}, 'FontSize', 18);
        hold off
        end
        
            subplot(2,4,5)
            plot(times,squeeze([(mean(fits_sr_dstr_con(:,contrasts{1,1},:),1)-mean(fits_sr_dstr_con(:,contrasts{1,2},:),1)), (mean(fits_sr_dstr_inc(:,contrasts{1,1},:),1)-mean(fits_sr_dstr_inc(:,contrasts{1,2},:),1))  ]), 'linewidth', 2) 
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MI-MC', 'FontSize', 18)
            set(gca,'FontSize', 18)
            legend('congruent non-cued conjunctions', 'incongruent non-cued conjunctions', 'location', 'southwest', 'FontSize', 18),
            title(titles(5), 'FontSize', 18);
                subplot(2,4,6)
            plot(times,squeeze([(mean(fits_sr_dstr_con(:,contrasts{2,1},:),1)-mean(fits_sr_dstr_con(:,contrasts{2,2},:),1:2)), (mean(fits_sr_dstr_inc(:,contrasts{2,1},:),1)-mean(fits_sr_dstr_inc(:,contrasts{2,2},:),1:2))  ]), 'linewidth', 2) 
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MI-NeuAVG', 'FontSize', 18)
            set(gca,'FontSize', 18)
            legend('congruent non-cued conjunctions', 'incongruent non-cued conjunctions', 'location', 'southwest', 'FontSize', 18),
            title(titles(6), 'FontSize', 18);
                subplot(2,4,7)
            plot(times,squeeze([(mean(fits_sr_dstr_con(:,contrasts{3,1},:),1)-mean(fits_sr_dstr_con(:,contrasts{3,2},:),1:2)), (mean(fits_sr_dstr_inc(:,contrasts{3,1},:),1)-mean(fits_sr_dstr_inc(:,contrasts{3,2},:),1:2))  ]), 'linewidth', 2)
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MC-NeuAVG', 'FontSize', 18)
            set(gca,'FontSize', 18)
            legend('congruent non-cued conjunctions', 'incongruent non-cued conjunctions', 'location', 'southwest', 'FontSize', 18),
            title(titles(7), 'FontSize', 18);
        
    
        RSA_permutest(fits_sr_dstr_con,times)    

       
    
 %% single distractors by cued task. Note: balanced cng by cued taks and PC.
        
        figure

        colors = num2cell([summer(4)],2);
        newtitles = {'C Neu1', 'C MI', 'C MC', 'C Neu2', ... % letter refers to cued task
                     'S Neu1', 'S MI', 'S MC', 'S Neu2', ...
                     'E Neu1', 'E MI', 'E MC', 'E Neu2', ...
                     'P Neu1', 'P MI', 'P MC', 'P Neu2', ...
                    };  

        if CLR
            ylims = [-0.3 0.3]; 
        else
            ylims = [-0.03 0.03]; 
        end
        
        sgtitle('RSA models fits for non-cued dimensions - single distractor by cued task', fontsize=20)
        for cued = 1:4
            for ses = 1:sesn
            subplot(4,4,ses+4*(cued-1))
            h = plot(times,squeeze([mean(fits_sr_dstr_bytask(:,ses,:,:,cued),1)  ])', 'linewidth', 2); % sbjs, ses, time, dist, task
            set(h, {'color'}, colors);
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta', 'FontSize', 18)
            set(gca,'FontSize', 18)
            if cued == 1 && ses == 1
                legend('dist C', 'dist S', 'dist E', 'dist P', 'location', 'southwest', 'FontSize', 18), 
            end
            title(newtitles{ses+4*(cued-1)}, 'FontSize', 18);
            hold off
            end
        end







        figure

        colors = num2cell([autumn(4)],2);
        newtitles = {'C Neu1', 'C MI', 'C MC', 'C Neu2', ... % letter refers to distractor
                     'S Neu1', 'S MI', 'S MC', 'S Neu2', ...
                     'E Neu1', 'E MI', 'E MC', 'E Neu2', ...
                     'P Neu1', 'P MI', 'P MC', 'P Neu2', ...
                    };  


        if CLR
            ylims = [-0.3 0.3]; 
        else
            ylims = [-0.03 0.03]; 
        end
        
        sgtitle('RSA models fits for non-cued dimensions - single distractor by cued task', fontsize=20)
        for dist = 1:4
            for ses = 1:sesn
            subplot(4,4,ses+4*(dist-1))
            h = plot(times,squeeze([mean(fits_sr_dstr_bytask(:,ses,:,dist,:),1)  ])', 'linewidth', 2); % sbjs, ses, time, dist, task
            set(h, {'color'}, colors);
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta', 'FontSize', 18)
            set(gca,'FontSize', 18)
            if dist == 1 && ses == 1
                legend('cued C', 'cued S', 'cued E', 'cued P', 'location', 'southwest', 'FontSize', 18), 
            end
            title(newtitles{ses+4*(dist-1)}, 'FontSize', 18);
            hold off
            end
        end



 %% congruency effect by single distractors      

colors = num2cell([autumn(2); 0 0 0],2);% Add black as the third color
newtitles = {'C Neu1', 'C MI', 'C MC', 'C Neu2', ...
             'S Neu1', 'S MI', 'S MC', 'S Neu2', ...
             'E Neu1', 'E MI', 'E MC', 'E Neu2', ...
             'P Neu1', 'P MI', 'P MC', 'P Neu2', ...
            };  

        if CLR
            ylims = [-0.5 1.0]; 
        else
            ylims = [-0.03 0.03]; 
        end
        
        figure
        sgtitle('RSA models fits for non-cued dimensions - when congruent vs incongruent', fontsize=20)
        for dist = 1:4
            for ses = 1:sesn
            subplot(4,4,ses+4*(dist-1))
            h = plot(times,squeeze([mean(fits_sr_singledstr_con(:,ses,:,dist),1); ...% sbjs, ses, time, dist
                                    mean(fits_sr_singledstr_inc(:,ses,:,dist),1); ...
                                    (mean(fits_sr_singledstr_con(:,ses,:,dist),1) - mean(fits_sr_singledstr_inc(:,ses,:,dist),1)) ])', 'linewidth', 2);
            set(h, {'color'}, colors);
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta', 'FontSize', 18)
            set(gca,'FontSize', 18)
            if dist == 1 && ses == 1
                legend('congruent non-cued conjunctions', 'incongruent non-cued conjunctions', 'con-inc', 'location', 'southwest', 'FontSize', 18), 
            end
            title(newtitles{ses+4*(dist-1)}, 'FontSize', 18);
            hold off
            end
        end

        


