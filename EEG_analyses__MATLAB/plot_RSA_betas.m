%% Prepare and plot model-1 RSA time courses
% Loads participant fits, aggregates cued model-1 RSA estimates across
% proportion-congruency conditions, and exports the arrays used for Figure 4.
% The script also contains optional CSE stratification.
% Author: Davide Gheza


%% 1. Load settings

clear

fits_sr = [];
fits_ru = [];
fits_si = [];

fits_sr_dstr = [];
fits_ru_dstr = [];
fits_si_dstr = [];

fits_sr_dstr_con = [];
fits_ru_dstr_con = [];
fits_si_dstr_con = [];
fits_sr_dstr_inc = [];
fits_ru_dstr_inc = [];
fits_si_dstr_inc = [];

fits_sr_singledstr = [];
fits_ru_singledstr = [];
fits_si_singledstr = [];
fits_sr_bytask = [];
fits_ru_bytask = [];
fits_si_bytask = [];

    
CLR = 1; % was classification profile transformed from [0 1] to centered log-ratio? 
% cued and non-cued are fitted with two sequential models. Non-cued fitted on residuals of fittig on cued?
residuals = 0;
% Stratify results by proportion-congruency condition.
stratify='PC';
% bootstrap (e.g., for controlling cng within each PC, match cng trials within PC, oversampling low cng trials)
bst = 1;
bootstrap_n = 1000;
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
if residuals == 1
    load(fullfile(sets.rsa_data_dir, ['RSA_pc/multiclass/singletrial/fits/fits_singletrials_residualsNC_sbj_' num2str(sbjIDX(sn)) '.mat'])) % filename as from HPC
else 
    % load(fullfile(sets.rsa_data_dir, ['RSA_pc/multiclass/singletrial/fits/fits_singletrials_m1_sbj_' num2str(sbjIDX(sn)) '.mat'])) % filename as from HPC
    load(fullfile(sets.rsa_data_dir, ['RSA_pc/multiclass/singletrial/fits/fits_singletrials_m1_RTnsc_sbj_' num2str(sbjIDX(sn)) '.mat'])) % fit with RT nuisance predictor (06/15/2025)
end
keep(sn).fits = fits; % all single trial data
% load behavioral data
behav = load([sets.sessions_data_dir 'V3_behav_seqinfo_trimmed_sbj_' num2str(sbjIDX(sn))]); % behavioral file with addidional info on subsequent and previous trials (10/01/2024)
behav = behav.behav_trim;

    if ~bst
    fits_sr = cat(3,fits_sr, squeeze(mean(fits.slope_sr,1,'omitnan')));% 1. compute mean over trials (trials not part of a session are NAN for that session) 2. catenate subjects over 3rd dim
    fits_ru = cat(3,fits_ru, squeeze(mean(fits.slope_ru,1,'omitnan')));
    fits_si = cat(3,fits_si, squeeze(mean(fits.slope_si,1,'omitnan')));
    
    fits_sr_dstr = cat(3,fits_sr_dstr, squeeze(mean(fits.slope_sr_dstr,1,'omitnan')));
    fits_ru_dstr = cat(3,fits_ru_dstr, squeeze(mean(fits.slope_ru_dstr,1,'omitnan')));
    fits_si_dstr = cat(3,fits_si_dstr, squeeze(mean(fits.slope_si_dstr,1,'omitnan')));
    
    fits_sr_dstr_con = cat(3,fits_sr_dstr_con, squeeze(mean(fits.slope_sr_dstr_con,1,'omitnan')));
    fits_ru_dstr_con = cat(3,fits_ru_dstr_con, squeeze(mean(fits.slope_ru_dstr_con,1,'omitnan')));
    fits_si_dstr_con = cat(3,fits_si_dstr_con, squeeze(mean(fits.slope_si_dstr_con,1,'omitnan')));
    fits_sr_dstr_inc = cat(3,fits_sr_dstr_inc, squeeze(mean(fits.slope_sr_dstr_inc,1,'omitnan')));
    fits_ru_dstr_inc = cat(3,fits_ru_dstr_inc, squeeze(mean(fits.slope_ru_dstr_inc,1,'omitnan')));
    fits_si_dstr_inc = cat(3,fits_si_dstr_inc, squeeze(mean(fits.slope_si_dstr_inc,1,'omitnan')));

    % added following m3 (03/20/2025)
    fits_sr_singledstr = cat(4,fits_sr_singledstr, squeeze(mean(fits.slope_sr_singledstr,1,'omitnan')));
    fits_ru_singledstr = cat(4,fits_ru_singledstr, squeeze(mean(fits.slope_ru_singledstr,1,'omitnan')));
    fits_si_singledstr = cat(4,fits_si_singledstr, squeeze(mean(fits.slope_si_singledstr,1,'omitnan')));
    fits_sr_bytask_sn = nan(4,230,4); % preallocate sessions, times, and cued task
    fits_ru_bytask_sn = nan(4,230,4);
    fits_si_bytask_sn = nan(4,230,4);
    for task=1:4
        fits_sr_bytask_sn(:,:,task) = squeeze(mean(fits.slope_sr(behav.task==task,:,:),1,'omitnan'));
        fits_ru_bytask_sn(:,:,task) = squeeze(mean(fits.slope_ru(behav.task==task,:,:),1,'omitnan'));
        fits_si_bytask_sn(:,:,task) = squeeze(mean(fits.slope_si(behav.task==task,:,:),1,'omitnan'));
    end
    fits_sr_bytask = cat(4,fits_sr_bytask, fits_sr_bytask_sn);
    fits_ru_bytask = cat(4,fits_ru_bytask, fits_ru_bytask_sn);
    fits_si_bytask = cat(4,fits_si_bytask, fits_si_bytask_sn);

    else
    
    % bootstrap results to match congruency (decimate and draw with replacement) # - define trial indices (03/28/2023-10/06/2023)
    idxs = cong_control_ST_function(behav,bootstrap_n);
    % preallocate sessions and times
    fits_sr_sn = nan(4,230);
    fits_ru_sn = nan(4,230);
    fits_si_sn = nan(4,230);
    
    fits_sr_dstr_sn = nan(4,230);
    fits_ru_dstr_sn = nan(4,230);
    fits_si_dstr_sn = nan(4,230);
    
    fits_sr_dstr_con_sn = nan(4,230);
    fits_ru_dstr_con_sn = nan(4,230);
    fits_si_dstr_con_sn = nan(4,230);
    fits_sr_dstr_inc_sn = nan(4,230);
    fits_ru_dstr_inc_sn = nan(4,230);
    fits_si_dstr_inc_sn = nan(4,230);
    
        for s=1:4 % loop over session (PC) 
            fits_sr_sn(s,:) = squeeze(mean(fits.slope_sr(idxs(s).bst,s,:),1,'omitnan')); % a few nan fits, due to not convergence of lm?
            fits_ru_sn(s,:) = squeeze(mean(fits.slope_ru(idxs(s).bst,s,:),1,'omitnan'));
            fits_si_sn(s,:) = squeeze(mean(fits.slope_si(idxs(s).bst,s,:),1,'omitnan'));
            
            fits_sr_dstr_sn(s,:) = squeeze(mean(fits.slope_sr_dstr(idxs(s).bst,s,:),1,'omitnan'));
            fits_ru_dstr_sn(s,:) = squeeze(mean(fits.slope_ru_dstr(idxs(s).bst,s,:),1,'omitnan'));
            fits_si_dstr_sn(s,:) = squeeze(mean(fits.slope_si_dstr(idxs(s).bst,s,:),1,'omitnan'));
            
            fits_sr_dstr_con_sn(s,:) = squeeze(mean(fits.slope_sr_dstr_con(idxs(s).bst,s,:),1,'omitnan'));
            fits_ru_dstr_con_sn(s,:) = squeeze(mean(fits.slope_ru_dstr_con(idxs(s).bst,s,:),1,'omitnan'));
            fits_si_dstr_con_sn(s,:) = squeeze(mean(fits.slope_si_dstr_con(idxs(s).bst,s,:),1,'omitnan'));
            fits_sr_dstr_inc_sn(s,:) = squeeze(mean(fits.slope_sr_dstr_inc(idxs(s).bst,s,:),1,'omitnan'));
            fits_ru_dstr_inc_sn(s,:) = squeeze(mean(fits.slope_ru_dstr_inc(idxs(s).bst,s,:),1,'omitnan'));
            fits_si_dstr_inc_sn(s,:) = squeeze(mean(fits.slope_si_dstr_inc(idxs(s).bst,s,:),1,'omitnan'));
        end
        fits_sr = cat(3,fits_sr, fits_sr_sn);
        fits_ru = cat(3,fits_ru, fits_ru_sn);
        fits_si = cat(3,fits_si, fits_si_sn);
        
        fits_sr_dstr = cat(3,fits_sr_dstr, fits_sr_dstr_sn);
        fits_ru_dstr = cat(3,fits_ru_dstr, fits_ru_dstr_sn);
        fits_si_dstr = cat(3,fits_si_dstr, fits_si_dstr_sn);
        
        fits_sr_dstr_con = cat(3,fits_sr_dstr_con, fits_sr_dstr_con_sn);
        fits_ru_dstr_con = cat(3,fits_ru_dstr_con, fits_ru_dstr_con_sn);
        fits_si_dstr_con = cat(3,fits_si_dstr_con, fits_si_dstr_con_sn);
        fits_sr_dstr_inc = cat(3,fits_sr_dstr_inc, fits_sr_dstr_inc_sn);
        fits_ru_dstr_inc = cat(3,fits_ru_dstr_inc, fits_ru_dstr_inc_sn);
        fits_si_dstr_inc = cat(3,fits_si_dstr_inc, fits_si_dstr_inc_sn);


    
    % added following m3 (03/20/2025)
    % bootstrap results to match congruency within each cued task and PC (for plotting each cued dimension separately)
    idxs_singlecued = cong_control_ST_singlecued_function(behav,bootstrap_n);
    % bootstrap results to match congruency levels within distractor and PC (independently on each non-cued dimension - already quite balanced by default) (for plotting each distractor independently)
    [idxs, idxs_singledstr] = cong_control_ST_singledstr_function(behav,bootstrap_n);
    

    fits_sr_bytask_sn = nan(4,230,4); % preallocate sessions, times, and cued task
    fits_ru_bytask_sn = nan(4,230,4);
    fits_si_bytask_sn = nan(4,230,4);
    fits_sr_singledstr_sn = nan(4,230,4);
    fits_ru_singledstr_sn = nan(4,230,4);
    fits_si_singledstr_sn = nan(4,230,4);
        for s=1:4 % loop over session (PC)             
            for dist=1:4
                fits_sr_singledstr_sn(s,:,dist) = squeeze(mean(fits.slope_sr_singledstr(idxs_singledstr(s,dist).bst,s,:,dist),1,'omitnan'));
                fits_ru_singledstr_sn(s,:,dist) = squeeze(mean(fits.slope_ru_singledstr(idxs_singledstr(s,dist).bst,s,:,dist),1,'omitnan'));
                fits_si_singledstr_sn(s,:,dist) = squeeze(mean(fits.slope_si_singledstr(idxs_singledstr(s,dist).bst,s,:,dist),1,'omitnan'));
            end
            for task=1:4
                fits_sr_bytask_sn(s,:,task) = squeeze(mean(fits.slope_sr(idxs_singlecued(s,task).bst,s,:),1,'omitnan'));
                fits_ru_bytask_sn(s,:,task) = squeeze(mean(fits.slope_ru(idxs_singlecued(s,task).bst,s,:),1,'omitnan'));
                fits_si_bytask_sn(s,:,task) = squeeze(mean(fits.slope_si(idxs_singlecued(s,task).bst,s,:),1,'omitnan'));
            end
        end
    fits_sr_singledstr = cat(4,fits_sr_singledstr, fits_sr_singledstr_sn);
    fits_ru_singledstr = cat(4,fits_ru_singledstr, fits_ru_singledstr_sn);
    fits_si_singledstr = cat(4,fits_si_singledstr, fits_si_singledstr_sn);
    fits_sr_bytask = cat(4,fits_sr_bytask, fits_sr_bytask_sn);
    fits_ru_bytask = cat(4,fits_ru_bytask, fits_ru_bytask_sn);
    fits_si_bytask = cat(4,fits_si_bytask, fits_si_bytask_sn);

    
    end



fprintf('\tLoaded sbj %g\n', sn) 
end






% plots for singletrials

% permute 
fits_sr = permute(fits_sr,[3,1,2]);
fits_ru = permute(fits_ru,[3,1,2]);
fits_si = permute(fits_si,[3,1,2]);

fits_sr_dstr = permute(fits_sr_dstr,[3,1,2]);
fits_ru_dstr = permute(fits_ru_dstr,[3,1,2]);
fits_si_dstr = permute(fits_si_dstr,[3,1,2]);

fits_sr_dstr_con = permute(fits_sr_dstr_con,[3,1,2]);
fits_ru_dstr_con = permute(fits_ru_dstr_con,[3,1,2]);
fits_si_dstr_con = permute(fits_si_dstr_con,[3,1,2]);
fits_sr_dstr_inc = permute(fits_sr_dstr_inc,[3,1,2]);
fits_ru_dstr_inc = permute(fits_ru_dstr_inc,[3,1,2]);
fits_si_dstr_inc = permute(fits_si_dstr_inc,[3,1,2]);

fits_sr_singledstr = permute(fits_sr_singledstr,[4,1,2,3]);
fits_ru_singledstr = permute(fits_ru_singledstr,[4,1,2,3]);
fits_si_singledstr = permute(fits_si_singledstr,[4,1,2,3]);
fits_sr_bytask = permute(fits_sr_bytask,[4,1,2,3]);
fits_ru_bytask = permute(fits_ru_bytask,[4,1,2,3]);
fits_si_bytask = permute(fits_si_bytask,[4,1,2,3]);


titles={'Neu1', 'MI', 'MC', 'Neu2', 'MI-MC', 'MI-NeuAVG', 'MC-NeuAVG'};
contrasts = {2,3;2,[1,4];3,[1,4]};
sesn = 4;

% reported 2025
% store and export to R, for plotting timeseries - assumes stratify='PC'; 
% single-trial fits are balanced with cong_control function (balance cng within PC), and then averaged within PC 
% fair to average timeseries across PCs
avg.fits_sr = fits_sr;
avg.fits_ru = fits_ru;
avg.fits_si = fits_si;
avg.fits_sr_dstr = fits_sr_dstr;
avg.fits_ru_dstr = fits_ru_dstr; 
avg.fits_si_dstr = fits_si_dstr;

% save Fig4_m1_RSA_betas.mat avg

% single-trial fits stratified by dimension (cued, or distractor) 
% since intent of this plot is to show RSA fit by dimension across PCs, bootstrapping is unnecessary and not reccomended.
avg2.fits_sr_bytask = fits_sr_bytask;
avg2.fits_ru_bytask = fits_ru_bytask;
avg2.fits_si_bytask = fits_si_bytask;
avg2.fits_sr_singledstr = fits_sr_singledstr;
avg2.fits_ru_singledstr = fits_ru_singledstr; 
avg2.fits_si_singledstr = fits_si_singledstr;

% save Fig4_m1_RSA_betas_bydim_bst0.mat avg2          % run with bst=0


        %% cued dim

        colors = num2cell([spring(3)],2);

        if CLR
            ylims = [-1 1]; %[-2 3.6]
        else
            ylims = [-0.04 0.06];        
        end

        
        figure
        sgtitle('RSA models fits for cued dimensions', fontsize=20)
        for ses = 1:sesn
        subplot(2,4,ses)
        h = plot(times,squeeze([mean(fits_sr(:,ses,:),1); mean(fits_ru(:,ses,:),1); mean(fits_si(:,ses,:),1)])', 'linewidth', 2);
        set(h, {'color'}, colors);
        hold on
        plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
        ylim(ylims)
        xlim([times(1) times(end)])
        xlabel('Time [ms]', 'FontSize', 18)
        ylabel('RSA beta', 'FontSize', 18)
        set(gca,'FontSize', 18)
        legend('SR', 'rule', 'side', 'location', 'south'), 
        title(titles{ses}, 'FontSize', 18);
        hold off
        end
            subplot(2,4,5)
            plot(times,squeeze([(mean(fits_sr(:,contrasts{1,1},:),1)-mean(fits_sr(:,contrasts{1,2},:),1)), (mean(fits_ru(:,contrasts{1,1},:),1)-mean(fits_ru(:,contrasts{1,2},:),1)), (mean(fits_si(:,contrasts{1,1},:),1)-mean(fits_si(:,contrasts{1,2},:),1))  ]), 'linewidth', 2) 
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]')
            ylabel('RSA beta - MI-MC')
            legend('SR', 'rule', 'side', 'location', 'north') 
            title(titles(5), 'FontSize', 18);
                subplot(2,4,6)
            plot(times,squeeze([(mean(fits_sr(:,contrasts{2,1},:),1)-mean(fits_sr(:,contrasts{2,2},:),1:2)), (mean(fits_ru(:,contrasts{2,1},:),1)-mean(fits_ru(:,contrasts{2,2},:),1:2)), (mean(fits_si(:,contrasts{2,1},:),1)-mean(fits_si(:,contrasts{2,2},:),1:2))  ]), 'linewidth', 2) 
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MI-NeuAVG', 'FontSize', 18)
            legend('SR', 'rule', 'side', 'location', 'north', 'FontSize', 18) 
            title(titles(6), 'FontSize', 18);
                subplot(2,4,7)
            plot(times,squeeze([(mean(fits_sr(:,contrasts{3,1},:),1)-mean(fits_sr(:,contrasts{3,2},:),1:2)), (mean(fits_ru(:,contrasts{3,1},:),1)-mean(fits_ru(:,contrasts{3,2},:),1:2)), (mean(fits_si(:,contrasts{3,1},:),1)-mean(fits_si(:,contrasts{3,2},:),1:2))  ]), 'linewidth', 2)
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MC-NeuAVG', 'FontSize', 18)
            legend('SR', 'rule', 'side', 'location', 'north', 'FontSize', 18)
            title(titles(7), 'FontSize', 18);
        
    
        RSA_permutest(fits_ru,times)    



        %% distractor
        
        colors = num2cell([autumn(3)],2);
        
        if CLR
            ylims = [-1 1]; 
        else
            ylims = [-0.04 0.06]; 
        end
        
        figure
        sgtitle('RSA models fits for non-cued dimensions', fontsize=20)
        for ses = 1:sesn
        subplot(2,4,ses)
        h = plot(times,squeeze([mean(fits_sr_dstr(:,ses,:),1); mean(fits_ru_dstr(:,ses,:),1); mean(fits_si_dstr(:,ses,:),1)])', 'linewidth', 2);
        set(h, {'color'}, colors);
        hold on
        plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
        ylim(ylims)
        xlim([times(1) times(end)])
        xlabel('Time [ms]', 'FontSize', 18)
        ylabel('RSA beta', 'FontSize', 18)
        set(gca,'FontSize', 18)
        legend('SR', 'rule', 'side', 'location', 'northeast'), 
        title(titles{ses}, 'FontSize', 18);
        hold off
        end
            subplot(2,4,5)
            plot(times,squeeze([(mean(fits_sr_dstr(:,contrasts{1,1},:),1)-mean(fits_sr_dstr(:,contrasts{1,2},:),1)), (mean(fits_ru_dstr(:,contrasts{1,1},:),1)-mean(fits_ru_dstr(:,contrasts{1,2},:),1)), (mean(fits_si_dstr(:,contrasts{1,1},:),1)-mean(fits_si_dstr(:,contrasts{1,2},:),1))  ]), 'linewidth', 2) 
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MI-MC', 'FontSize', 18)
            set(gca,'FontSize', 18)
            legend('SR', 'rule', 'side', 'location', 'northeast') 
            title(titles(5), 'FontSize', 18);
                subplot(2,4,6)
            plot(times,squeeze([(mean(fits_sr_dstr(:,contrasts{2,1},:),1)-mean(fits_sr_dstr(:,contrasts{2,2},:),1:2)), (mean(fits_ru_dstr(:,contrasts{2,1},:),1)-mean(fits_ru_dstr(:,contrasts{2,2},:),1:2)), (mean(fits_si_dstr(:,contrasts{2,1},:),1)-mean(fits_si_dstr(:,contrasts{2,2},:),1:2))  ]), 'linewidth', 2) 
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MI-NeuAVG', 'FontSize', 18)
            set(gca,'FontSize', 18)
            legend('SR', 'rule', 'side', 'location', 'northeast') 
            title(titles(6), 'FontSize', 18);
                subplot(2,4,7)
            plot(times,squeeze([(mean(fits_sr_dstr(:,contrasts{3,1},:),1)-mean(fits_sr_dstr(:,contrasts{3,2},:),1:2)), (mean(fits_ru_dstr(:,contrasts{3,1},:),1)-mean(fits_ru_dstr(:,contrasts{3,2},:),1:2)), (mean(fits_si_dstr(:,contrasts{3,1},:),1)-mean(fits_si_dstr(:,contrasts{3,2},:),1:2))  ]), 'linewidth', 2)
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MC-NeuAVG', 'FontSize', 18)
            set(gca,'FontSize', 18)
            legend('SR', 'rule', 'side', 'location', 'northeast')
            title(titles(7), 'FontSize', 18);
        
    
        RSA_permutest(fits_sr_dstr,times,stratify)    


        %% distractor congruent

        if CLR
            ylims = [-1 1]; 
        else
            ylims = [-0.04 0.06]; 
        end

        figure
        sgtitle('RSA models fits for non-cued dimensions - when congruent', fontsize=20)
        for ses = 1:sesn
        subplot(2,4,ses)
        h = plot(times,squeeze([mean(fits_sr_dstr_con(:,ses,:),1); mean(fits_ru_dstr_con(:,ses,:),1); mean(fits_si_dstr_con(:,ses,:),1)])', 'linewidth', 2);
        set(h, {'color'}, colors);
        hold on
        plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
        ylim(ylims)
        xlim([times(1) times(end)])
        xlabel('Time [ms]', 'FontSize', 18)
        ylabel('RSA beta', 'FontSize', 18)
        set(gca,'FontSize', 18)
        legend('SR', 'rule', 'side', 'location', 'southwest', 'FontSize', 18), 
        title(titles{ses}, 'FontSize', 18);
        hold off
        end
        
            subplot(2,4,5)
            plot(times,squeeze([(mean(fits_sr_dstr_con(:,contrasts{1,1},:),1)-mean(fits_sr_dstr_con(:,contrasts{1,2},:),1)), (mean(fits_ru_dstr_con(:,contrasts{1,1},:),1)-mean(fits_ru_dstr_con(:,contrasts{1,2},:),1)), (mean(fits_si_dstr_con(:,contrasts{1,1},:),1)-mean(fits_si_dstr_con(:,contrasts{1,2},:),1))  ]), 'linewidth', 2) 
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MI-MC', 'FontSize', 18)
            set(gca,'FontSize', 18)
            legend('SR', 'rule', 'side', 'location', 'southwest', 'FontSize', 18) 
            title(titles(5), 'FontSize', 18);
                subplot(2,4,6)
            plot(times,squeeze([(mean(fits_sr_dstr_con(:,contrasts{2,1},:),1)-mean(fits_sr_dstr_con(:,contrasts{2,2},:),1:2)), (mean(fits_ru_dstr_con(:,contrasts{2,1},:),1)-mean(fits_ru_dstr_con(:,contrasts{2,2},:),1:2)), (mean(fits_si_dstr_con(:,contrasts{2,1},:),1)-mean(fits_si_dstr_con(:,contrasts{2,2},:),1:2))  ]), 'linewidth', 2) 
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MI-NeuAVG', 'FontSize', 18)
            set(gca,'FontSize', 18)
            legend('SR', 'rule', 'side', 'location', 'southwest', 'FontSize', 18) 
            title(titles(6), 'FontSize', 18);
                subplot(2,4,7)
            plot(times,squeeze([(mean(fits_sr_dstr_con(:,contrasts{3,1},:),1)-mean(fits_sr_dstr_con(:,contrasts{3,2},:),1:2)), (mean(fits_ru_dstr_con(:,contrasts{3,1},:),1)-mean(fits_ru_dstr_con(:,contrasts{3,2},:),1:2)), (mean(fits_si_dstr_con(:,contrasts{3,1},:),1)-mean(fits_si_dstr_con(:,contrasts{3,2},:),1:2))  ]), 'linewidth', 2)
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MC-NeuAVG', 'FontSize', 18)
            set(gca,'FontSize', 18)
            legend('SR', 'rule', 'side', 'location', 'southwest', 'FontSize', 18)
            title(titles(7), 'FontSize', 18);
        
    
        RSA_permutest(fits_ru_dstr_con,times)    

        %% distractor incongruent

        if CLR
            ylims = [-1 1]; 
        else
            ylims = [-0.04 0.06]; 
        end
        
        figure
        sgtitle('RSA models fits for non-cued dimensions - when incongruent', fontsize=20)
        for ses = 1:sesn
        subplot(2,4,ses)
        h = plot(times,squeeze([mean(fits_sr_dstr_inc(:,ses,:),1); mean(fits_ru_dstr_inc(:,ses,:),1); mean(fits_si_dstr_inc(:,ses,:),1)])', 'linewidth', 2);
        set(h, {'color'}, colors);
        hold on
        plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
        ylim(ylims)
        xlim([times(1) times(end)])
        xlabel('Time [ms]', 'FontSize', 18)
        ylabel('RSA beta', 'FontSize', 18)
        set(gca,'FontSize', 18)
        legend('SR', 'rule', 'side', 'location', 'southwest', 'FontSize', 18), 
        title(titles{ses}, 'FontSize', 18);
        hold off
        end
            subplot(2,4,5)
            plot(times,squeeze([(mean(fits_sr_dstr_inc(:,contrasts{1,1},:),1)-mean(fits_sr_dstr_inc(:,contrasts{1,2},:),1)), (mean(fits_ru_dstr_inc(:,contrasts{1,1},:),1)-mean(fits_ru_dstr_inc(:,contrasts{1,2},:),1)), (mean(fits_si_dstr_inc(:,contrasts{1,1},:),1)-mean(fits_si_dstr_inc(:,contrasts{1,2},:),1))  ]), 'linewidth', 2) 
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MI-MC', 'FontSize', 18)
            set(gca,'FontSize', 18)
            legend('SR', 'rule', 'side', 'location', 'southwest', 'FontSize', 18) 
            title(titles(5), 'FontSize', 18);
                subplot(2,4,6)
            plot(times,squeeze([(mean(fits_sr_dstr_inc(:,contrasts{2,1},:),1)-mean(fits_sr_dstr_inc(:,contrasts{2,2},:),1:2)), (mean(fits_ru_dstr_inc(:,contrasts{2,1},:),1)-mean(fits_ru_dstr_inc(:,contrasts{2,2},:),1:2)), (mean(fits_si_dstr_inc(:,contrasts{2,1},:),1)-mean(fits_si_dstr_inc(:,contrasts{2,2},:),1:2))  ]), 'linewidth', 2) 
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MI-NeuAVG', 'FontSize', 18)
            set(gca,'FontSize', 18)
            legend('SR', 'rule', 'side', 'location', 'southwest', 'FontSize', 18) 
            title(titles(6), 'FontSize', 18);
                subplot(2,4,7)
            plot(times,squeeze([(mean(fits_sr_dstr_inc(:,contrasts{3,1},:),1)-mean(fits_sr_dstr_inc(:,contrasts{3,2},:),1:2)), (mean(fits_ru_dstr_inc(:,contrasts{3,1},:),1)-mean(fits_ru_dstr_inc(:,contrasts{3,2},:),1:2)), (mean(fits_si_dstr_inc(:,contrasts{3,1},:),1)-mean(fits_si_dstr_inc(:,contrasts{3,2},:),1:2))  ]), 'linewidth', 2)
            hold on
            plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta - MC-NeuAVG', 'FontSize', 18)
            set(gca,'FontSize', 18)
            legend('SR', 'rule', 'side', 'location', 'southwest', 'FontSize', 18)
            title(titles(7), 'FontSize', 18);
        
    
        RSA_permutest(fits_ru_dstr_inc,times)    
    











%% single cued dimensions' (based on plot_RSA_betas_m3.m)
        % Note: If bst==1, balanced cng by cued taks and PC. But, some missing conditions caues some cng levels by ses by task have 0 trial count 
        % Therefore, since intent of this plot is to show RSA fit by dimension across PCs, bootstrapping is unnecessary and not reccomended.

        colors = num2cell([autumn(4)],2);
        ylims = [-1 1.3]; 
        
        figure
        sgtitle('RSA models fits for cued dimensions by cued task', fontsize=20)
        subplot(1,3,1)
        h = plot(times,squeeze([mean(fits_sr_bytask(:,:,:,1),1:2); mean(fits_sr_bytask(:,:,:,2),1:2); mean(fits_sr_bytask(:,:,:,3),1:2);  mean(fits_sr_bytask(:,:,:,4),1:2)])', 'linewidth', 2);
        set(h, {'color'}, colors);
        hold on
        plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
        ylim(ylims)
        xlim([times(1) times(end)])
        xlabel('Time [ms]', 'FontSize', 18)
        ylabel('RSA beta', 'FontSize', 18)
        set(gca,'FontSize', 18)
        legend('C SR', 'S SR', 'E SR', 'P SR', 'location', 'southwest', 'FontSize', 18), 
        hold off

        subplot(1,3,2)
        h = plot(times,squeeze([mean(fits_ru_bytask(:,:,:,1),1:2); mean(fits_ru_bytask(:,:,:,2),1:2); mean(fits_ru_bytask(:,:,:,3),1:2);  mean(fits_ru_bytask(:,:,:,4),1:2)])', 'linewidth', 2);
        set(h, {'color'}, colors);
        hold on
        plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
        ylim(ylims)
        xlim([times(1) times(end)])
        xlabel('Time [ms]', 'FontSize', 18)
        ylabel('RSA beta', 'FontSize', 18)
        set(gca,'FontSize', 18)
        legend('C Rule', 'S Rule', 'E Rule', 'P Rule', 'location', 'southwest', 'FontSize', 18), 
        hold off

        subplot(1,3,3)
        h = plot(times,squeeze([mean(fits_si_bytask(:,:,:,1),1:2); mean(fits_si_bytask(:,:,:,2),1:2); mean(fits_si_bytask(:,:,:,3),1:2);  mean(fits_si_bytask(:,:,:,4),1:2)])', 'linewidth', 2);
        set(h, {'color'}, colors);
        hold on
        plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
        ylim(ylims)
        xlim([times(1) times(end)])
        xlabel('Time [ms]', 'FontSize', 18)
        ylabel('RSA beta', 'FontSize', 18)
        set(gca,'FontSize', 18)
        legend('C Side', 'S Side', 'E Side', 'P Side', 'location', 'southwest', 'FontSize', 18), 
        hold off
            
            






        %% single distractors
        % Note: If bst==1, cng is balanced for each distractor (but already ~50% by design). Bootstrapping is unnecessary and not reccomended
        
        colors = num2cell([autumn(4)],2);
        
        if CLR
            ylims = [-0.7 0.3]; 
        else
            ylims = [-0.03 0.03]; 
        end
        
        figure
        sgtitle('RSA models fits for non-cued dimensions', fontsize=20)
        subplot(1,3,1)
        h = plot(times,squeeze([mean(fits_sr_singledstr(:,:,:,1),1:2); mean(fits_sr_singledstr(:,:,:,2),1:2); mean(fits_sr_singledstr(:,:,:,3),1:2);  mean(fits_sr_singledstr(:,:,:,4),1:2)])', 'linewidth', 2);
        set(h, {'color'}, colors);
        hold on
        plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
        ylim(ylims)
        xlim([times(1) times(end)])
        xlabel('Time [ms]', 'FontSize', 18)
        ylabel('RSA beta', 'FontSize', 18)
        set(gca,'FontSize', 18)
        legend('C SR', 'S SR', 'E SR', 'P SR', 'location', 'southwest', 'FontSize', 18), 
        hold off

        subplot(1,3,2)
        h = plot(times,squeeze([mean(fits_ru_singledstr(:,:,:,1),1:2); mean(fits_ru_singledstr(:,:,:,2),1:2); mean(fits_ru_singledstr(:,:,:,3),1:2);  mean(fits_ru_singledstr(:,:,:,4),1:2)])', 'linewidth', 2);
        set(h, {'color'}, colors);
        hold on
        plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
        ylim(ylims)
        xlim([times(1) times(end)])
        xlabel('Time [ms]', 'FontSize', 18)
        ylabel('RSA beta', 'FontSize', 18)
        set(gca,'FontSize', 18)
        legend('C Rule', 'S Rule', 'E Rule', 'P Rule', 'location', 'southwest', 'FontSize', 18), 
        hold off

        subplot(1,3,3)
        h = plot(times,squeeze([mean(fits_si_singledstr(:,:,:,1),1:2); mean(fits_si_singledstr(:,:,:,2),1:2); mean(fits_si_singledstr(:,:,:,3),1:2);  mean(fits_si_singledstr(:,:,:,4),1:2)])', 'linewidth', 2);
        set(h, {'color'}, colors);
        hold on
        plot(times,zeros(length(times),1),'color','black','LineStyle',':','linewidth',2)
        ylim(ylims)
        xlim([times(1) times(end)])
        xlabel('Time [ms]', 'FontSize', 18)
        ylabel('RSA beta', 'FontSize', 18)
        set(gca,'FontSize', 18)
        legend('C Side', 'S Side', 'E Side', 'P Side', 'location', 'southwest', 'FontSize', 18), 
        hold off
            
