%% Plot congruency-sequence effects in model-1 RSA fits
% Generates the panel F in Figure 4.
% Stratifies single-trial model-1 RSA estimates by current and previous
% congruency using matched bootstrap samples.
% Author: Davide Gheza


%% 1. Load settings

clear


RSAadaptations_param_bytarget = [];



% Stratify results for the congruency-sequence analysis.
stratify='CSE';
% chose PC conditions
PC_idx = ([1,4]); % select Neu1 and Neu2 PC 
% PC_idx = (1:4); % all PCs
% PC_idx = ([2,3]); % MI and MC
% bootstrap RSA fits to balance congruency within each PC? (equal n of cng0,1,2,3 oversampling trials of less numerous conditions)
bst = 1;
bootstrap_n = 100;
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
load(fullfile(sets.rsa_data_dir, ['RSA_pc/multiclass/singletrial/fits/fits_singletrials_m1_sbj_' num2str(sbjIDX(sn)) '.mat'])) 
keep(sn).fits = fits; % all single trial data
% load behavioral data
behav = load([sets.sessions_data_dir 'V3_behav_seqinfo_trimmed_sbj_' num2str(sbjIDX(sn))]); % behavioral file with addidional info on subsequent and previous trials (10/01/2024)
behav = behav.behav_trim;

if strcmp(stratify,'PC') | strcmp(stratify,'CSE')
    if ~bst
    idxs = 1:length(behav.task); % index all the trials without bootstrap
    else    
    % bootstrap results to match congruency within PC (decimate and draw with replacement) # -> implemented for non-Neu PCs, but good also to ensure cng balance in Neu conditions
    idxs_bypc = cong_control_ST_precur_function(behav,bootstrap_n); % controlling for cng of current AND previous trial 
    idxs = vertcat(idxs_bypc(:).bst)'; % catenate indexes from all PCs
    end
    
    task_count = behav.taskcount(idxs);
    congruency = behav.congruency(idxs,:);
    prev_congruency = behav.prev_congruency(idxs,:);
    post_correct = behav.post_correct(idxs);
    prev_sum_con = behav.prev_sum_con(idxs);
    task = behav.task(idxs);
    
    
    % CSE
    % preallocate 
    RSAadapt_bytarget = nan(3,4,length(times),4); % CSE parametric (spaghetti shift) on cued dim, split by each cued dim. Firt dimension is axis


    % (cued) target adaptation - parametric congruency sequence effect - for M1 include Cue, Choice, and Feature coding
    fits_slope_sr = fits.slope_sr(idxs,:,:);
    fits_slope_ru = fits.slope_ru(idxs,:,:);
    fits_slope_si = fits.slope_si(idxs,:,:);

    meas = [];
    meas(1,:,:,:) = fits_slope_sr(:,PC_idx,:); % select PC 
    meas(2,:,:,:) = fits_slope_ru(:,PC_idx,:);
    meas(3,:,:,:) = fits_slope_si(:,PC_idx,:);
    % Indexes to compute adaptation on RSA fits. Loop over the four prev_sum_cong. 
    % Remove post-error trials and post-switch trials. With filter by task==cued_dim, trials are dropped when dim is not cued.
    prev_SC_idx = nan(4,length(task),4);
    for cued_dim = 1:4
        for prev_SC=1:4
            prev_SC_idx(prev_SC,:,cued_dim) = prev_sum_con==prev_SC-1 & task==cued_dim & task_count>1 & post_correct==1;
        end
    end
    % average over selected trials
    for cued_dim=1:4
        for prev_SC=1:4
            trial_idx = logical(prev_SC_idx(prev_SC,:,cued_dim));
            for m = 1:3
                RSAadapt_bytarget(m,prev_SC,:,cued_dim) = mean(meas(m,trial_idx,:,:), [2,3], 'omitnan'); % test adaptation on each cued dim. 3rd dim: include trials from both neu PCs.
            end
        end
    end

    RSAadaptations_param_bytarget = cat(5, RSAadaptations_param_bytarget, RSAadapt_bytarget);
    




    
end


fprintf('\tLoaded sbj %g\n', sn) 
end












        %% Target parametric adaptation - by axis (Cue, Choice, Feature)


        % Compute the mean and SEM
        RSAadaptations_target = squeeze(mean(RSAadaptations_param_bytarget,4)); % average of cued dims.
        mean_vals = squeeze(mean(RSAadaptations_target, 4)); 
        sem_vals = squeeze(std(RSAadaptations_target, 0, 4)) / sqrt(size(RSAadaptations_target, 4));
        
        % Create the figure and plot
        figure
        sgtitle('Parametric CSE on cued dimensions (target adaptation)', 'FontSize', 20)
        
        % Colors for the lines and shaded areas
        colors = winter(4);
        
        mean_vals = permute(mean_vals,[2,3,1]);
        sem_vals = permute(sem_vals,[2,3,1]);
        % Plot mean and SEM as shaded areas
        for sbp = 1:3
            subplot (3,1,sbp)
            hold on
            for i = 1:size(mean_vals, 1)
                % Plot the shaded area for SEM
                fill([times, fliplr(times)], ...
                    [squeeze(mean_vals(i,:,sbp)) - squeeze(sem_vals(i,:,sbp)), fliplr(squeeze(mean_vals(i,:,sbp)) + squeeze(sem_vals(i,:,sbp)))], ...
                    colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                % Plot the mean line
                plot(times, squeeze(mean_vals(i,:,sbp)), 'LineWidth', 2, 'Color', colors(i, :));
            end
            plot(times, zeros(length(times), 1), 'Color', 'black', 'LineStyle', ':', 'LineWidth', 2)        % Plot the zero line
            
            
            plot(times, zeros(length(times), 1), 'Color', 'black', 'LineStyle', ':', 'LineWidth', 2)        % Plot the zero line
            % ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta', 'FontSize', 18)
            set(gca, 'FontSize', 18)
            lgd = legend('', '0', '', '1', '', '2', '', '3', 'location', 'northwest');        % Add legend
            lgd.Title.String = '# previous congruent';
            lgd.Title.FontSize = 12;
            hold off
        end



   