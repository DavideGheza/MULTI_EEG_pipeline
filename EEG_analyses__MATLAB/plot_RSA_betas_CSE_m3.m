%% Prepare and plot congruency-sequence effects in model-3 RSA fits
% Loads single-trial model-3 conjunction estimates, applies matched
% bootstrap stratification, and prepares the CSE arrays used in Figure 5.
% The selected input determines whether ordinary or QR-orthogonalized fits
% are analyzed.
% Author: Davide Gheza

%% 1. Load settings

clear

RSAadaptations = [];
RSAadaptations_on_bydist = [];
RSAadaptations_param_bytarget = [];
RSAadaptations_param = [];
% for RT bin analyses
RSAadaptations_b = [];
RSAadaptations_on_bydist_b = [];
RSAadapt_htmp_b = [];
array_filenames = [    {'CS_Ix'}    {'CS_Cx'}    {'CS_II'}    {'CS_CI'}    {'CS_IC'}    {'CS_CC'}];
for arf = 1:numel(array_filenames)
    RSAadapts_b.(array_filenames{arf}) = [];
    RSAadapts_on_bydist_b.(array_filenames{arf}) = [];
    RSAadapts_htmp_b.(array_filenames{arf}) = [];
end

adapt_param_bytarget_byRT = [];
store = nan(4,4,3,26);


IX_ondiag = [];
IX_offdiag = [];
CX_ondiag = [];
CX_offdiag = [];
RSAadapt_htmp = [];

fits_sr_singledstr = [];


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
load(fullfile(sets.rsa_data_dir, ['RSA_pc/multiclass/singletrial/fits/fits_singletrials_m3_sbj_' num2str(sbjIDX(sn)) '.mat'])) 
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
    
    
    % CSE: adaptation by single distractor. Use Neu PC conditions
    % preallocate 
    IX_idx = nan(4,length(congruency));
    CX_idx = nan(4,length(congruency));
    RSAadapt_ondiag = nan(1,length(times));
    RSAadapt_offdiag = nan(1,length(times));
    RSAadapt_on_bydist = nan(4,length(times));
    CSE = nan(2,4,length(times),4);
    RSAadapt_bytarget = nan(4,length(times),4); % CSE parametric (spaghetti shift) on cued dim, split by each cued dim
    prev_SC_idx = []; % allocate before each use below
    RSAadapt_param = nan(4,length(times)); % CSE parametric on distractor
    IX_ondiag_sn = nan(1,length(times));
    IX_offdiag_sn = nan(1,length(times));
    CX_ondiag_sn = nan(1,length(times));
    CX_offdiag_sn = nan(1,length(times));

    % 1) (non-cued) distractor adaptation
    meas = fits.slope_sr_singledstr(idxs,PC_idx,:,:); % select PC
    % Indexes to compute adaptation on RSA fits. Post-switch trials are removed.
    % Loop over the four dimensions, considered with regard to previous cng. 
    % With filter by congruency==1|0, trials are dropped when dim is cued (unnecessary cause no singledistractor fit for those).
    % Remove post-error trials

    for prev_ncdim=1:4
        IX_idx(prev_ncdim,:) = prev_congruency(:,prev_ncdim)==0 & congruency(:,prev_ncdim)~=-1 & task_count>1 & post_correct==1;% note: when bst=1, exclusion (task_count>1 & post_correct==1) are redundant
        CX_idx(prev_ncdim,:) = prev_congruency(:,prev_ncdim)==1 & congruency(:,prev_ncdim)~=-1 & task_count>1 & post_correct==1;
    end
    CSE_idx = cat(3,IX_idx,CX_idx);

    
    
    % compute adaptation on RSA fit
    % average over selected trials
    for prev_ncdim=1:4
        for CSEtype=1:2
            trial_idx = logical(CSE_idx(prev_ncdim,:,CSEtype));
            CSE(CSEtype,prev_ncdim,:,:) = mean(meas(trial_idx,:,:,:), [1,2], 'omitnan'); % 4th dim: test adaptation on each singledistractor. 2nd dim: include trials from both neu PCs.
        end
    end

    % Prev Inc vs prev Con
    RSAadapt = squeeze(CSE(1,:,:,:) - CSE(2,:,:,:));
    RSAadapt = permute(RSAadapt,[1,3,2]); % previous dimension interference, current singledistractor, times
    CSE = permute(CSE,[1,2,4,3]);
    I = eye(4);
    for t = 1:size(RSAadapt,3)
        slice = RSAadapt(:,:,t);
        RSAadapt_ondiag(t) = mean(diag(slice));
        RSAadapt_offdiag(t) = mean(slice(~I));
        RSAadapt_on_bydist(:,t) = diag(slice);
        % save out also each post-congruency levels
        slice = squeeze(CSE(1,:,:,t));
        IX_ondiag_sn(t) = mean(diag(slice));
        IX_offdiag_sn(t) = mean(slice(~I));
        slice = squeeze(CSE(2,:,:,t));
        CX_ondiag_sn(t) = mean(diag(slice));
        CX_offdiag_sn(t) = mean(slice(~I));
    end

    RSAadaptations_sn = [RSAadapt_ondiag;RSAadapt_offdiag]; % on vs off diagonal adaptation over time
    RSAadaptations = cat(3, RSAadaptations, RSAadaptations_sn);
    RSAadaptations_on_bydist = cat(3, RSAadaptations_on_bydist, RSAadapt_on_bydist);
    IX_ondiag = cat(3, IX_ondiag, IX_ondiag_sn); % 1st dim is 1, for catenating on/off later
    IX_offdiag = cat(3, IX_offdiag, IX_offdiag_sn);
    CX_ondiag = cat(3, CX_ondiag, CX_ondiag_sn);
    CX_offdiag = cat(3, CX_offdiag, CX_offdiag_sn);
    RSAadapt_htmp = cat(4, RSAadapt_htmp, RSAadapt);


    % 1.1) Dynamcics: RT bins. (non-cued) distractor adaptation.
    % Indexes to compute adaptation on RSA fits. 
    % In addition to above, filter by RT bin.
    cuts = 6;
    rt = behav.rt(idxs);
    quantileAssignment = zeros(4, 2, size(rt, 1)); 
    quants = nan(4,2,cuts);
    for prev_ncdim = 1:4
        for prev_cng = 1:2
            idx = prev_congruency(:,prev_ncdim)==prev_cng-1 & congruency(:,prev_ncdim)~=-1 & task_count>1 & post_correct==1;
            if nnz(idx) > 0
                quants(prev_ncdim,prev_cng,:) = quantile(rt(idx),cuts);
                quantileAssignment(prev_ncdim,prev_cng,:) = discretize(rt, [-Inf, squeeze(quants(prev_ncdim,prev_cng,:))', Inf]);
            end
        end
    end
    
    % preallocate 
    IX_idx = nan(length(congruency),1);
    CX_idx = nan(length(congruency),1);
    RSAadapt_ondiag_b = nan(cuts+1,length(times));
    RSAadapt_offdiag_b = nan(cuts+1,length(times));
    RSAadapt_on_bydist_b = nan(cuts+1,4,length(times));
    % average over trials
    CSE = nan([2,4,cuts+1,size(meas,[3,4])]); % prev cng, prev nc dim, RT bins, times, current singledistractor
    store_counts = nan([2,4,cuts+1,26]);
    for bin = 1:cuts+1
        for prev_ncdim=1:4
            quantile_cond_IX = squeeze(quantileAssignment(prev_ncdim,1,:)==bin);
            quantile_cond_CX = squeeze(quantileAssignment(prev_ncdim,2,:)==bin);
            IX_idx = prev_congruency(:,prev_ncdim)==0 & congruency(:,prev_ncdim)~=-1 & task_count>1 & post_correct==1 & quantile_cond_IX;% note: when bst=1, exclusion (task_count>1 & post_correct==1) are redundant
            CX_idx = prev_congruency(:,prev_ncdim)==1 & congruency(:,prev_ncdim)~=-1 & task_count>1 & post_correct==1 & quantile_cond_CX;
            CSE(1,prev_ncdim,bin,:,:) = mean(meas(IX_idx,:,:,:), [1,2], 'omitnan'); % IX.  4th dim: test adaptation on each singledistractor. 2nd dim: include trials from both neu PCs.
            CSE(2,prev_ncdim,bin,:,:) = mean(meas(CX_idx,:,:,:), [1,2], 'omitnan'); % CX.  4th dim: test adaptation on each singledistractor. 2nd dim: include trials from both neu PCs.
            store_counts(1,prev_ncdim,bin,sn) = sum(IX_idx);
            store_counts(2,prev_ncdim,bin,sn) = sum(CX_idx);
        end
    end
    % Prev Inc vs prev Con, by RT bin
    RSAadapt = squeeze(CSE(1,:,:,:,:) - CSE(2,:,:,:,:));
    RSAadapt = permute(RSAadapt,[1,4,3,2]); % previous dimension interference, current singledistractor, times, RT bin
    CSE = permute(CSE,[1,2,5,4,3]);
    I = eye(4);
    for bin = 1:cuts+1
        for t = 1:size(RSAadapt,3)
            slice = RSAadapt(:,:,t,bin);
            RSAadapt_ondiag_b(bin,t) = mean(diag(slice));
            RSAadapt_offdiag_b(bin,t) = mean(slice(~I));
            RSAadapt_on_bydist_b(bin,:,t) = diag(slice);
        end
    end

    RSAadaptations_sn = cat(3, RSAadapt_ondiag_b,RSAadapt_offdiag_b); % on vs off diagonal adaptation over time, by RT bin
    RSAadaptations_b = cat(4, RSAadaptations_b, RSAadaptations_sn);
    RSAadaptations_on_bydist_b = cat(4, RSAadaptations_on_bydist_b, RSAadapt_on_bydist_b);
    RSAadapt_htmp_b = cat(5, RSAadapt_htmp_b, RSAadapt);
    

    % 1.2) Dynamcics: RT bins. (non-cued) distractor adaptation. 2x2 prev-cur (quantiles are computed WRT current congruency, which has stronger effect on RTs)
    cuts = 6;
    rt = behav.rt(idxs);
    quantileAssignment = zeros(4, 2, size(rt, 1)); 
    quants = nan(4,2,cuts);
    for cur_ncdim = 1:4 % cur_ncdim: compared to above contrast (pre cng), now compute quantiles on current cng (and stratify pre cng using same quantiles)
        for cur_cng = 1:2
            idx = congruency(:,cur_ncdim)==cur_cng-1 & task_count>1 & post_correct==1; % & congruency(:,cur_ncdim)~=-1
            if nnz(idx) > 0
                quants(cur_ncdim,cur_cng,:) = quantile(rt(idx),cuts);
                quantileAssignment(cur_ncdim,cur_cng,:) = discretize(rt, [-Inf, squeeze(quants(cur_ncdim,cur_cng,:))', Inf]);
            end
        end
    end
    
    % preallocate 
    II_idx = nan(length(congruency),1);
    CI_idx = nan(length(congruency),1);
    IC_idx = nan(length(congruency),1);
    CC_idx = nan(length(congruency),1);
    cs = [];
    % average over trials
    CSE = nan([2,4,cuts+1,size(meas,[3,4])]); % prev cng, prev nc dim, RT bins, times, current singledistractor
    for bin = 1:cuts+1
        for cur_ncdim=1:4
            quantile_cond_xI = squeeze(quantileAssignment(cur_ncdim,1,:)==bin); % again: quantile is based on current cng (xI and xC)
            quantile_cond_xC = squeeze(quantileAssignment(cur_ncdim,2,:)==bin);
            II_idx = prev_congruency(:,cur_ncdim)==0 & congruency(:,cur_ncdim)==0 & task_count>1 & post_correct==1 & quantile_cond_xI; % stratify by congruency in pre and cur trial WRT same dimension. I.e., congruency sequence within dimension
            CI_idx = prev_congruency(:,cur_ncdim)==1 & congruency(:,cur_ncdim)==0 & task_count>1 & post_correct==1 & quantile_cond_xI;
            IC_idx = prev_congruency(:,cur_ncdim)==0 & congruency(:,cur_ncdim)==1 & task_count>1 & post_correct==1 & quantile_cond_xC;
            CC_idx = prev_congruency(:,cur_ncdim)==1 & congruency(:,cur_ncdim)==1 & task_count>1 & post_correct==1 & quantile_cond_xC;
            % save single post-x conditions (used to compute CSE in previous 1.1 - but different quantile approach here)
            cs.CS_Ix(cur_ncdim,bin,:,:) = mean(meas(II_idx|IC_idx,:,:,:), [1,2], 'omitnan'); % IX.  
            cs.CS_Cx(cur_ncdim,bin,:,:) = mean(meas(CI_idx|CC_idx,:,:,:), [1,2], 'omitnan'); % CX.  
            % save 2x2 post-con/inc x curr-con/inc
            cs.CS_II(cur_ncdim,bin,:,:) = mean(meas(II_idx,:,:,:), [1,2], 'omitnan');
            cs.CS_CI(cur_ncdim,bin,:,:) = mean(meas(CI_idx,:,:,:), [1,2], 'omitnan');
            cs.CS_IC(cur_ncdim,bin,:,:) = mean(meas(IC_idx,:,:,:), [1,2], 'omitnan'); 
            cs.CS_CC(cur_ncdim,bin,:,:) = mean(meas(CC_idx,:,:,:), [1,2], 'omitnan'); 
        end
    end
    % Prev Inc vs prev Con, by RT bin
    fn = fieldnames(cs);
    for c = 1:numel(fn)
        RSAadapt_ondiag_b = nan(cuts+1,length(times));
        RSAadapt_offdiag_b = nan(cuts+1,length(times));
        RSAadapt_on_bydist_b = nan(cuts+1,4,length(times));
        cond = cs.(fn{c});
        RSAadapt = permute(cond,[1,4,3,2]); % previous dimension interference, current singledistractor, times, RT bin
        I = eye(4);
        for bin = 1:cuts+1
            for t = 1:size(RSAadapt,3)
                slice = RSAadapt(:,:,t,bin);
                RSAadapt_ondiag_b(bin,t) = mean(diag(slice));
                RSAadapt_offdiag_b(bin,t) = mean(slice(~I));
                RSAadapt_on_bydist_b(bin,:,t) = diag(slice);
            end
        end
    
        RSAadapts_sn = cat(3, RSAadapt_ondiag_b,RSAadapt_offdiag_b); % on vs off diagonal adaptation over time, by RT bin
        RSAadapts_b.(fn{c}) = cat(4, RSAadapts_b.(fn{c}), RSAadapts_sn);
        RSAadapts_on_bydist_b.(fn{c}) = cat(4, RSAadapts_on_bydist_b.(fn{c}), RSAadapt_on_bydist_b);
        RSAadapts_htmp_b.(fn{c}) = cat(5, RSAadapts_htmp_b.(fn{c}), RSAadapt);
    end




    
 
    
    % 2) (non-cued) parametric congruency sequence effect on distractors (average SR)
    meas = fits.slope_sr_dstr(idxs,PC_idx,:); % select PC 
    % Indexes to compute adaptation on RSA fits.
    % Equivalent of CSE spaghetti. As for CSE heatmap, only previous congruency matter - no need to cross prev. with cur. sum congruency.
    % Post-switch trials are removed.
    % Loop over the four prev_sum_cong. 
    % Remove post-error trials
    prev_SC_idx = nan(4,length(task));
    for prev_SC=1:4
        prev_SC_idx(prev_SC,:) = prev_sum_con==prev_SC-1 & task_count>1 & post_correct==1;
    end
        
    % average over selected trials
    for prev_SC=1:4
        trial_idx = logical(prev_SC_idx(prev_SC,:));
        RSAadapt_param(prev_SC,:) = mean(meas(trial_idx,:,:), [1,2], 'omitnan'); % 2nd dim: include trials from selected PCs.
    end
    
    RSAadaptations_param = cat(3, RSAadaptations_param, RSAadapt_param);


    
    % 3) (cued) target adaptation - parametric congruency sequence effect
    fits_slope_sr = fits.slope_sr(idxs,:,:);
    meas = fits_slope_sr(:,PC_idx,:); % select PC 
    % Indexes to compute adaptation on RSA fits. 
    % Equivalent of CSE spaghetti. As for CSE heatmap, only previous congruency matter - no need to cross prev. with cur. sum congruency.
    % Post-switch trials are removed.
    % Loop over the four prev_sum_cong. 
    % Remove post-error trials
    % With filter by task==cued_dim, trials are dropped when dim is not cued.
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
            RSAadapt_bytarget(prev_SC,:,cued_dim) = mean(meas(trial_idx,:,:), [1,2], 'omitnan'); % test adaptation on each cued dim. 2nd dim: include trials from both neu PCs.
        end
    end

    RSAadaptations_param_bytarget = cat(4, RSAadaptations_param_bytarget, RSAadapt_bytarget);



    % 3.1) Dynamcics: RT bins. (cued) target adaptation - parametric congruency sequence effect
    
    % Indexes to compute adaptation on RSA fits. 
    % In addition to above, filter by RT bin.
    cuts = 6;
    rt = behav.rt(idxs);
    quantileAssignment = zeros(4, 4, size(rt, 1)); 
    quants = nan(4,4,cuts);
    for cued_dim = 1:4
        for prev_SC = 1:4
            idx = prev_sum_con==prev_SC-1 & task==cued_dim & task_count>1 & post_correct==1;
            if nnz(idx) > 0
                quants(cued_dim,prev_SC,:) = quantile(rt(idx),cuts);
                quantileAssignment(cued_dim,prev_SC,:) = discretize(rt, [-Inf, squeeze(quants(cued_dim,prev_SC,:))', Inf]);
            end
        end
    end
    
    % average over selected trials
    RSAadapt_bytarget_byRT = nan([4,size(meas,3),4,cuts+1]);
    for bin = 1:cuts+1
        for cued_dim=1:4
            for prev_SC=1:4
                quantile_cond = squeeze(quantileAssignment(cued_dim,prev_SC,:)==bin);
                trial_idx = (prev_sum_con==prev_SC-1 & task==cued_dim & task_count>1 & post_correct==1 & quantile_cond);
                RSAadapt_bytarget_byRT(prev_SC,:,cued_dim,bin) = mean(meas(trial_idx,:,:), [1,2], 'omitnan'); 
                store(prev_SC,cued_dim,bin,sn) = sum(trial_idx(quantile_cond));
            end
        end
    end

    adapt_param_bytarget_byRT = cat(5, adapt_param_bytarget_byRT, RSAadapt_bytarget_byRT);






    
end


fprintf('\tLoaded sbj %g\n', sn) 
end








% reported 2025                                                 (08/08/25)
% store and export to R, for plotting timeseries - assumes stratify='CSE'; 
% single-trial fits are balanced with cong_control_ST_precur_function (balance cng within PC, controlling for cng of current AND previous trial), stratified by previous congruency, and then averaged within PC 
% fair to average timeseries across PCs
avg.RSAadaptations = RSAadaptations; %  on/off diagonal: 2   230    26
avg.RSAadaptations_on_bydist = RSAadaptations_on_bydist; %  on diagonal by distractor: 4   230    26
avg.RSAadapt_htmp = RSAadapt_htmp; % previous dimension interference, current singledistractor, times, sbj: 4     4   230    26
avg.RSAadaptations_param = RSAadaptations_param; % parametric adaptation: non-cued avg by previous sum con: 4   230    26
avg.on_off_byprev_cng = cat(4, cat(1, IX_ondiag, IX_offdiag), cat(1, CX_ondiag, CX_offdiag));
avg.on_off_byprev_cng = permute(avg.on_off_byprev_cng, [4,1,2,3]);                              % prev cng, on/off, time, sbj: 2 2 230 26

% save Fig5_SFig2_CSE_m3_RSA_betas.mat avg                % run with bst=1









        %% CSE on individual distractors: fixed-window heatmap
        figure

        timebnds = [100 400];
        timewin = dsearchn(times',timebnds(1)) : dsearchn(times',timebnds(2));

        d = mean(RSAadapt_htmp(:,:,timewin,:),[3,4],'omitnan');
        heatmap({'C','S','E','P'},{'C','S','E','P'},d')
        set(gca,'fontsize',14)
        colormap turbo
        ylabel('Current distractor RSA beta')
        xlabel('Previous distractor (congruent-incongruent)')

        customBlueColormap = [linspace(0, 1, 256)', linspace(0, 1, 256)', ones(256, 1)];
        colormap(customBlueColormap)

        % clim([-0.1 0.05])

        


        %% CSE on single distractors ON vs OFF diagonal

 
        fontsize = 18;
        % for poster
        fontsize = 22;


        % h 2.40
        
        % CSE on single distractors ON vs OFF diagonal
        % Compute the mean and SEM
        mean_vals = squeeze(mean(RSAadaptations, 3)); 
        sem_vals = squeeze(std(RSAadaptations, 0, 3)) / sqrt(size(RSAadaptations, 3));
        
        % Create the figure and plot
        figure
        % sgtitle('CSE on distractors ON vs OFF diagonal', 'FontSize', 20)
        
        % Colors for the lines and shaded areas
        colors = [0 0 1; 0.5 0 1]; % Example colors: blue for ON, violet for OFF
        
        ylims = [-0.08 0.031]; 

        hold on
        % Plot mean and SEM as shaded areas
        for i = 1:size(mean_vals, 1)
            % Plot the shaded area for SEM
            fill([times, fliplr(times)], ...
                [mean_vals(i,:) - sem_vals(i,:), fliplr(mean_vals(i,:) + sem_vals(i,:))], ...
                colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            % Plot the mean line
            plot(times, mean_vals(i,:), 'LineWidth', 2, 'Color', colors(i, :));
        end
        plot(times, zeros(length(times), 1), 'Color', 'black', 'LineStyle', ':', 'LineWidth', 2)        % Plot the zero line
        ylim(ylims)
        xlim([times(1) times(end)])
        xlabel('Time [ms]', 'FontSize', fontsize)
        ylabel('\Delta RSA beta', 'FontSize', fontsize)
        set(gca, 'FontSize', fontsize)
        legend('', 'ON diagonal', '', 'OFF diagonal', 'location', 'southeast')        % Add legend
        hold off

        CSE_permutest(RSAadaptations,times)  

        % CSE on single distractors ON diagonal, by distractor (diagonal elements)
        % Compute the mean and SEM
        mean_vals = squeeze(mean(RSAadaptations_on_bydist, 3)); 
        sem_vals = squeeze(std(RSAadaptations_on_bydist, 0, 3)) / sqrt(size(RSAadaptations_on_bydist, 3));
        
        % Create the figure and plot
        figure
        % sgtitle('CSE on diagonal, by single distractors', 'FontSize', 20)
        
        % Colors for the lines and shaded areas
        % colors = viridis(256);
        % quantile_indices = round(linspace(1, size(colors, 1), 5)); % 5 points to get 4 intervals
        % colors = colors(quantile_indices(2:end), :);
        darker_scale = 0.9;
        colors = parula(4);

        hold on
        % Plot mean and SEM as shaded areas
        for i = 1:size(mean_vals, 1)
            % Plot the shaded area for SEM
            fill([times, fliplr(times)], ...
                [mean_vals(i,:) - sem_vals(i,:), fliplr(mean_vals(i,:) + sem_vals(i,:))], ...
                colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            % Plot the mean line
            plot(times, mean_vals(i,:), 'LineWidth', 2, 'Color', colors(i, :) * darker_scale);
        end
        plot(times, zeros(length(times), 1), 'Color', 'black', 'LineStyle', ':', 'LineWidth', 2)        % Plot the zero line
        % ylim(ylims)
        xlim([times(1) times(end)])
        xlabel('Time [ms]', 'FontSize', fontsize)
        ylabel('\Delta RSA beta', 'FontSize', fontsize)
        set(gca, 'FontSize', fontsize)
        legend('', 'C', '', 'S', '', 'E', '', 'P', 'location', 'southeast')        % Add legend
        hold off


        % run permutest on each distractor's adaptation
        labels = {{'C'}, {'S'}, {'E'}, {'P'}};
        for didx = 1:4
            CSE_permutest(RSAadaptations_on_bydist(didx,:,:),times,labels{didx})
        end 

        %% Target parametric adaptation


        % Parametric congruency sequence effect on (averaged) cued dimensions (spaghetti shift)
        % Compute the mean and SEM
        RSAadaptations_target = squeeze(mean(RSAadaptations_param_bytarget,3)); % average of cued dims.
        mean_vals = squeeze(mean(RSAadaptations_target, 3)); 
        sem_vals = squeeze(std(RSAadaptations_target, 0, 3)) / sqrt(size(RSAadaptations_target, 3));
        
        % Create the figure and plot
        figure
        sgtitle('Parametric CSE on cued dimensions (target adaptation)', 'FontSize', 20)
        
        % Colors for the lines and shaded areas
        colors = winter(4);
        
        hold on
        % Plot mean and SEM as shaded areas
        for i = 1:size(mean_vals, 1)
            % Plot the shaded area for SEM
            fill([times, fliplr(times)], ...
                [mean_vals(i,:) - sem_vals(i,:), fliplr(mean_vals(i,:) + sem_vals(i,:))], ...
                colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            % Plot the mean line
            plot(times, mean_vals(i,:), 'LineWidth', 2, 'Color', colors(i, :));
        end
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


        %% Distractor parametric adaptation


        % Parametric congruency sequence effect on (averaged) cued dimensions (spaghetti shift)
        % Compute the mean and SEM
        mean_vals = squeeze(mean(RSAadaptations_param, 3)); 
        sem_vals = squeeze(std(RSAadaptations_param, 0, 3)) / sqrt(size(RSAadaptations_param, 3));
        
        % Create the figure and plot
        figure
        sgtitle('Parametric CSE on distractors', 'FontSize', 20)
        
        % Colors for the lines and shaded areas
        colors = summer(4);
        
        hold on
        % Plot mean and SEM as shaded areas
        for i = 1:size(mean_vals, 1)
            % Plot the shaded area for SEM
            fill([times, fliplr(times)], ...
                [mean_vals(i,:) - sem_vals(i,:), fliplr(mean_vals(i,:) + sem_vals(i,:))], ...
                colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            % Plot the mean line
            plot(times, mean_vals(i,:), 'LineWidth', 2, 'Color', colors(i, :));
        end
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









        %% Target parametric adaptation - by cued dimension (explored 03/31/25)


        % Parametric congruency sequence effect on (averaged) cued dimensions (spaghetti shift)
        % Compute the mean and SEM
        RSAadaptations_target = RSAadaptations_param_bytarget; % DO NOT average of cued dims.
        mean_vals = squeeze(mean(RSAadaptations_target, 4)); 
        sem_vals = squeeze(std(RSAadaptations_target, 0, 4)) / sqrt(size(RSAadaptations_target, 4));
        
        % Create the figure and plot
        figure
        sgtitle('Parametric CSE on cued dimensions (target adaptation)', 'FontSize', 20)
        
        % Colors for the lines and shaded areas
        colors = winter(4);
        
        % Plot mean and SEM as shaded areas
        for sbp = 1:4
            subplot (2,2,sbp)
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



        %% Target parametric adaptation - by RT bin (explored 03/31/25)



        % Parametric congruency sequence effect on (averaged) cued dimensions (spaghetti shift)
        % Compute the mean and SEM
        RSAadaptations_target = squeeze(mean(adapt_param_bytarget_byRT,3,'omitnan')); % fourth dimension is RT bin. Average of cued dims (some NANs in 3 sbjs/cued task due to decimation).
        mean_vals = squeeze(mean(RSAadaptations_target, 4)); 
        sem_vals = squeeze(std(RSAadaptations_target, 0, 4)) / sqrt(size(RSAadaptations_target, 4));
        
        % Create the figure and plot
        figure
        sgtitle('Parametric CSE on cued dimensions (target adaptation)', 'FontSize', 20)
        
        % Colors for the lines and shaded areas
        colors = winter(4);
        
        % Plot mean and SEM as shaded areas
        bins = size(mean_vals, 3);
        for sbp = 1:bins
            subplot (bins,1,sbp)
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
            ylim([-0.2 1])
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('RSA beta', 'FontSize', 18)
            set(gca, 'FontSize', 18)
            lgd = legend('', '0', '', '1', '', '2', '', '3', 'location', 'northwest');        % Add legend
            lgd.Title.String = '# previous congruent';
            lgd.Title.FontSize = 12;
            hold off
        end





         %% CSE on single distractors ON vs OFF diagonal - by RT bin (explored 04/01-03/25)

        fontsize = 18;
        % for poster
        fontsize = 22;

        % CSE on single distractors ON vs OFF diagonal
        % Compute the mean and SEM
        mean_vals = squeeze(mean(RSAadaptations_b, 4)); 
        sem_vals = squeeze(std(RSAadaptations_b, 0, 4)) / sqrt(size(RSAadaptations_b, 4));
        
        % Create the figure and plot
        figure        
        % Colors for the lines and shaded areas
        colors = [0 0 1; 0.5 0 1]; % Example colors: blue for ON, violet for OFF
        
        ylims = [-0.15 0.05]; 

        hold on
        % Plot mean and SEM as shaded areas
        bins = size(mean_vals, 1);
        for sbp = 1:bins
            subplot (bins,1,sbp)
            hold on
            for i = 1:size(mean_vals, 3)
                % Plot the shaded area for SEM
                fill([times, fliplr(times)], ...
                    [squeeze(mean_vals(sbp,:,i)) - squeeze(sem_vals(sbp,:,i)), fliplr(squeeze(mean_vals(sbp,:,i)) + squeeze(sem_vals(sbp,:,i)))], ...
                    colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                % Plot the mean line
                plot(times, squeeze(mean_vals(sbp,:,i)), 'LineWidth', 2, 'Color', colors(i, :));
            end
            plot(times, zeros(length(times), 1), 'Color', 'black', 'LineStyle', ':', 'LineWidth', 2)        % Plot the zero line
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', fontsize)
            ylabel('\Delta RSA beta', 'FontSize', fontsize)
            set(gca, 'FontSize', fontsize)
            legend('', 'ON diagonal', '', 'OFF diagonal', 'location', 'southeast')        % Add legend
            hold off
        end



        %%%%%%%%%%

        % CSE on single distractors ON diagonal, by distractor (diagonal elements)
        % Compute the mean and SEM
        mean_vals = squeeze(mean(RSAadaptations_on_bydist_b, 4)); 
        sem_vals = squeeze(std(RSAadaptations_on_bydist_b, 0, 4)) / sqrt(size(RSAadaptations_on_bydist_b, 4));
        
        % Create the figure and plot
        figure
        % sgtitle('CSE on diagonal, by single distractors', 'FontSize', 20)
        
        % Colors for the lines and shaded areas
        % colors = viridis(256);
        % quantile_indices = round(linspace(1, size(colors, 1), 5)); % 5 points to get 4 intervals
        % colors = colors(quantile_indices(2:end), :);
        darker_scale = 0.9;
        colors = parula(4);

        bins = size(mean_vals, 1);
        for sbp = 1:bins
            subplot (bins,1,sbp)
            hold on
            % Plot mean and SEM as shaded areas
            for i = 1:size(mean_vals, 2)
                shaded_area_x = [times, fliplr(times)];
                shaded_area_y = [squeeze(mean_vals(sbp, i, :))' - squeeze(sem_vals(sbp, i, :))', ...
                                 fliplr(squeeze(mean_vals(sbp, i, :))' + squeeze(sem_vals(sbp, i, :))')];
                fill(shaded_area_x, shaded_area_y, ...
                     colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                % Plot the mean line
                plot(times, squeeze(mean_vals(sbp, i, :)), 'LineWidth', 2, 'Color', colors(i, :) * darker_scale);
            end
            plot(times, zeros(length(times), 1), 'Color', 'black', 'LineStyle', ':', 'LineWidth', 2)        % Plot the zero line
            % ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', fontsize)
            ylabel('\Delta RSA beta', 'FontSize', fontsize)
            set(gca, 'FontSize', fontsize)
            legend('', 'C', '', 'S', '', 'E', '', 'P', 'location', 'southeast')        % Add legend
            hold off
        end


%% CSE on single distractors ON vs OFF diagonal - by RT bin - individual conditions (04/11/25)


        fontsize = 18;
        % for poster
        fontsize = 22;

titles = [    {'Ix'}    {'Cx'}    {'II'}    {'CI'}    {'IC'}    {'CC'}];
        
for arf = 1:numel(array_filenames)


        % CSE on single distractors ON vs OFF diagonal
        % Compute the mean and SEM
        mean_vals = squeeze(mean(RSAadapts_b.(array_filenames{arf}), 4)); 
        sem_vals = squeeze(std(RSAadapts_b.(array_filenames{arf}), 0, 4)) / sqrt(size(RSAadapts_b.(array_filenames{arf}), 4));
        
        % Create the figure and plot
        figure        
        % Colors for the lines and shaded areas
        colors = [0 0 1; 0.5 0 1]; % Example colors: blue for ON, violet for OFF
        if any(arf==[1,2])
            ylims = [-0.07 0.15]; 
        else
            ylims = [-0.65 0.8]; 
        end

        hold on
        % Plot mean and SEM as shaded areas
        bins = size(mean_vals, 1);
        for sbp = 1:bins
            subplot (bins,1,sbp)
            hold on
            for i = 1:size(mean_vals, 3)
                % Plot the shaded area for SEM
                fill([times, fliplr(times)], ...
                    [squeeze(mean_vals(sbp,:,i)) - squeeze(sem_vals(sbp,:,i)), fliplr(squeeze(mean_vals(sbp,:,i)) + squeeze(sem_vals(sbp,:,i)))], ...
                    colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                % Plot the mean line
                plot(times, squeeze(mean_vals(sbp,:,i)), 'LineWidth', 2, 'Color', colors(i, :));
            end
            plot(times, zeros(length(times), 1), 'Color', 'black', 'LineStyle', ':', 'LineWidth', 2)        % Plot the zero line
            ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', fontsize)
            ylabel('\Delta RSA beta', 'FontSize', fontsize)
            set(gca, 'FontSize', fontsize)
            % legend('', 'ON diagonal', '', 'OFF diagonal', 'location', 'best')        % Add legend
            sgtitle(titles{arf}, 'FontSize', 20)
            hold off
        end


end







        %%%%%%%%%%

        % CSE on single distractors ON diagonal, by distractor (diagonal elements)
        
        darker_scale = 0.9;
        colors = parula(4);
        bins = size(mean_vals, 1);

for arf = 1:numel(array_filenames)
        
        figure
        mean_vals = squeeze(mean(RSAadapts_on_bydist_b.(array_filenames{arf}), 4)); 
        sem_vals = squeeze(std(RSAadapts_on_bydist_b.(array_filenames{arf}), 0, 4)) / sqrt(size(RSAadapts_on_bydist_b.(array_filenames{arf}), 4));

        for sbp = 1:bins
            subplot (bins,1,sbp)
            hold on
            % Plot mean and SEM as shaded areas
            for i = 1:size(mean_vals, 2)
                shaded_area_x = [times, fliplr(times)];
                shaded_area_y = [squeeze(mean_vals(sbp, i, :))' - squeeze(sem_vals(sbp, i, :))', ...
                                 fliplr(squeeze(mean_vals(sbp, i, :))' + squeeze(sem_vals(sbp, i, :))')];
                fill(shaded_area_x, shaded_area_y, ...
                     colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                % Plot the mean line
                plot(times, squeeze(mean_vals(sbp, i, :)), 'LineWidth', 2, 'Color', colors(i, :) * darker_scale);
            end
            plot(times, zeros(length(times), 1), 'Color', 'black', 'LineStyle', ':', 'LineWidth', 2)        % Plot the zero line
            % ylim(ylims)
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', fontsize)
            ylabel('\Delta RSA beta', 'FontSize', fontsize)
            set(gca, 'FontSize', fontsize)
            legend('', 'C', '', 'S', '', 'E', '', 'P', 'location', 'southeast')        % Add legend
            sgtitle(titles{arf}, 'FontSize', 20)
            hold off
        end
end
