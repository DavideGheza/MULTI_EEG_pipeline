%% Prepare full-model encoding projections for CSE analyses
% Loads single-trial projections onto the stimulus-response (SR), spatial
% (SI), dimension-general conflict, dimension-specific conflict, previous
% conflict, and cue axes. Stratifies projections by current and previous
% congruency for SI Figure S3 and the response-locked SI Figure S4 analysis.
% Manuscript terminology: dimension-general and dimension-specific axes are
% described as dimension-shared and dimension-private, respectively.
% Author: Davide Gheza


%% 1. Load settings

clear

basic_saveout = 0; % save out basic stratification by task/distractor (all PCs, no bst), vs CSE stratification (Neu PCs, with bootstrapping) --> for exporting to R. For local plotting here, set to '0'

% Model selection - choose which axis to plot
model='SR'; % 'SR' 'SI' 'DGC' 'DSC' 'PrevC' 'Cue'
if strcmp('SR',model)
    ax_name = 'encoded_dim_spec_SR';
    file_modifier = 'dimSpec';
elseif strcmp('SI',model)
    ax_name = 'encoded_dim_gen_SI';
    file_modifier = 'dimGen';
elseif strcmp('DGC',model)
    ax_name = 'encoded_dim_gen_conf';
    file_modifier = 'DGC';
elseif strcmp('DSC',model)
    ax_name = 'encoded_dim_spec_conf';
    file_modifier = 'DSC';
elseif strcmp('PrevC',model)
    ax_name = 'encoded_prev_conf';
    file_modifier = 'PrevC';
elseif strcmp('Cue',model)
    ax_name = 'encoded_cue';
    file_modifier = 'Cue';
end


% chose PC conditions
if basic_saveout
    PC_idx = (1:4); % all PCs
else
    PC_idx = ([1,4]); % select Neu1 and Neu2 PC 
end
% PC_idx = ([2,3]); % MI and MC

% Bootstrap settings
if basic_saveout
    bst = 0; 
else
    bst = 1; % apply bootstrapping to balance congruency within PC
end
bootstrap_n = 100;

% Time or tf domain
tf = true;
onset = true;
if onset
    suffix = '';
else
    suffix = 'response_';   % response-locked
end

% Load settings
if ispc
    load('C:\Users\gdavide\OneDrive\Documenti\MATLAB\dg_workspace22\dimensions\analyses\scripts_and_settings\preprocessing_settings.mat');
else
    load('/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/scripts_and_settings/preprocessing_settings_mac.mat');
end


sets.rsa_data_dir=[sets.main_saving_dir 'rsa_data/'];

sbjIDX = [22,24:48];% subjects list

if tf
    times = -800:10:1498;
else
    times = -1000:2:1999;
end


% Preallocate
ENCadaptations = [];
ENCadaptations_on_bydist = [];
ENCadaptations_param_bytarget = [];
ENCadaptations_param = [];
% for RT bin analyses
ENCadaptations_b = [];
ENCadaptations_on_bydist_b = [];
ENCadapt_htmp_b = [];
array_filenames = [    {'CS_Ix'}    {'CS_Cx'}    {'CS_II'}    {'CS_CI'}    {'CS_IC'}    {'CS_CC'}];
for arf = 1:numel(array_filenames)
    ENCadapts_b.(array_filenames{arf}) = [];
    ENCadapts_on_bydist_b.(array_filenames{arf}) = [];
    ENCadapts_htmp_b.(array_filenames{arf}) = [];
end

adapt_param_bytarget_byRT = [];
store = nan(4,4,3,26);


IX_ondiag = [];
IX_offdiag = [];
CX_ondiag = [];
CX_offdiag = [];
ENCadapt_htmp = [];

fits_sr_singledstr = [];

cued_proj_all = nan(4,230,26); % store trial average per sbj to inpsect axes
noncued_proj_all = nan(4,230,4,26);

cued_bytask_all = nan(26,230,4); % for basic saveout
dstr_bytask_all = nan(26,230,4,4);

%% 2. Load single-trial projections and apply bootstrapping

for sn = 1:length(sbjIDX)
    
    % load behavioral data
    behav = load([sets.sessions_data_dir 'V3_behav_seqinfo_trimmed_' suffix 'sbj_' num2str(sbjIDX(sn))]); % behavioral file with addidional info on subsequent and previous trials (10/01/2024)
    behav = behav.behav_trim;
    % load single trial projections (~ RSA fits)
    data_st = load(fullfile(sets.rsa_data_dir, ['RSA_pc/multiclass/singletrial/encoding_fullmodel_proj/encoding_fullmodel_notempgen_' suffix 'sbj_' num2str(sbjIDX(sn)) '.mat']));    
     
    % Resolve single-trial source from result struct
    if isfield(data_st, 'result') && isfield(data_st.result, ax_name)
        result_ax = data_st.result.(ax_name);
    else
        error('Missing single-trial result field for axis %s in subject %d.', ax_name, sbjIDX(sn));
    end

    % result_ax is a 1x4 struct with field .s [trial x time x dim]
    projections = cat(4, result_ax(1).s, result_ax(2).s, result_ax(3).s, result_ax(4).s);% cat the 4 PCs 
    projections = permute(projections, [1, 4, 2, 3]); % [trial, pc, time, dist]


    
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
    
    % Parse cued from non-cued projections
    noncued_proj = projections(idxs,:,:,:);
    cued_proj = nan(size(noncued_proj,(1:3)));
    for t = 1:numel(idxs)
        d = task(t);
        cued_proj(t,:,:) = projections(idxs(t),:,:,d);  
        noncued_proj(t,:,:,d) = NaN;
    end
    
    % CSE: adaptation by single distractor. Use Neu PC conditions
    % preallocate 
    IX_idx = nan(4,length(congruency));
    CX_idx = nan(4,length(congruency));
    ENCadapt_ondiag = nan(1,length(times));
    ENCadapt_offdiag = nan(1,length(times));
    ENCadapt_on_bydist = nan(4,length(times));
    CSE = nan(2,4,length(times),4);
    ENCadapt_bytarget = nan(4,length(times),4); % CSE parametric (spaghetti shift) on cued dim, split by each cued dim
    prev_SC_idx = []; % allocate before each use below
    ENCadapt_param = nan(4,length(times)); % CSE parametric on distractor
    IX_ondiag_sn = nan(1,length(times));
    IX_offdiag_sn = nan(1,length(times));
    CX_ondiag_sn = nan(1,length(times));
    CX_offdiag_sn = nan(1,length(times));

    % 1) (non-cued) distractor adaptation
    meas = noncued_proj(:,PC_idx,:,:); % select PC
    % Indexes to compute adaptation on ENC projections. Post-switch trials are removed.
    % Loop over the four dimensions, considered with regard to previous cng. 
    % With filter by congruency==1|0, trials are dropped when dim is cued (unnecessary cause no singledistractor fit for those).
    % Remove post-error trials

    for prev_ncdim=1:4
        IX_idx(prev_ncdim,:) = prev_congruency(:,prev_ncdim)==0 & congruency(:,prev_ncdim)~=-1 & task_count>1 & post_correct==1;% note: when bst=1, exclusion (task_count>1 & post_correct==1) are redundant
        CX_idx(prev_ncdim,:) = prev_congruency(:,prev_ncdim)==1 & congruency(:,prev_ncdim)~=-1 & task_count>1 & post_correct==1;
    end
    CSE_idx = cat(3,IX_idx,CX_idx);

    
    
    % compute adaptation on ENC projections
    % average over selected trials
    for prev_ncdim=1:4
        for CSEtype=1:2
            trial_idx = logical(CSE_idx(prev_ncdim,:,CSEtype));
            CSE(CSEtype,prev_ncdim,:,:) = mean(meas(trial_idx,:,:,:), [1,2], 'omitnan'); % 4th dim: test adaptation on each singledistractor. 2nd dim: include trials from both neu PCs.
        end
    end

    % Prev Inc vs prev Con
    ENCadapt = squeeze(CSE(1,:,:,:) - CSE(2,:,:,:));
    ENCadapt = permute(ENCadapt,[1,3,2]); % previous dimension interference, current singledistractor, times
    CSE = permute(CSE,[1,2,4,3]);
    I = eye(4);
    for t = 1:size(ENCadapt,3)
        slice = ENCadapt(:,:,t);
        ENCadapt_ondiag(t) = mean(diag(slice));
        ENCadapt_offdiag(t) = mean(slice(~I));
        ENCadapt_on_bydist(:,t) = diag(slice);
        % save out also each post-congruency levels
        slice = squeeze(CSE(1,:,:,t));
        IX_ondiag_sn(t) = mean(diag(slice));
        IX_offdiag_sn(t) = mean(slice(~I));
        slice = squeeze(CSE(2,:,:,t));
        CX_ondiag_sn(t) = mean(diag(slice));
        CX_offdiag_sn(t) = mean(slice(~I));
    end

    ENCadaptations_sn = [ENCadapt_ondiag;ENCadapt_offdiag]; % on vs off diagonal adaptation over time
    ENCadaptations = cat(3, ENCadaptations, ENCadaptations_sn);
    ENCadaptations_on_bydist = cat(3, ENCadaptations_on_bydist, ENCadapt_on_bydist);
    IX_ondiag = cat(3, IX_ondiag, IX_ondiag_sn); % 1st dim is 1, for catenating on/off later
    IX_offdiag = cat(3, IX_offdiag, IX_offdiag_sn);
    CX_ondiag = cat(3, CX_ondiag, CX_ondiag_sn);
    CX_offdiag = cat(3, CX_offdiag, CX_offdiag_sn);
    ENCadapt_htmp = cat(4, ENCadapt_htmp, ENCadapt);


    % 1.1) Dynamcics: RT bins. (non-cued) distractor adaptation.
    % Indexes to compute adaptation on ENC projections.
    % In addition to above, filter by RT bin.
    cuts = 2;
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
    ENCadapt_ondiag_b = nan(cuts+1,length(times));
    ENCadapt_offdiag_b = nan(cuts+1,length(times));
    ENCadapt_on_bydist_b = nan(cuts+1,4,length(times));
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
    ENCadapt = squeeze(CSE(1,:,:,:,:) - CSE(2,:,:,:,:));
    ENCadapt = permute(ENCadapt,[1,4,3,2]); % previous dimension interference, current singledistractor, times, RT bin
    CSE = permute(CSE,[1,2,5,4,3]);
    I = eye(4);
    for bin = 1:cuts+1
        for t = 1:size(ENCadapt,3)
            slice = ENCadapt(:,:,t,bin);
            ENCadapt_ondiag_b(bin,t) = mean(diag(slice));
            ENCadapt_offdiag_b(bin,t) = mean(slice(~I));
            ENCadapt_on_bydist_b(bin,:,t) = diag(slice);
        end
    end

    ENCadaptations_sn = cat(3, ENCadapt_ondiag_b,ENCadapt_offdiag_b); % on vs off diagonal adaptation over time, by RT bin
    ENCadaptations_b = cat(4, ENCadaptations_b, ENCadaptations_sn);
    ENCadaptations_on_bydist_b = cat(4, ENCadaptations_on_bydist_b, ENCadapt_on_bydist_b);
    ENCadapt_htmp_b = cat(5, ENCadapt_htmp_b, ENCadapt);
    

    % 1.2) Dynamcics: RT bins. (non-cued) distractor adaptation. 2x2 prev-cur (quantiles are computed WRT current congruency, which has stronger effect on RTs)
    cuts = 2;
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
        ENCadapt_ondiag_b = nan(cuts+1,length(times));
        ENCadapt_offdiag_b = nan(cuts+1,length(times));
        ENCadapt_on_bydist_b = nan(cuts+1,4,length(times));
        cond = cs.(fn{c});
        ENCadapt = permute(cond,[1,4,3,2]); % previous dimension interference, current singledistractor, times, RT bin
        I = eye(4);
        for bin = 1:cuts+1
            for t = 1:size(ENCadapt,3)
                slice = ENCadapt(:,:,t,bin);
                ENCadapt_ondiag_b(bin,t) = mean(diag(slice));
                ENCadapt_offdiag_b(bin,t) = mean(slice(~I));
                ENCadapt_on_bydist_b(bin,:,t) = diag(slice);
            end
        end
    
        ENCadapts_sn = cat(3, ENCadapt_ondiag_b,ENCadapt_offdiag_b); % on vs off diagonal adaptation over time, by RT bin
        ENCadapts_b.(fn{c}) = cat(4, ENCadapts_b.(fn{c}), ENCadapts_sn);
        ENCadapts_on_bydist_b.(fn{c}) = cat(4, ENCadapts_on_bydist_b.(fn{c}), ENCadapt_on_bydist_b);
        ENCadapts_htmp_b.(fn{c}) = cat(5, ENCadapts_htmp_b.(fn{c}), ENCadapt);
    end




    


    
    % 2) (non-cued) parametric congruency sequence effect on distractors (average SR)
    meas = mean(noncued_proj(:,PC_idx,:,:),4,'omitnan'); % select PC, mean over non-cued dim
    % Indexes to compute adaptation on ENC projections.
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
        ENCadapt_param(prev_SC,:) = mean(meas(trial_idx,:,:), [1,2], 'omitnan'); % 2nd dim: include trials from selected PCs.
    end
    
    ENCadaptations_param = cat(3, ENCadaptations_param, ENCadapt_param);


    
    % 3) (cued) target adaptation - parametric congruency sequence effect
    meas = cued_proj(:,PC_idx,:); % select PC 
    % Indexes to compute adaptation on ENC projections.
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
            ENCadapt_bytarget(prev_SC,:,cued_dim) = mean(meas(trial_idx,:,:), [1,2], 'omitnan'); % test adaptation on each cued dim. 2nd dim: include trials from both neu PCs.
        end
    end

    ENCadaptations_param_bytarget = cat(4, ENCadaptations_param_bytarget, ENCadapt_bytarget);



    % 3.1) Dynamcics: RT bins. (cued) target adaptation - parametric congruency sequence effect
    
    % Indexes to compute adaptation on ENC projections. 
    % In addition to above, filter by RT bin.
    cuts = 2;
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
    ENCadapt_bytarget_byRT = nan([4,size(meas,3),4,cuts+1]);
    for bin = 1:cuts+1
        for cued_dim=1:4
            for prev_SC=1:4
                quantile_cond = squeeze(quantileAssignment(cued_dim,prev_SC,:)==bin);
                trial_idx = (prev_sum_con==prev_SC-1 & task==cued_dim & task_count>1 & post_correct==1 & quantile_cond);
                ENCadapt_bytarget_byRT(prev_SC,:,cued_dim,bin) = mean(meas(trial_idx,:,:), [1,2], 'omitnan'); 
                store(prev_SC,cued_dim,bin,sn) = sum(trial_idx(quantile_cond));
            end
        end
    end

    adapt_param_bytarget_byRT = cat(5, adapt_param_bytarget_byRT, ENCadapt_bytarget_byRT);


    
    cued_proj_all(:,:,sn) = squeeze(mean(cued_proj,1,'omitnan'));
    noncued_proj_all(:,:,:,sn) = squeeze(mean(noncued_proj,1,'omitnan'));

    
    % To save out basic stratification by task and distractor. Parse each dimension. 
    dstr_bytask_sn = nan(4,230,4,4); % preallocate PC, times, non-cued dim, and cued tasks
    cued_bytask_sn = nan(4,230,4); % preallocate PC, times, and cued task
    for task_id=1:4
        dstr_bytask_sn(:,:,:,task_id) = squeeze(mean(noncued_proj(task==task_id,:,:,:),1,'omitnan')); % averge over trials, create cued task dim
        cued_bytask_sn(:,:,task_id) = squeeze(mean(cued_proj(task==task_id,:,:),1,'omitnan'));
    end
    % average over PC
    dstr_bytask_all(sn,:,:,:) = mean(dstr_bytask_sn(:,:,:,:),1,'omitnan');
    cued_bytask_all(sn,:,:) = mean(cued_bytask_sn(:,:,:),1,'omitnan'); 

        
fprintf('\tLoaded sbj %g\n', sn) 
end




% same save-out as plot_RSA_betas_CSE_m3
% reported 2026                                                 (05/17/26)
% store and export to R, for plotting timeseries
if ~basic_saveout
    % single-trial fits are balanced with cong_control_ST_precur_function (balance cng within PC, controlling for cng of current AND previous trial), stratified by previous congruency, and then averaged within PC 
    % fair to average timeseries across PCs
    avg.ENCadaptations = ENCadaptations; %  on/off diagonal: 2   230    26
    avg.ENCadaptations_on_bydist = ENCadaptations_on_bydist; %  on diagonal by distractor: 4   230    26
    avg.ENCadapt_htmp = ENCadapt_htmp; % previous dimension interference, current singledistractor, times, sbj: 4     4   230    26
    avg.ENCadaptations_param = ENCadaptations_param; % parametric adaptation: non-cued avg by previous sum con: 4   230    26
    avg.on_off_byprev_cng = cat(4, cat(1, IX_ondiag, IX_offdiag), cat(1, CX_ondiag, CX_offdiag));
    avg.on_off_byprev_cng = permute(avg.on_off_byprev_cng, [4,1,2,3]);                              % prev cng, on/off, time, sbj: 2 2 230 26
    
    if strcmp(model,'DSC')
        fname = 'notreported_CSE_encoding_DSC.mat';
    elseif strcmp(model,'SR') || strcmp(model,'SI')
        fname = ['SFig3_CSE_encoding_' file_modifier '.mat'];
    end
    save(fname, 'avg')                 % run with bst=1

else
    % for plotting basic stratification (run including all PCs, no bst)
    % single-trial projections are stratified by dimension (cued, or distractor) 
    avg.cued_bytask_all = cued_bytask_all; % 3d: sbjs, time, task (projection)
    avg.dstr_bytask_all = dstr_bytask_all; % 4d: sbjs, time, dist (projection), task
    
    if onset
        fname = ['SFig3_basic_encoding_' file_modifier '_bst0.mat']; % run with bst=0. Average across all PCs
    else
        fname = ['SFig4_resplocked_encoding_' file_modifier '_bst0.mat'];
    end
    save(fname, 'avg')
end

%% plot simple projections (across PC, for cued / each non cued dimension)

fontsize = 18;
colors = parula(4);

cued_proj_squezed = squeeze(mean(cued_proj_all,1, 'omitnan'));       % average over PCs
noncued_proj_squezed = squeeze(mean(noncued_proj_all,1, 'omitnan')); % average over PCs
mean_vals = squeeze(mean(cued_proj_squezed, 2)); 
mean_vals_2 = squeeze(mean(noncued_proj_squezed, 3)); 
sem_vals = squeeze(std(cued_proj_squezed, 0, 2)) / sqrt(size(cued_proj_squezed, 2));
sem_vals_2 = squeeze(std(noncued_proj_squezed, 0, 3)) / sqrt(size(noncued_proj_squezed, 3));

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



% Save out pdf.
inches_out = [9 6];
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1 1 inches_out]);          % on screen
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 inches_out]);     % for printed/exported output
set(gcf, 'PaperSize', inches_out);        

% print(gcf, '-dpdf', 'encoding_??.pdf');
        

        %% CSE on single distractors - Heatmap                              
        figure

        timebnds = [100 400];
        timewin = dsearchn(times',timebnds(1)) : dsearchn(times',timebnds(2));

        d = mean(ENCadapt_htmp(:,:,timewin,:),[3,4],'omitnan');
        heatmap({'C','S','E','P'},{'C','S','E','P'},d')
        set(gca,'fontsize',14)
        colormap turbo
        ylabel('Current distractor')
        xlabel('Previous distractor (congruent-incongruent)')

        customBlueColormap = [linspace(0, 1, 256)', linspace(0, 1, 256)', ones(256, 1)];
        colormap(customBlueColormap)

        % clim([-0.1 0.05])

        
 
        %% CSE on single distractors ON vs OFF diagonal

 
        fontsize = 18;
        % for poster
        fontsize = 22;

        
        % CSE on single distractors ON vs OFF diagonal
        % Compute the mean and SEM
        mean_vals = squeeze(mean(ENCadaptations, 3)); 
        sem_vals = squeeze(std(ENCadaptations, 0, 3)) / sqrt(size(ENCadaptations, 3));
        
        % Create the figure and plot
        figure
        % sgtitle('CSE on distractors ON vs OFF diagonal', 'FontSize', 20)
        
        % Colors for the lines and shaded areas
        colors = [0 0 1; 0.5 0 1]; % Example colors: blue for ON, violet for OFF
        
        ylims = [-0.02 0.01]; 

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
        ylabel('\Delta', 'FontSize', fontsize)
        set(gca, 'FontSize', fontsize)
        legend('', 'ON diagonal', '', 'OFF diagonal', 'location', 'southeast')        % Add legend
        hold off
        
        CSE_permutest(ENCadaptations,times)  





        % CSE on single distractors ON diagonal, by distractor (diagonal elements)
        % Compute the mean and SEM
        mean_vals = squeeze(mean(ENCadaptations_on_bydist, 3)); 
        sem_vals = squeeze(std(ENCadaptations_on_bydist, 0, 3)) / sqrt(size(ENCadaptations_on_bydist, 3));
        
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
        ylabel('\Delta', 'FontSize', fontsize)
        set(gca, 'FontSize', fontsize)
        legend('', 'C', '', 'S', '', 'E', '', 'P', 'location', 'southeast')        % Add legend
        hold off

        
        % run permutest on each distractor's adaptation
        labels = {{'C'}, {'S'}, {'E'}, {'P'}};
        for didx = 1:4
            CSE_permutest(ENCadaptations_on_bydist(didx,:,:),times,labels{didx})
        end


        %% Target parametric adaptation


        % Parametric congruency sequence effect on (averaged) cued dimensions (spaghetti shift)
        % Compute the mean and SEM
        ENCadaptations_target = squeeze(mean(ENCadaptations_param_bytarget,3)); % average of cued dims.
        mean_vals = squeeze(mean(ENCadaptations_target, 3)); 
        sem_vals = squeeze(std(ENCadaptations_target, 0, 3)) / sqrt(size(ENCadaptations_target, 3));
        
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
        ylabel('Projection', 'FontSize', 18)
        set(gca, 'FontSize', 18)
        lgd = legend('', '0', '', '1', '', '2', '', '3', 'location', 'northwest');        % Add legend
        lgd.Title.String = '# previous congruent';
        lgd.Title.FontSize = 12;
        hold off


        %% Distractor parametric adaptation


        % Parametric congruency sequence effect on (averaged) cued dimensions (spaghetti shift)
        % Compute the mean and SEM
        mean_vals = squeeze(mean(ENCadaptations_param, 3)); 
        sem_vals = squeeze(std(ENCadaptations_param, 0, 3)) / sqrt(size(ENCadaptations_param, 3));
        
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
        ylabel('Projection', 'FontSize', 18)
        set(gca, 'FontSize', 18)
        lgd = legend('', '0', '', '1', '', '2', '', '3', 'location', 'northwest');        % Add legend
        lgd.Title.String = '# previous congruent';
        lgd.Title.FontSize = 12;
        hold off











        %% Target parametric adaptation - by RT bin (explored 03/31/25)



        % Parametric congruency sequence effect on (averaged) cued dimensions (spaghetti shift)
        % Compute the mean and SEM
        ENCadaptations_target = squeeze(mean(adapt_param_bytarget_byRT,3,'omitnan')); % fourth dimension is RT bin. Average of cued dims (some NANs in 3 sbjs/cued task due to decimation).
        mean_vals = squeeze(mean(ENCadaptations_target, 4)); 
        sem_vals = squeeze(std(ENCadaptations_target, 0, 4)) / sqrt(size(ENCadaptations_target, 4));
        
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
            ylim([-0.1 0.2])
            xlim([times(1) times(end)])
            xlabel('Time [ms]', 'FontSize', 18)
            ylabel('Encoding projection', 'FontSize', 18)
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
        mean_vals = squeeze(mean(ENCadaptations_b, 4)); 
        sem_vals = squeeze(std(ENCadaptations_b, 0, 4)) / sqrt(size(ENCadaptations_b, 4));
        
        % Create the figure and plot
        figure        
        % Colors for the lines and shaded areas
        colors = [0 0 1; 0.5 0 1]; % Example colors: blue for ON, violet for OFF
        
        ylims = [-0.03 0.03]; 

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
            ylabel('\Delta', 'FontSize', fontsize)
            set(gca, 'FontSize', fontsize)
            legend('', 'ON diagonal', '', 'OFF diagonal', 'location', 'southeast')        % Add legend
            hold off
        end



        %%%%%%%%%%

        % CSE on single distractors ON diagonal, by distractor (diagonal elements)
        % Compute the mean and SEM
        mean_vals = squeeze(mean(ENCadaptations_on_bydist_b, 4)); 
        sem_vals = squeeze(std(ENCadaptations_on_bydist_b, 0, 4)) / sqrt(size(ENCadaptations_on_bydist_b, 4));
        
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
            ylabel('\Delta', 'FontSize', fontsize)
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
        mean_vals = squeeze(mean(ENCadapts_b.(array_filenames{arf}), 4)); 
        sem_vals = squeeze(std(ENCadapts_b.(array_filenames{arf}), 0, 4)) / sqrt(size(ENCadapts_b.(array_filenames{arf}), 4));
        
        % Create the figure and plot
        figure        
        % Colors for the lines and shaded areas
        colors = [0 0 1; 0.5 0 1]; % Example colors: blue for ON, violet for OFF
        ylims = [-0.05 0.05]; 

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
            ylabel('\Delta', 'FontSize', fontsize)
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
        mean_vals = squeeze(mean(ENCadapts_on_bydist_b.(array_filenames{arf}), 4)); 
        sem_vals = squeeze(std(ENCadapts_on_bydist_b.(array_filenames{arf}), 0, 4)) / sqrt(size(ENCadapts_on_bydist_b.(array_filenames{arf}), 4));

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
            ylabel('\Delta', 'FontSize', fontsize)
            set(gca, 'FontSize', fontsize)
            legend('', 'C', '', 'S', '', 'E', '', 'P', 'location', 'southeast')        % Add legend
            sgtitle(titles{arf}, 'FontSize', 20)
            hold off
        end
end









%% Double difference ON vs OFF diagonal --- exploring any asymmetry in adaptaiton - not comparable with behavioral double-difference! (05/20/2026)

fontsize = 18;

% ENCadapts_b fields: bins x time x on/off x sbj
% double difference: (CI - CC) - (II - IC)
% = prevC_incongruent - prevC_congruent - prevI_incongruent + prevI_congruent

CI = ENCadapts_b.CS_CI; % bins x time x on/off x sbj
CC = ENCadapts_b.CS_CC;
II = ENCadapts_b.CS_II;
IC = ENCadapts_b.CS_IC;

DD = (CI - CC) - (II - IC); % bins x time x on/off x sbj

% average over bins
DD_avg = squeeze(mean(DD, 1)); % time x on/off x sbj

mean_vals = squeeze(mean(DD_avg, 3)); % time x on/off
sem_vals = squeeze(std(DD_avg, 0, 3)) / sqrt(size(DD_avg, 3));

colors = [0 0 1; 0.5 0 1];

figure
hold on
for i = 1:2
    fill([times, fliplr(times)], ...
        [mean_vals(:,i)' - sem_vals(:,i)', fliplr(mean_vals(:,i)' + sem_vals(:,i)')], ...
        colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(times, mean_vals(:,i), 'LineWidth', 2, 'Color', colors(i,:));
end
plot(times, zeros(length(times),1), 'k:', 'LineWidth', 2)
xlim([times(1) times(end)])
xlabel('Time [ms]', 'FontSize', fontsize)
ylabel('Double difference (CI-CC)-(II-IC)', 'FontSize', fontsize)
title('Asymmetry in adaptation? Neural dobule difference - On vs Off diagonal', 'FontSize', fontsize)
legend('', 'ON diagonal', '', 'OFF diagonal', 'location', 'southeast')
set(gca, 'FontSize', fontsize)
hold off
