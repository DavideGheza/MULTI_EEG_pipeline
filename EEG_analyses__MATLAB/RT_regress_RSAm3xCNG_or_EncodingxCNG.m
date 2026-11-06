% Fit trial-level RT models from RSA or encoding estimates.
% Intended for HPC execution. `RSA = 1` loads model-3 RSA estimates;
% `RSA = 0` loads full-model encoding projections.
% Fitting model (multilevel), at each time point:
%   RT ~ cued*congruency + non-cued average*congruency (parametric) + task*side
% Fitted separately post- low vs high- congruency (0,1 vs. 2,3 congruent
% dimensions on the previous trial), and on all trials.
% Author: Davide Gheza


function reg = RT_regress_RSAm3xCNG_or_EncodingxCNG

RSA=1;
local = 0;

% Select stimulus- or response-locked inputs.
onset = true;
if onset
    suffix = '';
else
    suffix = 'response_';   % response-locked
end

if local==1 
    if ispc
        load('C:\Users\gdavide\OneDrive\Documenti\MATLAB\dg_workspace22\dimensions\analyses\scripts_and_settings\preprocessing_settings.mat');
    else
        load('/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/scripts_and_settings/preprocessing_settings_mac.mat');
    end
end

data = [];  % initialize data table
sn=[22 24:48];
for si = 1:length(sn)
    
    % load behavioral data
    if local==0
        behav = load(['/scratch/gheza/EEG_data/behav/V3_behav_seqinfo_trimmed_' suffix 'sbj_' num2str(sn(si))]);
    else
        behav = load([sets.sessions_data_dir 'V3_behav_seqinfo_trimmed_' suffix 'sbj_' num2str(sn(si))]); % behavioral file with addidional info on subsequent and previous trials (10/01/2024)
    end
    behav = behav.behav_trim;
    task = behav.task;
    
    if RSA==0
        % load single trial projections
        if local==0
            load(['/scratch/gheza/RSA_fits/encoding_dssr_' suffix 'sbj_' num2str(sn(si)) '.mat']);
        else
            load(fullfile(sets.rsa_data_dir, ['encoding_model/encoding_dssr_' suffix 'sbj_' num2str(sn(si)) '.mat'])) % = all single trial data    (trials, ses, times)
        end
        projections = result.encoded_dim_spec_SR;
        projections = cat(4,projections(1).s,... % cat the 4 PCs 
                                       projections(2).s,...
                                       projections(3).s,...
                                       projections(4).s);
        projections = permute(projections,[1,4,2,3]); % permute as RSA fits by single-distractor: trial, pc, time, dimension
        % parse cued from non-cued projections
        noncued_proj = projections;
        cued_proj = nan(size(noncued_proj,(1:3)));
        for t = 1:size(noncued_proj,1)
            cued_proj(t,:,:) = projections(t,:,:,task(t));  
            noncued_proj(t,:,:,task(t)) = nan;
        end

        % organize data in long format - cued and each-non cued dims
        reshaped.cued = reshape(cued_proj, [], 1);
        reshape_field = @(field) reshape(squeeze(field(:,:,:,1)), [], 1);
        fnames = {'nc_1','nc_2','nc_3','nc_4'};
        for i = 1:length(fnames)
            reshaped.(fnames{i}) = reshape_field(noncued_proj(:,:,:,i));
        end

    else
        % load single trial RSA (m3 QR)
        if local==0
            load(['/scratch/gheza/RSA_fits/fits_singletrials_m3_sbj_' num2str(sn(si)) '.mat']);
        else
            load(fullfile(sets.rsa_data_dir, ['fits/fits_singletrials_m3_sbj_' num2str(sn(si)) '.mat'])) % fits = all single trial data    (trials, ses, times)
        end
        cued_proj = fits.slope_sr; % rename just to size(cued_proj) later
        % organize data in long format
        reshaped.cued = reshape(fits.slope_sr,[],1); % here proj is misnomer for RSA cued SR, but using it for consistency (adapting script built for encoding results)
        slope_sr_singledstr_1 = squeeze(fits.slope_sr_singledstr(:,:,:,1));
        reshaped.nc_1 = reshape(slope_sr_singledstr_1,[],1);
        slope_sr_singledstr_2 = squeeze(fits.slope_sr_singledstr(:,:,:,2));
        reshaped.nc_2 = reshape(slope_sr_singledstr_2,[],1);
        slope_sr_singledstr_3 = squeeze(fits.slope_sr_singledstr(:,:,:,3));
        reshaped.nc_3 = reshape(slope_sr_singledstr_3,[],1);
        slope_sr_singledstr_4 = squeeze(fits.slope_sr_singledstr(:,:,:,4));
        reshaped.nc_4 = reshape(slope_sr_singledstr_4,[],1);
        reshaped.dstr_avg_original = reshape(fits.slope_sr_dstr,[],1);% storing average of singledstr, to compare to average done below
    end




    trials_to_reshape = repmat(1:size(cued_proj,1), [size(cued_proj,2),size(cued_proj,3)]);
    trials = reshape(trials_to_reshape',[],1);
    ses_to_reshape = repmat(1:size(cued_proj,2), [size(cued_proj,1),size(cued_proj,3)]);
    ses = reshape(ses_to_reshape,1,[]);
    ses = reshape(ses,[],1);
    times_to_reshape = repmat(1:size(cued_proj,3), [size(cued_proj,1)*size(cued_proj,2),1]);
    times = reshape(times_to_reshape,[],1);
    % log tranform and prewhiten RTs (linear and quadratic trends over trials)
    RT = log(behav.rt);  % log transform
    n = length(RT);
    X = [(1:n)', (1:n)'.^2];
    mdl = fitlm(X, RT, 'linear');
    RT = mdl.Residuals.Raw; % take residuals
    % repeat and reshape RT like trial index above
    RT = repmat(RT', [size(cued_proj,2),size(cued_proj,3)]);
    RT = reshape(RT',[],1);
    % add sbj
    sbj = repmat(si, [size(cued_proj,1)*size(cued_proj,2)*size(cued_proj,3),1]);
    % convert ses to PC ordinal (1-3 from MI to Neu to MC)
    PCo = nan(size(ses));
    PCo(ses==1|ses==4)=2; % Neu 1 and 2
    PCo(ses==2)=1;        % MI
    PCo(ses==3)=3;        % MC
    % add congruency level
    cng = sum(behav.congruency==1,2);  
    % repeat and reshape cng like trial index above
    cng = repmat(cng', [size(cued_proj,2),size(cued_proj,3)]);
    cng = reshape(cng',[],1);
    % add previous congruency level
    prev_cng = behav.prev_sum_con;  
    % repeat and reshape cng like trial index above
    prev_cng = repmat(prev_cng', [size(cued_proj,2),size(cued_proj,3)]);
    reshaped.prev_cng = reshape(prev_cng',[],1);

    % add filtering factors
    reshaping = repmat(behav.taskcount', [size(cued_proj,2),size(cued_proj,3)]);
    reshaped.taskcount = reshape(reshaping',[],1);
    reshaping = repmat(behav.post_correct', [size(cued_proj,2),size(cued_proj,3)]);
    reshaped.post_correct = reshape(reshaping',[],1);
    % add info about cued task and side of response (to model out RT difference by dimension/side)
    reshaping = repmat(behav.task', [size(cued_proj,2),size(cued_proj,3)]);
    reshaped.task = reshape(reshaping',[],1);
    reshaping = repmat(behav.response', [size(cued_proj,2),size(cued_proj,3)]); % coding for response side. Equivalent to labeling used in RSA fit. To see: behav.response(behav.task==1) matches behav.colors(behav.task==1)
    reshaping = ([reshaping(:)]=='j') + 1;
    reshaped.side = reshape(reshaping',[],1);
    % store in table
    all_fnames = fieldnames(reshaped);
    data_s = [];
    for k = 1:length(all_fnames)
        data_s = [data_s, reshaped.(all_fnames{k})];
    end
    dstr_avg = mean(data_s(:, 2:5), 2, 'omitnan');
    data_cols = [data_s, dstr_avg, sbj, trials, ses, times, RT, PCo, cng];
    tbl = array2table(data_cols, 'VariableNames', {all_fnames{:}, 'dstr_avg', 'sbj', 'trials', 'ses', 'times', 'RT', 'PCo', 'cng'});
    % append subject-specific tables into one
    data = [data; tbl];
    fprintf('\tLoaded sbj %g\n', si) 

end

data.task = categorical(data.task);
data.side = categorical(data.side);
    
    % preallocate
    len = length(unique(data.times));
    
    % all field names to be initialized
    % fnames = {'intercept', 'cued', 'dstr_avg', 'parm_cng', 'cued_and_cng', 'dstr_avg_and_cng'};
    fnames = {'intercept', 'cued', 'task_2', 'task_3', 'task_4', 'side_2', 'dstr_avg', 'parm_cng', 'task_2and_side_2', 'task_3and_side_2', 'task_4and_side_2', 'cued_and_cng', 'dstr_avg_and_cng'};
    
    % Initialize the structure reg with NaN vectors
    reg = struct();
    for i = 1:length(fnames)
        reg.(fnames{i}) = nan(len, 3); % for betas % 2nd dimension for storing post_low vs. high cng, or both
        fn = [fnames{i},'_t'];
        reg.(fn) = nan(len, 3); % for tStat
    end
    
    
    
        % filter trials: given we model sequential effects, let's drop post-incorrect and post-switch trials
        include = data.taskcount>1 & data.post_correct==1;
        data = data(include,:);
        % filter out NaN data - if cued proj is NaN, no data for that PC (mind that each PC is applied to some trials only - NaN used as filler)
        data = data(~isnan(data.cued),:); 

        store.post_inc = data(data.prev_cng < 2, :);
        store.post_con = data(data.prev_cng >= 2, :);
        store.all = data;
        store_fn = fieldnames(store);
        for prev_cng = 1:3 
            data = store.(store_fn{prev_cng});
            
            for t = 1:length(unique(data.times))    
                timesdata = data(data.times==t,:); % filter single time data
                % predictors
                lm = fitlme(timesdata, 'RT~ cued*cng + dstr_avg*cng + task*side + ( 1 + cued + dstr_avg + cng + task*side | sbj)');
                % store beta and t-values
                for b_idx = 1:length(fnames)
                fn = fnames{b_idx};
                reg.(fn)(t,prev_cng) = lm.Coefficients.Estimate(b_idx);
    
                fn = [fnames{b_idx},'_t'];
                reg.(fn)(t,prev_cng) = lm.Coefficients.tStat(b_idx);
                end
         
        
                fprintf('\tFitted t %g\n', t)    
            end
        end
    
                
                
    if RSA
        fielname = sprintf('RT_RSA_cued_distAVG_XparamCNG_m3_nonQR');
    else
        fielname = sprintf('RT_Encoding_cued_distAVG_XparamCNG');
    end
    if local == 0
        fname = sprintf(['/scratch/gheza/RSA_fits/' fielname '.mat']);
    else
        fname = sprintf([sets.rsa_data_dir, 'fits/' fielname '.mat']);
    end

    save(fname, 'reg');



end
