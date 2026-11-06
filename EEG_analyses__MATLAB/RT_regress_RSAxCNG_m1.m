% Fit trial-level RT models from model-1 RSA estimates.
% Intended for HPC execution. Fits cued stimulus-response, rule, and
% dimension-general side estimates across participants, including their
% interactions with congruency. Non-cued RSA estimates are not modeled.
% Author: Davide Gheza


function reg = RT_regress_RSAxCNG_m1

local = 0;
model = 1;

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

    % load single trial RSA fits
    if local==0
        load(['/scratch/gheza/RSA_fits/fits_singletrials_m' num2str(model) '_sbj_' num2str(sn(si)) '.mat']);
    else
        load(fullfile(sets.rsa_data_dir, ['fits/fits_singletrials_m' num2str(model) '_sbj_' num2str(sn(si)) '.mat'])) % fits = all single trial data    (trials, ses, times)
    end
    % organize data in long format
    reshape_field = @(field) reshape(field, [], 1);
    fnames = {'slope_sr', 'slope_ru', 'slope_si'};
    reshaped = struct();
    for i = 1:length(fnames)
        reshaped.(fnames{i}) = reshape_field(fits.(fnames{i}));
    end
    reshape_field = @(field) reshape(squeeze(field(:,:,:,1)), [], 1);
    fnames = {'slope_ru_singledstr_1','slope_ru_singledstr_2','slope_ru_singledstr_3','slope_ru_singledstr_4'};
    for i = 1:length(fnames)
        reshaped.(fnames{i}) = reshape_field(fits.slope_ru_singledstr(:,:,:,i));
    end


    trials_to_reshape = repmat(1:size(fits.slope_sr,1), [size(fits.slope_sr,2),size(fits.slope_sr,3)]);
    trials = reshape(trials_to_reshape',[],1);
    ses_to_reshape = repmat(1:size(fits.slope_sr,2), [size(fits.slope_sr,1),size(fits.slope_sr,3)]);
    ses = reshape(ses_to_reshape,1,[]);
    ses = reshape(ses,[],1);
    times_to_reshape = repmat(1:size(fits.slope_sr,3), [size(fits.slope_sr,1)*size(fits.slope_sr,2),1]);
    times = reshape(times_to_reshape,[],1);

    % log tranform and prewhiten RTs (linear and quadratic trends over trials)
    RT = log(behav.rt);  % log transform
    n = length(RT);
    X = [(1:n)', (1:n)'.^2];
    mdl = fitlm(X, RT, 'linear');
    RT = mdl.Residuals.Raw; % take residuals
    % repeat and reshape RT like trial index above
    RT = repmat(RT', [size(fits.slope_sr,2),size(fits.slope_sr,3)]);
    RT = reshape(RT',[],1);
    % add sbj
    sbj = repmat(si, [size(fits.slope_sr,1)*size(fits.slope_sr,2)*size(fits.slope_sr,3),1]);
    % convert ses to PC ordinal (1-3 from MI to Neu to MC)
    PCo = nan(size(ses));
    PCo(ses==1|ses==4)=2; % Neu 1 and 2
    PCo(ses==2)=1;        % MI
    PCo(ses==3)=3;        % MC
    % add congruency level
    cng = sum(behav.congruency==1,2);  
    % repeat and reshape cng like trial index above
    cng = repmat(cng', [size(fits.slope_sr,2),size(fits.slope_sr,3)]);
    cng = reshape(cng',[],1);
    % add info about cued task and side of response (to model out RT difference by dimension/side)
    reshaping = repmat(behav.task', [size(fits.slope_sr,2),size(fits.slope_sr,3)]);
    reshaped.task = reshape(reshaping',[],1);
    reshaping = repmat(behav.response', [size(fits.slope_sr,2),size(fits.slope_sr,3)]); % coding for response side. Equivalent to labeling used in RSA fit. To see: behav.response(behav.task==1) matches behav.colors(behav.task==1)
    reshaping = ([reshaping(:)]=='j') + 1;
    reshaped.side = reshape(reshaping',[],1);
    % store in table
    all_fnames = fieldnames(reshaped);
    data_s = [];
    for k = 1:length(all_fnames)
        data_s = [data_s, reshaped.(all_fnames{k})];
    end
    ru_dstr_avg = mean(data_s(:, 4:7), 2, 'omitnan');
    data_cols = [data_s, ru_dstr_avg, sbj, trials, ses, times, RT, PCo, cng];
    tbl = array2table(data_cols, 'VariableNames', {all_fnames{:}, 'ru_dstr_avg', 'sbj', 'trials', 'ses', 'times', 'RT', 'PCo', 'cng'});
    % append subject-specific tables into one
    data = [data; tbl];
    fprintf('\tLoaded sbj %g\n', si) 

end

data.task = categorical(data.task);
data.side = categorical(data.side);

% preallocate
len = length(unique(data.times));

% all field names to be initialized
fnames = {'intercept', 'sr', 'ru', 'si', 'task2', 'task3', 'task4', 'side2', 'cng', 'task2Xside2', 'task3Xside2', 'task4Xside2', 'slope_srXcng', 'slope_ruXcng', 'slope_siXcng'};

% Initialize the structure reg with NaN vectors
reg = struct();
for i = 1:length(fnames)
    reg.(fnames{i}) = nan(len, 1);
    fn = [fnames{i},'_t'];
    reg.(fn) = nan(len, 1); % for tStat
end




    
    data = data(~isnan(data.slope_sr),:); % filter out NaN data - check with SR model fits (mind that each PC is applied to some trials only - NaN used as filler)

    for t = 1:length(unique(data.times))    
        timesdata = data(data.times==t,:); % filter single time data
        % predictors are cued dimension fits
        lm = fitlme(timesdata, 'RT~ slope_sr*cng + slope_ru*cng + slope_si*cng + task*side + ( 1 + slope_sr*cng + slope_ru*cng + slope_si*cng + task*side | sbj)'); % individual intercepts and slopes
        % store beta and t-values
        for b_idx = 1:length(fnames)
        fn = fnames{b_idx};
        reg.(fn)(t) = lm.Coefficients.Estimate(b_idx);

        fn = [fnames{b_idx},'_t'];
        reg.(fn)(t) = lm.Coefficients.tStat(b_idx);
        end        


        fprintf('\tFitted t %g\n', t)    
    end

            
            



    fielname = sprintf(['RT_RSAxCNG_M' num2str(model) '_cued.mat']);

    if local == 0
        fname = sprintf(['/scratch/gheza/RSA_fits/' fielname '.mat']);
    else
        fname = sprintf([sets.rsa_data_dir, 'fits/' fielname '.mat']);
    end

    save(fname, 'reg');



end
