%% Analyze the relationship between time-frequency power and RT
% Loads single-trial normalized power, aligns behavioral RTs, and computes
% or visualizes the selected stimulus- or response-locked association.
% Author: Davide Gheza

clc
clear

%% 1. Load settings

freqshown = [3 7];
singletrial = true;
rt_corr=1;
%Select onset or response events
onset = true;
if onset
    suffix_eeg = 'onset_';
    suffix = '';
else
    suffix = 'response_';   % response-locked requires precomputing stimulus-locked to get a fixed pre-stimulus baseline
    suffix_eeg = 'response_';
end

%Struct 'sets' is loaded from the 'scripts_and_settings' folder:
if ispc
    load('C:\Users\gdavide\OneDrive\Documenti\MATLAB\dg_workspace22\dimensions\analyses\scripts_and_settings\preprocessing_settings.mat');
else
    load('/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/scripts_and_settings/preprocessing_settings_mac.mat');
end

% load an eeg dataset to fetch info
sampleEEG = pop_loadset('filename',['22' suffix_eeg 'trialrej.set'],...
                'filepath',[sets.clean_data_dir]...
                );
                
chanlocs = sampleEEG.chanlocs;
labels = {sampleEEG.chanlocs.labels};

%% Load normalized power 
sn=[22 24:48]; 

tf_times = (-800:10:1490); % hardcoded WARNING - as from t/f 2025 settings
if onset
    times_regression = [-50 700]; % for stim-locked data
else
    times_regression = [-700 50]; % for response-locked data. Useless to go much beyond response (no conceptual sense)
end
t_idx_reg = dsearchn(tf_times',times_regression(1)):dsearchn(tf_times',times_regression(2));

% sbj, chans, freqs, times
power = zeros(length(sn),61,30,length(t_idx_reg));
power_si = zeros(61,30,length(t_idx_reg));
intrcpt = zeros(length(sn),61,30,length(t_idx_reg));
intrcpt_si = zeros(61,30,length(t_idx_reg));

for si = 1:length(sn)%(sets.sub_idx)
    
    load([sets.sessions_data_tf_dir...
            'data_' suffix 'tf_singletrial_NB_sbj_'...
            sets.sub_idx{sn(si)}...
            '.mat']);
    
    % load behavioral data
    behav = load([sets.sessions_data_dir '/V3_behav_seqinfo_trimmed_' suffix 'sbj_' num2str(sn(si))]);
    behav = behav.behav_trim;
    


    % RT-power correlations

    power_temp = data.ses1(:,:,:,t_idx_reg); % subset times
    rt_selection = behav.rt;
    power_temp = zscore(power_temp,[],[1,3]); % z score over frequencies, and over trials (comparison of predictive power across frequencies - which frequencies are most predictive of RT, independent of amplitude scaling)
    % power_temp = zscore(power_temp,[],1); % z score over trials only (NOT TESTED - for a more physiological comparison of bands, but not directly comparable bands)
    rt_selection = zscore(rt_selection,[],1); % z score over trials
    for c = 1:size(power_temp,2)
        for f = 1:size(power_temp,3)
            for t = 1:size(power_temp,4)
                lm = fitlm(power_temp(:,c,f,t),rt_selection,'RobustOpts','on'); % fit linear model, for each channel, frequency, time - variance over trials per condition
                slope = lm.Coefficients.Estimate(2);
                intercept = lm.Coefficients.Estimate(1);
                power_si(c,f,t) = slope;
                intrcpt_si(c,f,t) = intercept;
            end
            fprintf('\tLMed sbj %g chan %g freq %g\n', si, c, f)                      
        end
    end

        power(si,:,:,:) = power_si; % catenate sbj data along 1st dim (sbj, chan, freq, time) - NB correlations run on power without baseline normalization 
        intrcpt(si,:,:,:) = intrcpt_si; 
fprintf('\tLoaded sbj %g / %g\n', si, length(sn))
end

cd(sets.sessions_data_tf_dir)

if onset
    fname1 = ('rt_power_regression_stimlocked2025_slope');
    fname2 = ('rt_power_regression_stimlocked2025_intercept');
else
    fname1 = ('rt_power_regression_resplocked2025_slope');
    fname2 = ('rt_power_regression_resplocked2025_intercept');
end

save (fname1, 'power')
save (fname2, 'intrcpt')


