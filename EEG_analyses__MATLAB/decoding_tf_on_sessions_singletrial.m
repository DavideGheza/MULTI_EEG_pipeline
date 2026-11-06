%% Compute and save single-trial time-frequency power by session
% Produces the single-trial power arrays used by CSE and RT-regression
% analyses.
% Author: Davide Gheza

%% 1. Load settings
clear

%Select onset or response events
onset = false;
if onset
    suffix = '';
else
    suffix = 'response_';   % response-locked requires precomputing stimulus-locked to get a fixed pre-stimulus baseline
end

%Struct 'sets' is loaded from the 'scripts_and_settings' folder:
if ispc
    load('C:\Users\gdavide\OneDrive\Documenti\MATLAB\dg_workspace22\dimensions\analyses\scripts_and_settings\preprocessing_settings.mat');
else
    load('/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/scripts_and_settings/preprocessing_settings_mac.mat');
end

%Add path to sub-folder of 'main_dir' to store datasets:
sets.sessions_data_tf_dir=[sets.main_saving_dir 'sessions_data_tf/'];
%If the folder doesn't exist yet, create one:
if ~exist(sets.sessions_data_tf_dir,'dir')
    mkdir(sets.sessions_data_tf_dir);
    %Save updated 'sets' struct:
    if ispc
        save([sets.settings_dir 'preprocessing_settings.mat'],'sets');
    else
        save([sets.settings_dir 'preprocessing_settings_mac.mat'],'sets');
    end
end



root = sets.sessions_data_dir; % data

% behavioral data file path
cd (sets.raw_dir);
sbj = dir('*PC*');

%% 2. t-f convolution

% specify params:
% wavelet parameters
min_freq = 1.5;
max_freq = 40;
num_frex = 30;             
wavelet_cycles_lims = [2,10];
% define baseline period            consider not applying
baselinetime = [ -800 -500 ]; % in ms                                      % 07/01/2022 task V3: task-cue onset = -200  
% downsampled times
times2save = -800:10:1498;  % from 5 to 10 to save space                   % 07/01/2022 task V3: increased iti range: 1300-1500 

% EEG data info
EEG.srate = 500;
EEG.times = (-950:1000/EEG.srate:2049); 


% Loop over participant datasets.
for si=[22,24:length(sets.sub_idx)]
    
    
    load(fullfile([root 'sessions_data_' suffix 'sbj_' num2str(sets.sub_idx{si}) '.mat']));
    dn = fieldnames(data);
    if ~onset   % if response locked: load baseline computed for stimulus-locked 
        stimulus_locked_tf = load(fullfile([sets.sessions_data_tf_dir 'data_tf_singletrial_NB_sbj_' num2str(sets.sub_idx{si}) '.mat']));
        baseline_sl = stimulus_locked_tf.data.baseline;
        stimulus_locked_tf = [];
    end
    
    
    % other wavelet parameters
    frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
    time = -2:1/EEG.srate:2;    %define length in time of wavelet: change according to desired time window (e.g. -2 2 sec) - keep in ceneterd, odd lenght in datapoints, sampling rate of EEG data
    half_of_wavelet_size = (length(time)-1)/2;
    s = logspace(log10(wavelet_cycles_lims(1)),log10(wavelet_cycles_lims(2)),num_frex)./(2*pi*frequencies);          % s = std of gaussian: wavelet_cycles /(2*pi*f) --> instead of fixed n cycles, here define a vector of cycles n, varying across frex with a log scaling
    
    for se = 1:numel(dn)/3                                                           % loop through sessions (only 1 in task V3) (use first fields in data struct)
        % session specific EEG data info
        EEG.data = data.(dn{se});
        EEG.pnts = size(EEG.data,3);
        EEG.trials = size(EEG.data,1);
        EEG.nbchan = size(EEG.data,2);
        EEG.data = permute(EEG.data,[2,3,1]);        %  chans, times, trials - back to EEGlab order
        
        % FFT parameters (use next-power-of-2)
        n_wavelet     = length(time);
        n_data        = EEG.pnts*EEG.trials;                                    %if avoid chan loop and concatenate chans: n_data        = EEG.pnts*EEG.trials*EEG.nbchan;
        n_convolution = n_wavelet+n_data-1;
        n_conv_pow2   = pow2(nextpow2(n_convolution));

        % convert downsampled times to indices
        times2save_idx = zeros(size(times2save));
        for i=1:length(times2save)
            [junk,times2save_idx(i)] = min(abs(EEG.times-times2save(i)));
        end
        % convert baseline window time to indices (based on times2save)
        [~,baselineidx(1)] = min(abs(times2save-baselinetime(1)));
        [~,baselineidx(2)] = min(abs(times2save-baselinetime(2)));

        % allocate tf matrix (chan, freq, time, trial)
        tf_data = zeros(EEG.nbchan,length(frequencies),EEG.pnts,EEG.trials,'single');               % single precision to save ram and hdd
        tf_data_sub = zeros(EEG.nbchan,length(frequencies),numel(times2save),EEG.trials,'single');
        baseline = zeros(EEG.nbchan,length(frequencies),EEG.trials);                                % save single-trial baseline (19 Sept 2025)

        % fft of data 
        % originally avoided loop over chans ----------------> it is not a bad idea to concatenate chans, but requires a lot of RAM
        for chan_x = 1:EEG.nbchan
            fft_data = fft(reshape(EEG.data(chan_x,:,:),1,[]),n_conv_pow2);     % EEG.data: chans,times,trials. Reshape to 1 timeserie (concatenate times, trials). Apply FFT with zero-padding to match n_conv_pow2

            for fi=1:length(frequencies)

                % create wavelet and get its FFT
            %     wavelet = exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi);     %(pi*frequencies(fi)*sqrt(pi))^-.5 scaling factor
                wavelet = exp(2*1i*pi*frequencies(fi).*time) .* ...          % complex sine wave
                          exp(-time.^2./(2*(s(fi))^2));                      % gaussian - std depends on frequency (defined on varying number of cycles)

                fft_wavelet = fft(wavelet,n_conv_pow2);

                % run convolution
                convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
                convolution_result_fft = convolution_result_fft(1:n_convolution);                                   % cut out 0-padding from FFT function (applied to match n_conv_pow2)
                convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);   % cut out 0-padding for convolution of entire signal (applied to match n_convolution)
                convolution_result_fft = reshape(convolution_result_fft,EEG.pnts,EEG.trials);                       % if you want to avoid chan loop, here reshape including 1st dim EEG.nbchan


                % put power data into time-frequency matrix (chan, freq, time, trial)
        %        tf_data(chan_x,fi,:,:) = convolution_result_fft;            % save complex numbers (include power and phase)
                tf_data(chan_x,fi,:,:) = abs(convolution_result_fft).^2;  % power only
        %         tf_data(2,fi,:) = median(abs(convolution_result_fft).^2,2);

            % show counter
            formatSpec = 'done with sbj %d, session %d, chan % 2.0f, freq %2.0f corresponding to %4.2f Hz\n';
            fprintf(formatSpec, si, se, chan_x, fi, frequencies(fi)); 

            end
        end
        % subsample times
        tf_data_sub(:,:,:,:) = tf_data(:,:,times2save_idx,:);
        % save pre-cue baseline (chan, freq, trials)       
        if onset
            baseline(:,:,:) = mean(tf_data_sub(:,:,(baselineidx(1):baselineidx(2)),:),3);        % mean power across times, KEEP trials
        else    % response locked: use baseline computed for stimulus-locked (pre-stumulus)
            baseline = baseline_sl;
        end

            
        % permute
        tf_data_sub = permute(tf_data_sub,[4,1,2,3]);      %  trials, chans, freq, time - back to RSA order ---> saving tf_data without normalization
        baseline = permute(baseline,[3,1,2]);              % (trials, chans, freq) 

    % save session tf
    data.(dn{se}) = tf_data_sub;
    end    

    data.times = times2save;
    data.freqs = frequencies;
    data.baseline = baseline;
    
    
    cd(sets.sessions_data_tf_dir)
%     fname = sprintf('sessions_data_tf_decoding_NB_Z_sbj_%d.mat', si);
    fname = sprintf(['data_' suffix 'tf_singletrial_NB_sbj_%d.mat'], si); % baseline normalization not yet applied
    save(fname, 'data', '-v7.3')    
end
