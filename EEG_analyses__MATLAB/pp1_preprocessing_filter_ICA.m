%% Preprocessing of 'Dimensions' experiment EEG data
% Imports BrainVision `.vhdr`, `.eeg`, and `.vmrk` data, applies the continuous
% EEG filtering pipeline, and estimates ICA weights. Requires EEGLAB, ERPLAB,
% and the EEGLAB BVA-IO plugin.
% Author: Davide Gheza

%% 1. Settings
%
% All preprocessing parameters and paths are stored in the `sets` structure.


% --- Paths and directories ---
%
%1-'exp_name': [char] experiment name.
sets.exp_name='dimensions';

%2-'sub_num': [double] number of experimental subjects.
sets.sub_num=48;                                                           % note:   1:8 pilots 1st version.    9:21 pilots 2nd version.    22:end 3rd version (reported experiment).

%3-'settings_dir': [char] directory for loading and saving various files used in preprocessing (e.g. channels location
%file). A copy of struct 'sets' is saved in this direcctory, as .mat file.
if ispc
    sets.settings_dir='C:\Users\gdavide\OneDrive\Documenti\MATLAB\dg_workspace22\dimensions\analyses\scripts_and_settings\';
else
    sets.settings_dir='/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/scripts_and_settings/';
end
%4-'raw_dir': [char] location of raw eeg data directory.
if ispc
    sets.raw_dir='C:\Users\gdavide\OneDrive\Documenti\MATLAB\dg_workspace22\dimensions\data\EEG_data';
else
    sets.raw_dir='/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/data/EEG_data';
end

%5-'main_saving_dir': [char] main directory for saving preprocessed data.
if ispc
    sets.main_saving_dir='C:\Users\gdavide\OneDrive\Documenti\MATLAB\dg_workspace22\dimensions\analyses\preprocessed_data\';
else
    sets.main_saving_dir='/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/preprocessed_data/';
end

% %6-'raw_data_dir': [char] subfolder of 'main_saving_dir', to store raw data as .set and .fdt files.
% sets.raw_data_dir=[sets.main_saving_dir 'raw_data_dir\'];
% %If the folder doesn't exist yet, create one:
% if ~exist(sets.raw_data_dir,'dir')
%     mkdir(sets.raw_data_dir);
% end

%7-'continuous_data_dir': [char] subfolder of 'main_saving_dir', to store and load conitnuous filtered data (data right
%before epoching).
sets.continuous_data_dir=[sets.main_saving_dir 'continuous_filtered_data/'];
%If the folder doesn't exist yet, create one:
if ~exist(sets.continuous_data_dir,'dir')
    mkdir(sets.continuous_data_dir);
end

%8-'epoch_data_dir': [char] subfolder of 'main_saving_dir', to store and load epoched data.
sets.epoch_data_dir=[sets.main_saving_dir 'epoch_data/'];
%If the folder doesn't exist yet, create one:
if ~exist(sets.epoch_data_dir,'dir')
    mkdir(sets.epoch_data_dir);
end


% --- Participants, conditions and triggers ---

%9-'sub_idx': [cell array] subjects indexes, each index stored in a different cell as string. Strings correpond to names
%given to raw .vhdr files.
%In this specific case, file names are simply numbers corresponding to participants' IDs (from 1 to 44, e.g.: '01.bdf'). 
%So they can be created in a loop. 
sets.sub_idx=cell(1,sets.sub_num);
count=1;
for i= 1:sets.sub_num %
%     if i<10
%         sets.sub_idx{count}=['0' num2str(i)];
%     else
        sets.sub_idx{count}=num2str(i);
%     end
    count=count+1;
end

%10-trigger categories: [struct] sub-struct containing 2 fields:
%   -names: [char] name of condition:
%   -markers: [double] corresponding trigger value
sets.taskcue.names={'tcC', 'tcS', 'tcE', 'tcP'};                               % task-cue markers (181, 182, 183, 184) 
sets.taskcue.markers={'S181' 'S182' 'S183' 'S184'};

sets.onset_dim.names={'C', 'S', 'E', 'P'};                            % Task-cue exceptions: sbj 22 used O ('orientation') instead of P (Pattern); sbj 22 and 23 used B ('border') instead of E (Edge).
sets.onset_dim.markers={{'S  1', 'S 11', 'S 51', 'S101', 'S 61', 'S151', 'S111', 'S161'},...
                        {'S  2', 'S  7', 'S 52', 'S102', 'S 57', 'S107', 'S152', 'S157'},...
                        {'S  3', 'S  8', 'S 13', 'S103', 'S 18', 'S108', 'S113', 'S118'},...
                        {'S  4', 'S  9', 'S 14', 'S 54', 'S 19', 'S 59', 'S 64', 'S 69'}};
                
sets.onset_con.names={'con3', 'con2', 'con1', 'con0'};
sets.onset_con.markers={{'S  1', 'S  2', 'S  3', 'S  4'},...
                        {'S 11', 'S 51', 'S101', 'S  7', 'S 52', 'S102', 'S  8', 'S 13', 'S103', 'S  9', 'S 14', 'S 54'},...
                        {'S 61', 'S151', 'S111', 'S 57', 'S107', 'S152', 'S 18', 'S108', 'S113', 'S 19', 'S 59', 'S 64'},...
                        {'S161', 'S157', 'S118', 'S 69'}};




sets.resp_pc.names={'MI', 'MC', 'N'};
sets.resp_pc.markers={'S 71' 'S 72' 'S 75'};


% --- Preprocessing parameters ---

%14-'remove_chan': [cell array] channels to remove before preprocessing. One channel for cell, as string;
sets.remove_chan={};

%15-'filters': [struct] filters' parameters. Each filter is a different entry in the struct. Filters must be entered IN
% Store filters in execution order so they can be applied in one loop.
%Fields:
%  -'locutoff': [double] lower passband edge frequency (if [] or 0 -> lowpass filter);
%  -'hicutoff': [double] higher passband edge frequency (if [] or 0 -> highpass filter);
%  -'order': [double] filter order (filter length - 1).
%   High-pass:  Hamming windowed sinc FIR filter, passband edge 0.1 Hz, filter order 16500 (estimated), transition bandwidth 0.1 Hz, cutoff frequency (-6 dB) 0.05 Hz
sets.filters(1).locutoff=0.1;        
sets.filters(1).hicutoff=[];
sets.filters(1).order=16500;
%   Low-pass:   Hamming windowed sinc FIR filter, passband edge 80 Hz, filter order 84 (estimated), transition bandwidth 20 Hz, cutoff frequency (-6 dB) 90 Hz
sets.filters(2).locutoff=[];
sets.filters(2).hicutoff=80;
sets.filters(2).order=84;

%16-'downsample_freq': [double] number specifying downsample frequency.
sets.downsample_freq=500;


% --- Epoching ---

%17-'epoch_length': [double] array with start and end point of epoch in seconds. Baseline is defined as negative number
%and 0 correspond to the locked trigger.
sets.epoch_length=[-0.950, 2.050]; % extended beyond 2 to include markers S200 (incorrect) after deadline (2sec) = missing response.

%18-'baseline': [double] start and end point of baseline to remove, in milliseconds.
sets.baseline=[-500 -300];          % note: baseline used for preprocessing


%Save a copy of 'sets' in 'settings_dir':
if ispc
    save([sets.settings_dir 'preprocessing_settings.mat'],'sets');
else
   save([sets.settings_dir 'preprocessing_settings_mac.mat'],'sets');
end



%% 2. Import files, filtering and cleaning operations

%This section loops through subjects, loading the corresponding .vhdr file from 'sets.raw_dir'.
%It then performs the following preprocessing steps:
%   -Removes unused channels;
%   -Loads a EEGlab .ced file into EEG.chanlocs;
%   -Applies filters;
%   -Applies higher HP filter for ICA training data + epochs and baseline;
%   -Run ICA on training data;
%   -Transfer ICA weights to original EEG;

%Data are saved individually as .set and .fdt files in the folder specified in 'sets.continuous_data_dir'.
%Saved data and corresponding EEG struct's field 'set name' are defined as 'SUB_ID_continuous_filtered'.

% Loop over participant datasets.
for si=47:length(sets.sub_idx) 
    
    sprintf('Working on subject %d...',si)
    

    %Import current subject EEG stored in 'sets.raw_dir':
    EEG=pop_loadbv(sets.raw_dir, ...
                   ['dimensions_EEG_s0' ...
                    sets.sub_idx{si}...
                    '.vhdr']);            
    
    
    %Apply filters by looping through entries in 'sets.filter' struct and applying corresponding
    %parameters:
    for f=1:length(sets.filters)
        
        EEG=pop_eegfiltnew(EEG,...
            sets.filters(f).locutoff,...
            sets.filters(f).hicutoff,...
            sets.filters(f).order,...
            0,...
            0,...
            0);
    end
    
    % Create a copy of EEG to be used as ICA training data, and adapt:      (cf. Dimigen, O. (2020). https://doi.org/10.1016/j.neuroimage.2019.116117) 
    % 1) Apply additional HP filter: Hamming windowed sinc FIR filter, passband edge 2.0 Hz, filter order 826 (estimated), transition bandwidth 2.0 Hz, cutoff frequency (-6 dB) 1.0 Hz
    EEG_train = pop_eegfiltnew(EEG,...
            2.0,...                     % lower passband edge
            [],...
            826,...                     % filter order (estimated)
            0,...
            0,...
            0);
    % 2) Cut into epochs (short epochs after serialobj markers)
    EEG_train = pop_epoch(EEG_train,...
        [sets.onset_dim.markers{:}],... 
        [-0.500 1.5]);
    % 3) Baseline correction. [] -> subtract entire epoch mean
    EEG_train = pop_rmbase(EEG_train, [], []); 
   
    
    % Run ICA
    EEG_train = pop_runica(EEG_train,'extended',1,'interupt','on');%,'icatype','binica'
                
    % Remove any existing ICA solutions from original EEG dataset
    EEG.icaact      = [];
    EEG.icasphere   = [];
    EEG.icaweights  = [];
    EEG.icachansind = [];
    EEG.icawinv     = [];

    % Pass ICA info (unmixing matrix) from training EEG data to original EEG
    EEG.icaweights  = EEG_train.icaweights;
    EEG.icasphere   = EEG_train.icasphere;
    EEG.icachansind = EEG_train.icachansind;
    EEG = eeg_checkset(EEG); % let EEGLAB re-compute EEG.icaact & EEG.icawinv
    
    
    %Update current subject dataset name:
    EEG.setname=[sets.sub_idx{si} '_continuous_filtered'];
    %Update current subject .set file name:
    EEG.filename=[sets.sub_idx{si} '_continuous_filtered.set'];
    %Update .set file path:
    EEG.filepath=sets.continuous_data_dir;
    %Save current dataset as .set file in 'sets.continuous_data_dir' folder:
    EEG=pop_saveset(EEG,...
                    'filename', EEG.filename,...
                    'filepath', EEG.filepath);
        
end
