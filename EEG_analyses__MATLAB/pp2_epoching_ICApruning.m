%% Epoch data, remove ICA components, and baseline-correct
% Loads the continuous filtered datasets from preprocessing step 1, creates
% stimulus- or response-locked epochs, removes selected ICA components, and
% saves the epoched datasets.
% Author: Davide Gheza

%% 1. Load settings

%Struct 'sets' is loaded from the 'scripts_and_settings' folder:
if ispc
    load('C:\Users\gdavide\OneDrive\Documenti\MATLAB\dg_workspace22\dimensions\analyses\scripts_and_settings\preprocessing_settings.mat');
else
    load('/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/scripts_and_settings/preprocessing_settings_mac.mat');
end
%% 2. Epoching and baseline correction

% Extract epochs around the selected event markers and apply the configured
% baseline. Each participant is saved as EEGLAB `.set`/`.fdt` data in
% `sets.epoch_data_dir`.

%Select onset or response events
onset = false;
if onset
    events = [sets.onset_dim.markers{:}];
    suffix = 'onset';
else
    events = [sets.classresp_pc.markers, sets.resp_pc.markers]; % NB: no epoches for missing response trials
    suffix = 'response';
end

% Loop over participant datasets.
for si=22:48%[14:21]%:length(sets.sub_idx)
     
    eeglab %open it
    
    %Load current subject EEG struct stored in 'sets.continuous_data_dir':
    EEG=pop_loadset('filename',[sets.sub_idx{si}...
                    '_continuous_filtered.set'],...
                    'filepath',[sets.continuous_data_dir]...
                    );
                
    %Epoch extraction, using time window defined in 'settings.epoch_length'. 
    %Epoched data are time-locked to cue presentation markers, defined in 'sets.xxx.markers':
    EEG=pop_epoch(EEG,...
                         events,...
                         sets.epoch_length,...
                         'epochinfo', 'yes');
                     
    %Remove baseline:
    EEG=pop_rmbase(EEG,sets.baseline);
    
    
                     

%% 3. Clean data 1/2: reject ICA components
 
    %Check if a .csv file with component indexes exists. If it exists, load it. If it doesn't, preallocate a vector to 
    %store flagged components. 
    %Initial values are all set to 0 ('unflagged').
    if exist([sets.epoch_data_dir 'ICA_comp_archive_' sets.sub_idx{si} '.csv'],'file')
        comparchive=dlmread([sets.epoch_data_dir 'ICA_comp_archive_' sets.sub_idx{si} '.csv']);
        EEG.reject.gcompreject = comparchive;
    else
        comparchive=zeros(1,length(EEG.icachansind)); 
    

                
    % ICA component selection for pruning:
    % 1 - IClabel: automatic IC classifier
    EEG = iclabel(EEG, 'default');
    %EEG = eeg_checkset( EEG );
    % 2 - Manual ICA component selection
    pop_eegplot( EEG, 1, 1, 1);
    pop_eegplot( EEG, 0, 1, 1);
    pop_selectcomps(EEG, 1:32);            
    eeglab redraw;                  % (debug until next line and select components here before continuing)
    
    % Saving selected components indexes
    comparchive (1, 1:length(EEG.reject.gcompreject)) = EEG.reject.gcompreject;     %%     store comp removed 
    %save rej components index
    csvwrite([sets.epoch_data_dir,'ICA_comp_archive_',sets.sub_idx{si},'.csv'],comparchive);
    comparchive=[];
    end

    %Update current subject dataset name:
    EEG.setname=[sets.sub_idx{si}...
                      suffix '_epoched_ICApruned'];

    %Update current subject .set file name:
    EEG.filename=[sets.sub_idx{si}...
                      suffix '_epoched_ICApruned.set'];  
                  
    %Update .set file path:
    EEG.filepath=sets.epoch_data_dir;

    %Save current dataset as .set file in 'sets.epoch_data_dir' subfolder:
    EEG=pop_saveset(EEG,...
                    'filename', EEG.filename,...
                    'filepath', EEG.filepath);
                
end
