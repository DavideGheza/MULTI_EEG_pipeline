%% Reject artifactual trials and apply average rereferencing
% Loads the epoched data from preprocessing step 2, applies the configured
% trial-rejection criteria, rereferences retained epochs, and saves cleaned
% participant datasets and a preprocessing log.
% Author: Davide Gheza



%% 1. Load settings

%Struct 'sets' is loaded from the 'scripts_and_settings' folder:
if ispc
    load('C:\Users\gdavide\OneDrive\Documenti\MATLAB\dg_workspace22\dimensions\analyses\scripts_and_settings\preprocessing_settings.mat');
else
    load('/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/scripts_and_settings/preprocessing_settings_mac.mat');
end

%Add path to sub-folder of 'main_dir' to store datasets after trial rejection (clean data):
sets.clean_data_dir=[sets.main_saving_dir 'clean_data/'];
%If the folder doesn't exist yet, create one:
if ~exist(sets.clean_data_dir,'dir')
    mkdir(sets.clean_data_dir);
end

%Save updated 'sets' struct:
if ispc
    save([sets.settings_dir 'preprocessing_settings.mat'],'sets');
else
   save([sets.settings_dir 'preprocessing_settings_mac.mat'],'sets');
end

%Select onset or response events
onset = false;
if onset
    events = [sets.taskcue.markers]; % count trials based on cued dimension    
    suffix = 'onset';
else
    events = [sets.classresp_pc.markers, sets.resp_pc.markers];
    suffix = 'response';
end

%% 2. Trial rejection

% Loop over participant datasets.
for si=22:length(sets.sub_idx)
     
    eeglab %open it
    
    %Load current subject EEG struct stored in 'sets.epoch_data_dir':
    EEG=pop_loadset('filename',[sets.sub_idx{si}...
                suffix '_epoched_ICApruned.set'],...
                'filepath',[sets.epoch_data_dir]...
                );
                
    %Reject ICs marked for pruning:
    EEG=pop_subcomp( EEG, [], 0);
                
    %Remove baseline, again (to be sure ICA pruning did not affect baseline):
    EEG=pop_rmbase(EEG,sets.baseline);
    
    %Count the number of epochs per condition, useful for log file (see below)(adapted from Schettino 2016 NI)
    trigtotepochs=[];                                                                   % reset vector
    for irej=1:numel(EEG.epoch)                                                         % loop through epochs
            if ischar(EEG.epoch(1,irej).eventtype{1})                                       %convert eventtype from str to double (ignore first character = S)
                EEG.epoch(1,irej).eventtype{1} = sscanf(EEG.epoch(1,irej).eventtype{1}(2:end), '%d');
%                 EEG.epoch(1,irej).eventtype{2} = sscanf(EEG.epoch(1,irej).eventtype{2}(2:end), '%d');
            end
        trigtotepochs(irej,:)=EEG.epoch(1,irej).eventtype{1};                           % create vector with the first trigger (serialobj) of each epoch (total epochs before trial rejection)
    end
    for k=1:numel(events)                                                     % loop through conditions
        totepochs(:,k)=sum(trigtotepochs'==sscanf(events{k}(2:end), '%d'));     % find all triggers belonging to each condition and sum them
    end
    
    %Clean data 2/2: trial rejection
    
    % Interpolate bad chans first (first spot datasets with bad channels...)
    if si==999
    EEG = pop_interp(EEG, [6], 'spherical'); % interpolate (e.g. F5 (chan 6) for si==999).     
    end

    rejarchive=[]; %delete marked trials archive
    %Automatic artifact rejection
    EEG = pop_eegthresh(EEG,1,[1:5 7:29 32:63] ,-120,120,-0.500,1.500,1,0); % extreme values +- 120uv         
%     EEG = pop_rejtrend(EEG,1,[1:5 7:29 32:63] ,500,120,0.3,1,0,0);         % abnodmal trends (500 = winsize in datapoints (~half 1155 of epoch), 90 = max slope in uV/epoch, 0.3 = R^2 limit)
    EEG = pop_jointprob(EEG,1,[1:5 7:29 32:63] ,6,3,0,0,0,[],0);          % improbable data. SD single ch 6, all ch 3
    %Manual artifact rejection 
    eeglab redraw               %(debug until next line and visual inspect data + select bad trials before continuing)

    %Merge automatic and manual selections
    EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);

    % print total marked trials
    sum(EEG.reject.rejglobal)       
    % store marked trials (logical) indexes
    rejarchive = EEG.reject.rejglobal;
    % reject marked trials
    EEG = pop_rejepoch( EEG, EEG.reject.rejglobal ,0);
    
    % calculate the percentage of rejected epochs per condition, after artifact rejection (adapted from Schettino 2016 NI)
    trigepochsleft=[];                                                                     % reset vector
    for irej2=1:numel(EEG.epoch)
            if ischar(EEG.epoch(1,irej2).eventtype{1})                                       % convert eventtype from str to double
                EEG.epoch(1,irej2).eventtype{1} = sscanf(EEG.epoch(1,irej2).eventtype{1}(2:end), '%d');
%                 EEG.epoch(1,irej2).eventtype{2} = sscanf(EEG.epoch(1,irej2).eventtype{2}(2:end), '%d');
            end
        trigepochsleft(irej2,:)=EEG.epoch(1,irej2).eventtype{1};                           % create vector with the first trigger (serialobj) of each epoch left (after artifact rejection)
    end
    for j=1:numel(events)                                                        % loop through conditions
        epochsleft(:,j)=sum(trigepochsleft'==sscanf(events{j}(2:end), '%d'));      % find all triggers belonging to each condition and sum them
        rejepxcond(:,j)=100-(epochsleft(:,j)*100)/totepochs(:,j);                          % calculate percentage of rejected epochs separately for each condition
    end
    
    % save log file in a structure variable, to report the percentage of rejected epochs
    logfile(:,si)=struct('subject',cell2mat(sets.sub_idx(si)),... % participant number
                         'rej_trials',rejarchive,...              % marked & rejected trials indexes
                         'perc_rejepochs',rejepxcond);            % percentage of rejected epochs per condition
                     
    % edit chan location, add Cz
    EEG=pop_chanedit(EEG, 'append',63,'changefield',{64,'labels','Cz'},'changefield',{64,'theta','0'},'changefield',{64,'radius','0'},'convert',{'topo2all'});
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
                     
    %Re-reference the data to all 63 scalp-electrodes (average-reference, ocular electrodes 60 30 31 are excluded): 
    EEG = pop_reref( EEG, [],'refloc',struct('labels',{'Cz'},'sph_radius',{[]},'sph_theta',{0},'sph_phi',{90},'theta',{0},'radius',{0},'X',{5.8317e-17},'Y',{0},'Z',{0.95238},'type',{''},'ref',{''},'urchan',{[]},'datachan',{0}),'exclude',[6 30 31] );

    %Discard EOG channels
    EEG = pop_select( EEG, 'nochannel',{'HEOGL','VEOG','HEOGR'});


    % Deal with events that must be from the following trial. Edit EEG.epoch: look for taskcues, rename (to negative) events from any second taskcue onward 
    % assumes that exogenous markers can leak from the following trial, but not from the previous one (given epoching -1 2, and ITI => 1.3)
    for d=1:numel(EEG.epoch)                                                         % loop through epochs
            if sum( ismember(EEG.epoch(1,d).eventtype, events) ) > 1                 % find epochs containing more than one taskcue maker (first marker in any trial)
               extra_trial_taskcue = find(ismember(EEG.epoch(1,d).eventtype, events));
               extra_trial_taskcue = extra_trial_taskcue(2);
               for e = extra_trial_taskcue:length(EEG.epoch(1,d).eventtype)          % loop through the following events (markers from the following trial)
                   EEG.epoch(1,d).eventtype{e} = num2str( sscanf(EEG.epoch(1,d).eventtype{e}(2:end), '%d') *-1 ); % rename those with their marker number (without S) * -1
                   EEG.event(EEG.epoch(1,d).event(e)).type = EEG.epoch(1,d).eventtype{e};                   % rename event also in EEG.event.type field
               end
            end
    end
   
    %Update current subject dataset name:
    EEG.setname=[sets.sub_idx{si}...
                     suffix '_trialrej'];

    %Update current subject .set file name:
    EEG.filename=[sets.sub_idx{si}...
                     suffix '_trialrej.set'];  

    %Update .set file path:
    EEG.filepath=sets.clean_data_dir;

    %Save current dataset as .set file in 'sets.epoch_data_dir' subfolder:
    EEG=pop_saveset(EEG,...
                    'filename', EEG.filename,...
                    'filepath', EEG.filepath);
                
end

% save log file
cd(sets.clean_data_dir)
log_filename = ['log_preproc_' suffix '.mat'];
save(log_filename, 'logfile')

% load log_preproc.mat % load log file

% create a .csv file with percentage of rejected epochs per condition
% (to extract the rejected trials indexes, refer to log_preproc.mat)
ssj=[]; temp=[]; temp_abs=[];
for itemp=22:numel(logfile)                                  % loop through participants
    ssj(itemp,1)=itemp;                                     % participant number
    temp(itemp,:)=logfile(1,itemp).perc_rejepochs;          % percentage of rejected epochs per condition
    temp_abs(itemp,:)=sum(logfile(1,itemp).rej_trials);     % absolute number of rejected trials
end
if onset == 1
    head_summary=[{'ssj'}, sets.taskcue.names, {'abs_total_rej'}]; % file header
else
    head_summary=[{'ssj'}, sets.classresp_pc.names, sets.resp_pc.names, {'abs_total_rej'}]; % file header
end

summary=[ssj temp temp_abs]; % merge data
summary=[head_summary;cellfun(@num2str,num2cell(summary),'UniformOutput',false)]; % merge header and data
clear ssj temp_abs temp % clear temporary variables
cell2csv([sets.clean_data_dir 'logfile_trialrej' suffix '.csv'],summary); % save as .csv
