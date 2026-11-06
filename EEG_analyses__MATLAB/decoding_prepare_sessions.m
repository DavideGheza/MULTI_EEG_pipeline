%% Prepare V3 cleaned EEG data for session-level decoding
% Retains correct, unrejected V3 trials; extracts task and stimulus markers;
% sorts EEG trials by task; and applies the identical selection and ordering
% to the behavioral variables used by decoding and RSA analyses.
%
% Pipeline summary:
% 1. Load each participant's cleaned EEGLAB dataset.
% 2. Remove incorrect epochs identified by S200. For stimulus-locked data,
%    the selection window extends beyond the response deadline so it also
%    includes an S200 sent after a missing response. Missing responses have
%    no response-locked epoch and are removed from the behavioral data.
% 3. Recover task and stimulus identities from EEG events: directly from
%    epoch events when stimulus-locked, or from the preceding urevents when
%    response-locked. These EEG-derived labels provide an independent check
%    that EEG trials and behavioral trials remain aligned.
% 4. Sort EEG trials by task identity and apply the same order to behavior.
%    Later decoding scripts sort again using the eight cue-by-side classes,
%    so this first sort is largely redundant but harmless and retained for
%    compatibility with the established workflow.
% 5. Optionally average every five time samples for faster exploratory runs.
% 6. Load raw behavior and the preprocessing archive that flags rejected EEG
%    trials. Noisy EEG trials were rejected in pp03; the rejection archive is
%    therefore required to apply the same exclusions to the behavioral data.
% 7. Compute sequential variables before removing trials, so previous- and
%    next-trial variables refer to the original experimental sequence.
% 8. Apply EEG rejection, accuracy filtering, and EEG trial ordering to the
%    behavioral variables; verify task-label alignment; and save the aligned,
%    augmented behavioral dataset (V3_behav_seqinfo_trimmed_[response_]sbj_N.mat)
%    and the sorted EEG data (sessions_data_[response_]sbj_N.mat).
% Author: Davide Gheza

% Select stimulus- or response-locked inputs.
onset = false;
if onset
    suffix = 'onset';
    saveout_suffix = '';
else
    suffix = 'response';
    saveout_suffix = 'response_';
end

downsample = false;

%% 1. Load settings

if ispc
    settings_file = 'C:\Users\gdavide\OneDrive\Documenti\MATLAB\dg_workspace22\dimensions\analyses\scripts_and_settings\preprocessing_settings.mat';
    settings_name = 'preprocessing_settings.mat';
else
    settings_file = '/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/scripts_and_settings/preprocessing_settings_mac.mat';
    settings_name = 'preprocessing_settings_mac.mat';
end
load(settings_file);

sets.sessions_data_dir = [sets.main_saving_dir 'sessions_data/'];
if ~exist(sets.sessions_data_dir,'dir')
    mkdir(sets.sessions_data_dir);
end
save([sets.settings_dir settings_name],'sets');

%% 2. Prepare each V3 participant

for si = 22:length(sets.sub_idx)
    eeglab
    data = [];

    EEG = pop_loadset( ...
        'filename',[sets.sub_idx{si} suffix '_trialrej.set'], ...
        'filepath',sets.clean_data_dir);

    % V3 contains one test session. Remove epochs containing an incorrect-
    % response marker in the appropriate post-locking latency window.
    if onset
        EEG_test = pop_selectevent(EEG,'latency','0 <= 2050','type',{'S200'}, ...
            'deleteevents','off','deleteepochs','on','invertepochs','on');
    else
        EEG_test = pop_selectevent(EEG,'latency','0 <= 200','type',{'S200'}, ...
            'deleteevents','off','deleteepochs','on','invertepochs','on');
    end

    % Participant 23's stimulus-locked recording ends with an incomplete
    % trial after a MATLAB crash; S199 identifies that final trial.
    if si == 23 && onset == 1
        EEG_test = pop_selectevent(EEG_test,'type',{'S199'}, ...
            'deleteevents','off','deleteepochs','on','invertepochs','on');
    end

    data.ses1 = permute(EEG_test.data,[3,1,2]); % trials x channels x time

    %% 3. Recover task and stimulus markers

    n_epochs = numel(EEG_test.epoch);
    labels = nan(n_epochs,1);
    onsettriggers = nan(n_epochs,1);

    for i = 1:n_epochs
        if ischar(EEG_test.epoch(1,i).eventtype{1})
            if onset
                % For V3 stimulus-locked epochs, the first two events are
                % the task cue and stimulus identity, respectively.
                labels(i) = sscanf(EEG_test.epoch(1,i).eventtype{1}(2:end),'%d');
                onsettriggers(i) = sscanf(EEG_test.epoch(1,i).eventtype{2}(2:end),'%d');
            else
                % For response-locked epochs, recover the two immediately
                % preceding events from the continuous urevent sequence.
                response_event = [EEG_test.epoch(1,i).eventlatency{:}] == 0;
                response_urevent = EEG_test.epoch(1,i).eventurevent{response_event};
                labels(i) = sscanf(EEG_test.urevent(response_urevent-2).type(2:end),'%d');
                onsettriggers(i) = sscanf(EEG_test.urevent(response_urevent-1).type(2:end),'%d');
            end
        end
        labels(i) = mod(labels(i),round(labels(i)/10)); % final digit is task ID
    end

    [labels,I] = sort(labels,'ascend');
    data.ses1 = data.ses1(I,:,:);
    data.labels_ses1 = labels;
    data.onsettriggers_ses1 = onsettriggers(I);

    if downsample
        data_size = size(data.ses1);
        data.ses1 = squeeze(mean(reshape(data.ses1, ...
            data_size(1),data_size(2),5,[]),3));
    end

    %% 4. Match behavioral data to the retained EEG trials

    behavioral_files = dir([sets.raw_dir '/*PC*']);
    if onset
        rejection_files = dir([sets.clean_data_dir '/*log_preproc.mat*']);
    else
        rejection_files = dir([sets.clean_data_dir '/*log_preproc_' suffix '.mat*']);
    end
    load([sets.clean_data_dir rejection_files.name]);

    behav_file = load([sets.raw_dir '/' behavioral_files(si).name]);
    behav = behav_file.data;
    rejarchive = logfile(si).rej_trials;
    if si == 23 && onset == 1
        rejarchive(end) = [];
    end

    % Early V3 files used different names for edge and pattern features.
    if si == 22 || si == 23
        behav.edges = behav.borders;
    end
    if si == 22
        behav.patterns = behav.orientations;
    end

    blkwise_PCs = reshape([behav.pcs{ceil(behav.block/8)}],4,[])';

    % Sequential variables are computed before any trial removal so their
    % values continue to describe neighbors in the original task sequence.
    behav.n(:,1) = 1:length(behav.rt);
    behav.newblock = [1;diff(behav.block)];
    behav.endblock = [diff(behav.block);1];
    behav.RTtp1 = [behav.rt(2:end);nan];
    behav.sum_con = sum(behav.congruency==1,2);
    behav.prev_sum_con = [nan;behav.sum_con(1:end-1)];
    behav.cse = behav.sum_con + behav.prev_sum_con*10;
    behav.prev_congruency = [nan(1,4);behav.congruency(1:end-1,:)];
    behav.post_correct = [1;behav.accuracy(1:end-1)];
    behav.fut_congruency = [behav.congruency(2:end,:);nan(1,4)];
    behav.fut_correct = [behav.accuracy(2:end);nan];
    behav.fut_taskcount = [behav.taskcount(2:end);nan];

    vector_fields = {'accuracy','rt','response','taskcount','pc','task','block', ...
        'n','newblock','endblock','RTtp1','sum_con','prev_sum_con','cse', ...
        'post_correct','fut_correct','fut_taskcount'};
    matrix_fields = {'congruency','colors','shapes','edges','patterns', ...
        'prev_congruency','fut_congruency'};

    % Participant 32 lacks the first four EEG trials.
    if si == 32
        blkwise_PCs(1:4,:) = [];
        for field = vector_fields
            behav.(field{1})(1:4) = [];
        end
        for field = matrix_fields
            behav.(field{1})(1:4,:) = [];
        end
    end

    % Missing-response trials were never response-locked and must therefore
    % be removed before applying the preprocessing rejection archive.
    if onset == 0
        missing_response = behav.rt == -1;
        blkwise_PCs(missing_response,:) = [];
        for field = vector_fields
            behav.(field{1})(missing_response) = [];
        end
        for field = matrix_fields
            behav.(field{1})(missing_response,:) = [];
        end
    end

    if size(blkwise_PCs,1) ~= size(rejarchive,2)
        error("mismatch behavioral data - trial rejection list")
    end

    % Apply preprocessing rejection, behavioral accuracy, and the EEG task
    % sort in that order. This is the same ordering used for data.ses1.
    retained = ~rejarchive;
    correct = logical(behav.accuracy(retained));
    blkwise_PCs = blkwise_PCs(retained,:);
    blkwise_PCs = blkwise_PCs(correct,:);
    behav_trim.blkwise_PCs = blkwise_PCs(I,:);

    output_fields = {'rt','response','taskcount','congruency','pc','task', ...
        'block','colors','shapes','edges','patterns','n','newblock', ...
        'endblock','RTtp1','sum_con','prev_sum_con','cse', ...
        'prev_congruency','post_correct','fut_congruency','fut_correct', ...
        'fut_taskcount'};
    for field = output_fields
        if ismember(field,matrix_fields)
            values = behav.(field{1})(retained,:);
            values = values(correct,:);
            behav_trim.(field{1}) = values(I,:);
        else
            values = behav.(field{1})(retained);
            values = values(correct);
            behav_trim.(field{1}) = values(I);
        end
    end

    if any(diff(data.labels_ses1) ~= diff(behav_trim.task))
        error("mismatch behavioral data - EEG data")
    end

    cd(sets.sessions_data_dir)
    fname = sprintf(['V3_behav_seqinfo_trimmed_' saveout_suffix ...
        'sbj_%d.mat'],si);
    save(fname,'behav_trim')

    % Save the task-sorted EEG data, labels, and stimulus triggers; input to
    % decoding_tf_on_sessions*.m.
    fname = sprintf(['sessions_data_' saveout_suffix 'sbj_%d.mat'],si);
    save(fname,'data')
end
