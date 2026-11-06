%% Repeated-CV multiclass decoding for single-trial RSA
% Trains an eight-class probabilistic discriminant classifier for the
% cue-by-side conjunction. Trials are balanced across proportion-congruency conditions and
% classes, whitened using training data, and evaluated with repeated six-fold
% cross-validation. Shrinkage is optimized on balanced training data only.
%
% The output retains each trial in its original order and stores class
% probabilities separately for the four proportion-congruency sessions. These
% probability profiles are the dependent data for the model-1 and model-3 RSA
% fits and provide the decoding-accuracy panel in Figure 4.
% Author: Davide Gheza

function [result] = RSA_distances_classprob_kfold_SIDEnCue_multiclass_repCV_HPC(si)

local=0;
% addAttachedFiles(Pool,'/mnt/beegfs/scratch/gheza/')
%% 1. Load settings

% Select stimulus- or response-locked inputs.
onset = true;
if onset
    suffix = '';
else
    suffix = 'response_';   % response-locked requires precomputing stimulus-locked to get a fixed pre-stimulus baseline
end


%%
sn = [22,24:48];

if local==0
    d = load(['/mnt/beegfs/scratch/gheza/EEG_data/sessions_data_' suffix 'tf_decoding_BN_sbj_'...% for tf!       If time-domain, delete 'tf_decoding_BN_'
            num2str(sn(si))...
            '.mat']);
    behav = load(['/mnt/beegfs/scratch/gheza/EEG_data/behav/V3_behav_seqinfo_trimmed_' suffix 'sbj_' num2str(sn(si))]);    % load behavioral data
else
    root = '/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/preprocessed_data/sessions_data_tf/'; % for local
    d=load(fullfile(root, ['sessions_data_' suffix 'tf_decoding_BN_sbj_' num2str(sn(si)) '.mat']));                                                             % for local
    behav = load(['/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/preprocessed_data/sessions_data/V3_behav_seqinfo_trimmed_' suffix 'sbj_' num2str(sn(si))]); % for local
end
data = d.data;
behav = behav.behav_trim;


% define vectors for indexing trials by PC
Neu_1_2 = nan(length(behav.blkwise_PCs),2);     %preallocate - two neutral dimensions per trial
MI = nan(length(behav.blkwise_PCs),1);   
MC = nan(length(behav.blkwise_PCs),1);
% find task based on PC
for t = 1:length(behav.blkwise_PCs)      
    Neu_1_2(t,:) = find(behav.blkwise_PCs(t,:) == 1/2);
    MI(t,:) = find(behav.blkwise_PCs(t,:) == 1/3);
    MC(t,:) = find(behav.blkwise_PCs(t,:) == 2/3);
    if Neu_1_2(t,1) == 1 && Neu_1_2(t,2) == 4               % Neu_1 and Neu_2 are artificially segregated as neutral dimension on left vs right column (of behav.blkwise_PCs). Hence, if Neu1 happens to be contiguous to Neu2, Neu1 never contains cond 4, and Neu2 never contains cond 1.
        Neu_1_2(t,1) = 4;                                   % To restore original balance, swap neu left and neu right, when they are respectively dimension 1 or 4 (i.e., follow right-shift logic used for associating taskID-PCs over blocks)
        Neu_1_2(t,2) = 1;
    elseif Neu_1_2(t,1) == 1 && Neu_1_2(t,2) == 3 || Neu_1_2(t,1) == 2 && Neu_1_2(t,2) == 4 % Alternatively, if Neu1 and Neu2 are not contiguous, Neu1 will always contain cond 1 and 2, and Neu2 3 and 4. 
        if Neu_1_2(t,1) == 2 && MI(t) == 1
            Neu_1_2(t,1) = 4;                                   % To restore original balance, swap neu left and neu right, based on relative position of N1 and MI (i.e., follow right-shift logic used for associating taskID-PCs over blocks)
            Neu_1_2(t,2) = 2;
        elseif Neu_1_2(t,1) == 1 && MI(t) == 4
            Neu_1_2(t,1) = 3;                                   % To restore original balance, swap neu left and neu right, based on relative position of N1 and MI  (i.e., follow right-shift logic used for associating taskID-PCs over blocks)
            Neu_1_2(t,2) = 1;
        end
    end
end
Neu_1 = Neu_1_2(:,1);
Neu_2 = Neu_1_2(:,2);
PC_matrix = [MI,MC,Neu_1,Neu_2];

% convert feature levels to standard labels [1 2]
behav.colors = behav.colors;
behav.shapes = behav.shapes+1;
behav.edges = behav.edges+1;
behav.patterns = ((behav.patterns-90)==90)+1;

% load original behavioral data - to get correct (target) features          (taken from 8feat script as of 11/12/2022)
if local==0
    full_behav = load(['/mnt/beegfs/scratch/gheza/EEG_data/behav/dimensions_PC_s' num2str(sn(si)) '_0' num2str(sn(si))]);
else
    sbj = dir('/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/data/EEG_data/*PC*'); % for local
    full_behav = load(['/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/data/EEG_data', '/', sbj(sn(si)).name]);
end
correct_feature = full_behav.data.correct_feature;
% Align feature-level coding so that levels 1 and 2 consistently represent
% the same response sides across color, shape, edge, and pattern. This makes
% the side component of the eight cue-by-side classes comparable across
% dimensions and participants; it does not remove side-related information.
correct_feature(1)=correct_feature(1)-1; % make color coding 0-1
correct_feature(4)=correct_feature(4)-90==90; % make pattern coding 0-1
flip_shape = correct_feature(1)~=correct_feature(2);
flip_edge = correct_feature(1)~=correct_feature(3);
flip_pattern = correct_feature(1)~=correct_feature(4);
behav.shapes = behav.shapes*(-1*flip_shape+(1-flip_shape))+(3*flip_shape); % swap 1 and 2 if flip_shape
behav.edges = behav.edges*(-1*flip_edge+(1-flip_edge))+(3*flip_edge); % swap 1 and 2 if flip_edge
behav.patterns = behav.patterns*(-1*flip_pattern+(1-flip_pattern))+(3*flip_pattern); % swap 1 and 2 if flip_pattern

order = (1:length(data.labels_ses1))'; % index of trials as from data (matching behav)
% Session 1 = PC neu 1
sessions(1).data = data.ses1(data.labels_ses1==Neu_1, :,:);
sessions(1).labels = data.labels_ses1(data.labels_ses1==Neu_1);
sessions(1).colors = behav.colors(data.labels_ses1==Neu_1);
sessions(1).shapes = behav.shapes(data.labels_ses1==Neu_1);
sessions(1).edges = behav.edges(data.labels_ses1==Neu_1);
sessions(1).patterns = behav.patterns(data.labels_ses1==Neu_1);
sessions(1).ord = order(data.labels_ses1==Neu_1);
% Session 4 = PC neu 2
sessions(4).data = data.ses1(data.labels_ses1==Neu_2, :,:);
sessions(4).labels = data.labels_ses1(data.labels_ses1==Neu_2);
sessions(4).colors = behav.colors(data.labels_ses1==Neu_2);
sessions(4).shapes = behav.shapes(data.labels_ses1==Neu_2);
sessions(4).edges = behav.edges(data.labels_ses1==Neu_2);
sessions(4).patterns = behav.patterns(data.labels_ses1==Neu_2);
sessions(4).ord = order(data.labels_ses1==Neu_2);
% Session 2 = PC MI
sessions(2).data = data.ses1(data.labels_ses1==MI, :,:);
sessions(2).labels = data.labels_ses1(data.labels_ses1==MI);
sessions(2).colors = behav.colors(data.labels_ses1==MI);
sessions(2).shapes = behav.shapes(data.labels_ses1==MI);
sessions(2).edges = behav.edges(data.labels_ses1==MI);
sessions(2).patterns = behav.patterns(data.labels_ses1==MI);
sessions(2).ord = order(data.labels_ses1==MI);
% Session 3 = PC MC
sessions(3).data = data.ses1(data.labels_ses1==MC, :,:);
sessions(3).labels = data.labels_ses1(data.labels_ses1==MC);
sessions(3).colors = behav.colors(data.labels_ses1==MC);
sessions(3).shapes = behav.shapes(data.labels_ses1==MC);
sessions(3).edges = behav.edges(data.labels_ses1==MC);
sessions(3).patterns = behav.patterns(data.labels_ses1==MC);
sessions(3).ord = order(data.labels_ses1==MC);

% Parameters
n_perm = 6;  % number of cross-validation folds
n_conditions = length(unique(data.labels_ses1))*2; % cue * side = 8 conditions
n_sensors = size(data.ses1, 2);
n_time = size(data.ses1, 3);
n_sessions = length(sessions);
n_meta = 8; % number of repetitions of k-fold after random shuffling
gamma = nan(n_time); % store LDA gamma optimized for n times

% Cross-validated decoding steps:
% 1 Randomize trial order independently within each condition.
% 2 Whiten the training and test data using the training-data covariance.
% 3 Fit the LDA classifier to the training data.
% 4 Compute classification accuracy and class probabilities on test data.


    % preallocate LDA results
    result.lda_cond_acc = nan(n_sessions, n_meta, n_perm, n_conditions, n_time);
    result.lda = nan(n_sessions, n_meta, n_perm, n_time);
    result.lda_prob_full = nan(n_sessions, n_meta, n_perm, n_conditions, n_conditions, n_time);

    for s = 1:n_sessions

        X = [];
        % X1 = cued dimension is COLOR      
        X(1).data = sessions(s).data(sessions(s).labels==1, :,:);
        X(1).labels = sessions(s).labels(sessions(s).labels==1);
        X(1).colors = sessions(s).colors(sessions(s).labels==1);                % keep only left feature (right is complementary)
        X(1).shapes = sessions(s).shapes(sessions(s).labels==1);                % keep only left feature (right is complementary)
        X(1).edges = sessions(s).edges(sessions(s).labels==1);                  % keep only left feature (right is complementary)
        X(1).patterns = sessions(s).patterns(sessions(s).labels==1);            % keep only left feature (right is complementary)
        X(1).labels = X(1).labels+((X(1).colors-1)*4); 
        X(1).ord = sessions(s).ord(sessions(s).labels==1);   
        % X2 = cued dimension is SHAPE      
        X(2).data = sessions(s).data(sessions(s).labels==2, :,:);
        X(2).labels = sessions(s).labels(sessions(s).labels==2);
        X(2).colors = sessions(s).colors(sessions(s).labels==2);                % keep only left feature (right is complementary)
        X(2).shapes = sessions(s).shapes(sessions(s).labels==2);                % keep only left feature (right is complementary)
        X(2).edges = sessions(s).edges(sessions(s).labels==2);                  % keep only left feature (right is complementary)
        X(2).patterns = sessions(s).patterns(sessions(s).labels==2);            % keep only left feature (right is complementary)
        X(2).labels = X(2).labels+((X(2).shapes-1)*4);        
        X(2).ord = sessions(s).ord(sessions(s).labels==2);   
        % X3 = cued dimension is EDGE      
        X(3).data = sessions(s).data(sessions(s).labels==3, :,:);
        X(3).labels = sessions(s).labels(sessions(s).labels==3);
        X(3).colors = sessions(s).colors(sessions(s).labels==3);                % keep only left feature (right is complementary)
        X(3).shapes = sessions(s).shapes(sessions(s).labels==3);                % keep only left feature (right is complementary)
        X(3).edges = sessions(s).edges(sessions(s).labels==3);                  % keep only left feature (right is complementary)
        X(3).patterns = sessions(s).patterns(sessions(s).labels==3);            % keep only left feature (right is complementary)
        X(3).labels = X(3).labels+((X(3).edges-1)*4);         
        X(3).ord = sessions(s).ord(sessions(s).labels==3);   
        % X4 = cued dimension is PATTERN      
        X(4).data = sessions(s).data(sessions(s).labels==4, :,:);
        X(4).labels = sessions(s).labels(sessions(s).labels==4);
        X(4).colors = sessions(s).colors(sessions(s).labels==4);                % keep only left feature (right is complementary)
        X(4).shapes = sessions(s).shapes(sessions(s).labels==4);                % keep only left feature (right is complementary)
        X(4).edges = sessions(s).edges(sessions(s).labels==4);                  % keep only left feature (right is complementary)
        X(4).patterns = sessions(s).patterns(sessions(s).labels==4);            % keep only left feature (right is complementary)
        X(4).labels = X(4).labels+((X(4).patterns-1)*4); 
        X(4).ord = sessions(s).ord(sessions(s).labels==4);   
        % catenate dimensions
        XX = [X(1).data;X(2).data;X(3).data;X(4).data];
        y = [X(1).labels;X(2).labels;X(3).labels;X(4).labels];
        ord = [X(1).ord;X(2).ord;X(3).ord;X(4).ord];
        [y,i] = sort(y,'ascend');
        X = XX(i,:,:);
        ord = ord(i);
        Xs(s).X = X;      % collect sessions in a structure
        ys(s).y = y;
        ords(s).ord = ord;

        % Preallocate trial-resolved results in the original experiment order.
        result.lda_prob_st(s).s = nan(n_meta, n_perm, length(data.labels_ses1), n_time, n_conditions);
        result.lda_acc_st(s).s = nan(n_meta, n_perm, length(data.labels_ses1), n_time);
        result.y_test(s).s = nan(n_meta,n_perm,length(data.labels_ses1));
    end
    
    for s=1:size(Xs,2)
        shuffle(s).sh = randperm(size(Xs(s).X,1));  % used when n_meta==2. Indices of shuffled trials - to ensure each fold contains items from each label   
    end
    for m=1:n_meta                                       
        for s=1:size(Xs,2)
            if n_meta==2 && m==2
                shuffle(s).sh = flip(shuffle(s).sh); % order as first shuffle, but at m=2 it is reversed. % one meta-perm is to ensure most trials enters training in some kfold. (Otherwise, with single shuffle, the trials beyond PC cap do not)
            elseif n_meta>2
                shuffle(s).sh = randperm(size(Xs(s).X,1)); % new random re-shuffle every time
            end
            Xs(s).X=Xs(s).X(shuffle(s).sh,:,:);
            ys(s).y=ys(s).y(shuffle(s).sh);
            ords(s).ord = ords(s).ord(shuffle(s).sh);
            ur(s).order(shuffle(s).sh) = ys(s).y; %check restoring order of labels
        end
        
    
        conditions = unique(ys(1).y); % take one 'session' as example
        for s=1:size(Xs,2)
            n_trials(s,:) = histc(ys(s).y, conditions);
        end
        % Use a common cap across PC sessions and classes to balance training.
        minPC = min(n_trials,[],'all');
        
        % get trial indices by PC and cond
        for s=1:size(Xs,2)
            for c = 1:length(conditions)
                    tmp_idx = 1:(length(ys(s).y));               % all of data indices for a 'session' (PC)
                    c_idx(s,c).idx(:) = tmp_idx(ys(s).y==conditions(c));        % store indices in a structure
            end
        end
    
        for f = 1:n_perm
            
            y_trainALL = [];
            X_trainALL = [];
            for s=1:size(Xs,2)
                y_tesT(s).y = [];
                X_test(s).X = [];
                ord_test(s).ord=[];
                for c = 1:length(conditions)                        % kfold tailored on condition-wise minPC
                    % Hold out the current fold from the balanced trial cap.
                    voidPC = ((floor(minPC/n_perm) * (f-1))+1 : (floor(minPC/n_perm) * f));
                    trainPC = find(~ismember(1:minPC,voidPC));
                    trainPCcap(s,c).idx = c_idx(s,c).idx(trainPC); % loop through cond (c) and PC (s), store indices capped 
                    y_trainALL = [y_trainALL, ys(s).y(trainPCcap(s,c).idx)'];        % composition of labels (8 cond, 4 ses) - already ascending within PCs (nb: PC are stacked)
                    X_trainALL = cat(1, X_trainALL, Xs(s).X(trainPCcap(s,c).idx,:,:));
                end
                %allocate to test - includes voidPC + others beyond PCmin 
                tmp_idx = 1:(length(ys(s).y));
                testPC(s).idx = find(~ismember(tmp_idx,[trainPCcap(s,:).idx])); % for each session, select indices not used for train (any cond)
                y_tesT(s).y = [y_tesT(s).y, ys(s).y(testPC(s).idx)'];           % composition of labels (8 cond, 4 ses) - already ascending (nb: PC are stacked)
                X_test(s).X = cat(1, X_test(s).X, Xs(s).X(testPC(s).idx,:,:));
                ord_test(s).ord = [ord_test(s).ord, ords(s).ord(testPC(s).idx)'];
                [y_tesT(s).y, i] = sort(y_tesT(s).y, 'ascend');
                X_test(s).X = X_test(s).X(i,:,:); 
                ord_test(s).ord = ord_test(s).ord(i);
                n_trials_test(s,:) = hist(y_tesT(s).y, conditions);
            end
            [y_trainALL, i] = sort(y_trainALL, 'ascend');
            X_trainALL = X_trainALL(i,:,:); % keep condition blocks contiguous
            n_trials_train = hist(y_trainALL, conditions);
    
            % Randomize trials independently within each condition. The
            % condition-block offsets preserve the established permutation scheme.
            fprintf('metaperm %g, kfold %g / %g\n', m, f, n_perm)
            train_indices = cell(1, n_conditions);
            for c1 = 1:n_conditions
                train_indices{c1} = randperm(n_trials_train(c1)) + ...
                    sum(n_trials_train(1:c1-1));
            end
            for s = 1:size(Xs,2)
                test_indices(s).idx = cell(1, n_conditions);
                for c1 = 1:n_conditions
                    test_indices(s).idx{c1} = randperm(n_trials_test(s,c1)) + ...
                        sum(n_trials_test(s,1:c1-1));
                end
            end

            X_train_permuted = X_trainALL([train_indices{:}], :, :);
            for s = 1:size(Xs,2)
                X_test_permuted(s).X = X_test(s).X([test_indices(s).idx{:}], :, :);
            end

            % Whiten using the Epoch covariance estimated from training data.
            sigma_conditions = reshape(y_trainALL([train_indices{:}]),[],1);
            sigma_ = nan(n_conditions, n_sensors, n_sensors);
            for c = 1:n_conditions
                tmp_ = nan(n_time, n_sensors, n_sensors);
                for t = 1:n_time
                    tmp_(t, :, :) = cov1para(X_train_permuted(sigma_conditions==c, :, t));
                end
                sigma_(c, :, :) = mean(tmp_, 1);
            end
            sigma = squeeze(mean(sigma_, 1));
            sigma_inv = sigma^-0.5;
            for t = 1:n_time
                X_train_permuted(:, :, t) = squeeze(X_train_permuted(:, :, t)) * sigma_inv;
                for s = 1:size(Xs,2)
                    X_test_permuted(s).X(:, :, t) = squeeze(X_test_permuted(s).X(:, :, t)) * sigma_inv;
                end
            end

            y_train = reshape(y_trainALL([train_indices{:}]),[],1);
            for t = 1:n_time
                fprintf('time %g\n', t)
                data_train = X_train_permuted(:, :, t);
                if isnan(gamma(t))
                    model_lda_alldata = fitcdiscr(data_train, y_train, ...
                        'OptimizeHyperparameters', 'Gamma', ...
                        'HyperparameterOptimizationOptions', struct('ShowPlots', 0));
                    gamma(t) = model_lda_alldata.Gamma;
                end
                model_lda = fitcdiscr(data_train, y_train, 'Gamma', gamma(t));

                for s = 1:size(Xs,2)
                    data_test = X_test_permuted(s).X(:, :, t);
                    y_test = reshape(y_tesT(s).y([test_indices(s).idx{:}]),[],1);
                    [lda.predictedlabels,lda.prob_estimates,lda.costs] = predict(model_lda,data_test);
                    for i = 1:n_conditions
                        for j = 1:n_conditions
                            result.lda_prob_full(s,m,f,i,j,t) = mean(lda.prob_estimates(y_test==conditions(i),conditions(j)),'omitnan');
                        end
                        result.lda_cond_acc(s,m,f,i,t) = mean(lda.predictedlabels(y_test==conditions(i))==conditions(i),'omitnan');
                    end
                    result.lda(s,m,f,t) = mean(lda.predictedlabels==y_test)-0.5;

                    % Map the within-condition permutation back to experiment order.
                    test_order = ord_test(s).ord([test_indices(s).idx{:}]);
                    result.lda_prob_st(s).s(m,f,test_order,t,:) = lda.prob_estimates;
                    result.lda_acc_st(s).s(m,f,test_order,t) = lda.predictedlabels==y_test;
                    result.y_test(s).s(m,f,test_order) = y_test;
                end
            end
        end
    end
    % average across repetitions and folds
    result_.lda_prob_full = squeeze(mean(result.lda_prob_full, [2,3],'omitnan'));
    result_.lda = squeeze(mean(result.lda, [2,3],'omitnan'));
    result_.lda_cond_acc = squeeze(mean(result.lda_cond_acc, [2,3],'omitnan'));
    result_.behav = behav;
    result_.gamma = gamma;
    for s=1:size(Xs,2)
        result_.lda_prob_st(s).s = squeeze(mean(result.lda_prob_st(s).s, [1,2],'omitnan'));
        result_.lda_acc_st(s).s = squeeze(mean(result.lda_acc_st(s).s, [1,2],'omitnan'));
        result_.y_test(s).s = mean(result.y_test(s).s,(1:2),'omitnan'); % consolidate labels across repetitions and folds
    end
    result = result_;

if local==0
    fname = sprintf(('/mnt/beegfs/home/gheza/EEG_data/RSA/distances_mltclass_sbj_%d.mat'), sn(si));
else
    fname = sprintf(('/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/preprocessed_data/rsa_data/distances_mltclass_sbj_%d.mat'), sn(si));
end
save (fname, 'result')


end
