%% Estimate cross-condition similarity between full-model encoding axes
% Uses cross-validation and covariance-weighted cosine similarity on the
% original, pre-orthogonalization axes. Outputs support the location and
% conflict similarity matrices in Figure 7C and SI Figure S5.
% Author: Davide Gheza

function [result] = encoding_fullmodel_repCV__similarity_HPC(si)

local=1;

%% 1. Load settings


%Select onset or response events
onset = false;
if onset
    suffix = '';
else
    suffix = 'response_';   % response-locked requires precomputing stimulus-locked to get a fixed pre-stimulus baseline
end


%%
sn = [22,24:48];

if local==0
    d = load(['/scratch/gheza/EEG_data/sessions_data_' suffix 'tf_decoding_BN_sbj_'...% for tf!       If time-domain, delete 'tf_decoding_BN_'
            num2str(sn(si))...
            '.mat']);
    behav = load(['/scratch/gheza/EEG_data/behav/V3_behav_seqinfo_trimmed_' suffix 'sbj_' num2str(sn(si))]);    % load behavioral data - Updated to latest version, to use in extracting trial aggregates (should not matter for main encoding part - 12/16/2025)
else
    root = '/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/preprocessed_data/sessions_data_tf/'; % for local
    d=load(fullfile(root, ['sessions_data_' suffix 'tf_decoding_BN_sbj_' num2str(sn(si)) '.mat']));                                                             % for local
    behav = load(['/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/preprocessed_data/sessions_data/V3_behav_seqinfo_trimmed_' suffix 'sbj_' num2str(sn(si))]); % for local
end
data = d.data;
behav = behav.behav_trim;

% field names - to be used to define design matrix (encoding model)
fns = fieldnames(behav);
fns = fns(9:end); % dimensions' l r labels
% mapping between labels (1:8, with 1:4=left and 5:8=right) with RSA model matrix (L-R color, L-R shape, L-R edge, L-R pattern). (encoding model)
column_map = [1,5,2,6,3,7,4,8];
% define congruency axes (dimension specific conflict regressors)
all_conflict_regs = behav.congruency;
all_conflict_regs = all_conflict_regs + (all_conflict_regs == -1) - (all_conflict_regs == 0); % set cued dimension as 0, incongruent as -1, congruent as 1
all_prev_conflict_regs = behav.prev_congruency; % note: do not just lag 'all_conflict_regs' because behav and data have been reordered blockwise. Use behav prev_congruency instead
all_prev_conflict_regs = all_prev_conflict_regs + (all_prev_conflict_regs == -1) - (all_prev_conflict_regs == 0); % set cued dimension as 0, incongruent as -1, congruent as 1
all_prev_conflict_regs(isnan(all_prev_conflict_regs)) = 0;  % previous congruency for 1st trial: set to zero (in principle it's NaN, and we should remove trial from encoding data. However, since contrast coded, zero means no contrast contribution on this trial)
% define cue dummy regressors (4 cols c,s,e,p, taking value 1 if dimension is cued, 0 otherwise)
all_cue_regs = (1:4)==behav.task;
% to check if any fold has ever a singular design
sing_counter = 0;

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
    full_behav = load(['/scratch/gheza/EEG_data/behav/dimensions_PC_s' num2str(sn(si)) '_0' num2str(sn(si))]);
else
    sbj = dir('/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/data/EEG_data/*PC*'); % for local
    full_behav = load(['/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/data/EEG_data', '/', sbj(sn(si)).name]);
end
correct_feature = full_behav.data.correct_feature;
% force side coding - align SIDE coding across Color, Edge, Shape, and Pattern
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
n_conditions = length(unique(data.labels_ses1))*2; % cue * side = 8 conditions
n_sensors = size(data.ses1, 2);
n_time = size(data.ses1, 3);
n_sessions = length(sessions);
n_perm = 6;  % number of kfold 
n_meta = 8; % number of repetition of k-fold, after random shuffling

% Cross-validated encoding steps (for each permutation step)
% 1 Randomize single trials within class for the training and test data.
% 2 Whiten the training and test data using the Epoch method.
% For each timepoint:
% 3 Encode axes with training data (EEG features)
% 4 Encode axes also with test data (combine trials from all PCs)
% 5 Compute pairwise axes similarity (cross-condition alignment)

% preallocate
m_f_count = 0; % counter of repeated-CV iterations
similarity = zeros(16,16,n_time); % similarity matrix (16 axes)

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
    minPC = min(n_trials,[],'all');% balance the training set across proportion-congruency conditions and classes
    
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
        ord_train.ord=[];
        for s=1:size(Xs,2)
            y_tesT(s).y = [];
            X_test(s).X = [];
            ord_test(s).ord=[];
            for c = 1:length(conditions)                        % kfold tailored on condition-wise minPC
                %allocate to train - to control for PC # (and condition #, if singletrial approach), a subselection of trials = minPC goes into train
                voidPC = ( (floor(minPC/n_perm) * (f-1))+1 : (floor(minPC/n_perm) * (f)) );  
                trainPC = find(~ismember((1:minPC),voidPC));           % k-fold index for cond indexes
                
                trainPCcap(s,c).idx = c_idx(s,c).idx(trainPC); % loop through cond (c) and PC (s), store indices capped 
                y_trainALL = [y_trainALL, ys(s).y(trainPCcap(s,c).idx)'];        % composition of labels (8 cond, 4 ses) - already ascending within PCs (nb: PC are stacked)
                X_trainALL = cat(1, X_trainALL, Xs(s).X(trainPCcap(s,c).idx,:,:));
                ord_train.ord = [ord_train.ord, ords(s).ord(trainPCcap(s,c).idx)'];
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
        X_trainALL = X_trainALL(i,:,:);                                % sort trials by class before within-class permutation
        ord_train.ord = ord_train.ord(i);
        n_trials_train = hist(y_trainALL, conditions);

        % Randomize single trials independently within each class.
        fprintf('metaperm %g, kfold %g / %g\n', m, f, n_perm)
        train_indices = cell(1, n_conditions);
        for c1 = 1:n_conditions
            train_indices{c1} = randperm(n_trials_train(c1)) + sum(n_trials_train(1:c1-1));
        end

        for s = 1:size(Xs,2)
            test_indices(s).idx = cell(1, n_conditions);
            for c1 = 1:n_conditions
                test_indices(s).idx{c1} = randperm(n_trials_test(s,c1)) + sum(n_trials_test(s,1:c1-1));
            end
        end

        % Apply the within-class permutations without averaging trials.
        X_train_permuted = X_trainALL([train_indices{:}], :, :);
        for s = 1:size(Xs,2)
            X_test_permuted(s).X = X_test(s).X([test_indices(s).idx{:}], :, :);
        end

            % 2. Whitening using the Epoch method    
            sigma_conditions = reshape(y_trainALL([train_indices{:}]),[],1); %  single trial labels
            sigma_ = nan(n_conditions, n_sensors, n_sensors);
            for c = 1:n_conditions
                % compute sigma for each time point, then average across time
                tmp_ = nan(n_time, n_sensors, n_sensors);
                for t = 1:n_time
                    tmp_(t, :, :) = cov1para(X_train_permuted(sigma_conditions==c, :, t));
                end
                sigma_(c, :, :) = mean(tmp_, 1);
            end
            sigma = squeeze(mean(sigma_, 1));  % average across conditions
            sigma_inv = sigma^-0.5;
            for t = 1:n_time
                X_train_permuted(:, :, t) = squeeze(X_train_permuted(:, :, t)) * sigma_inv;
                for s = 1:size(Xs,2)
                    X_test_permuted(s).X(:, :, t) = squeeze(X_test_permuted(s).X(:, :, t)) * sigma_inv;
                end
            end
            





            for t = 1:n_time % TRAIN time
                fprintf('time %g\n', t)
                % 3. Fit the classifier using training data
                data_train = X_train_permuted(:, :, t);

                y_train = y_trainALL([train_indices{:}]);
                y_train = reshape(y_train,[],1);

                % build X matrix of regressors 
                % cued side regressors
                side_regs = (y_train == (column_map)) .* (2 * ((column_map)>4) - 1); 
                % side_regs = [ones(length(y_train),1), (side_regs(:,1)+side_regs(:,2)), (side_regs(:,3)+side_regs(:,4)), (side_regs(:,5)+side_regs(:,6)), (side_regs(:,7)+side_regs(:,8))]; % added intercept
                side_regs = [(side_regs(:,1)+side_regs(:,2)), (side_regs(:,3)+side_regs(:,4)), (side_regs(:,5)+side_regs(:,6)), (side_regs(:,7)+side_regs(:,8))]; % no intercept
                % confict regressors (congruency of each dimension)
                train_order = ord_train.ord([train_indices{:}]); % recover original trial indices (to match behavioral file)
                train_order = train_order';
                conflict_regs = all_conflict_regs(train_order,:);
                p_conflict_regs = all_prev_conflict_regs(train_order,:);
                % cue dummies regressors
                cue_regs = all_cue_regs(train_order,:);
                % assemble and inspect rank deficiency:
                xmat = [cue_regs,side_regs,conflict_regs,p_conflict_regs];    
                rankX = rank(xmat);
                if rankX < size(xmat,2)
                    [U,S,V] = svd(xmat,'econ');
                    singVals = diag(S);  
                    warning('Design matrix is rank-deficient!');
                    sing_counter = sing_counter+1;
                end
                % fit encoding model on first fold (train data) using the normal equation: B = (X' * X) \ (X' * Y)
                B1 = (xmat' * xmat) \ (xmat' * data_train);
                data1_fitted = xmat * B1;
                res1 = data_train-data1_fitted;
                % remove mean from cue regressors (one-hot coded) 
                B1(1:4,:) = B1(1:4,:) - mean(B1(1:4,:), 1); % Remove shared mean across cue regressors to eliminate intercept-like (common-mode) variance that would artificially inflate similarity

                % fit encoding model on second fold (test data). 
                % Aggregate trials from all PCs subsets
                ymat = [];
                all_data_test = [];
                for s=1:size(Xs,2)
                    data_test = X_test_permuted(s).X(:, :, t);
                    y_test = y_tesT(s).y([test_indices(s).idx{:}]);
                    y_test = reshape(y_test,[],1);
                    test_order = ord_test(s).ord([test_indices(s).idx{:}]); % recover trial indices to pull all features (cued and non-cued) from behavioral file
                    % build Y matrix of regressors 
                    % cued side regressors
                    side_regs = (y_test == (column_map)) .* (2 * ((column_map)>4) - 1); 
                    side_regs = [(side_regs(:,1)+side_regs(:,2)), (side_regs(:,3)+side_regs(:,4)), (side_regs(:,5)+side_regs(:,6)), (side_regs(:,7)+side_regs(:,8))]; % no intercept
                    % confict regressors (congruency of each dimension)
                    conflict_regs = all_conflict_regs(test_order,:);
                    p_conflict_regs = all_prev_conflict_regs(test_order,:);
                    % cue dummies regressors
                    cue_regs = all_cue_regs(test_order,:);
                    % assemble and inspect rank deficiency:
                    tmp = [cue_regs,side_regs,conflict_regs,p_conflict_regs];
                    ymat = [ymat;tmp];
                    all_data_test = [all_data_test; data_test];
                end
                rankY = rank(ymat);
                if rankY < size(ymat,2)
                    [U,S,V] = svd(ymat,'econ');
                    singVals = diag(S);  
                    warning('Design matrix is rank-deficient!');
                    sing_counter = sing_counter+1;
                end    
                % fit
                B2 = (ymat' * ymat) \ (ymat' * all_data_test);
                data2_fitted = ymat * B2;
                res2 = all_data_test-data2_fitted;
                B2(1:4,:) = B2(1:4,:) - mean(B2(1:4,:), 1); % Remove shared mean across cue regressors to eliminate intercept-like (common-mode) variance that would artificially inflate similarity


                % estimate a single noise covariance matrix
                % first, concatenate the two sets of residuals observation-wise 
                residuals = [res1;res2];
                % covariance estimation, shrinkage, and inversion
                % W_est = inv(cov1para(residuals));  % Ledoit Wolf to shrink/regularize covariance, in case matrix is rank deficient
                % BBt = B1 * W_est * B2';
                Sigma = cov1para(residuals);  % Avoid explicitly computing inverse
                % whiten and estimate similarity
                BBt = B1 * (Sigma \ B2');           % raw similarity between pairs of axes
                % Normalize by axes length (Mahalanobis cosine similarity)
                L1 = sqrt(diag(B1 * (Sigma \ B1'))); % Vector. Norm (“length”) of each axis in B1 (in whitened space)
                L2 = sqrt(diag(B2 * (Sigma \ B2')));
                den = (L1*L2');       % Outer product. Matrix of pairwise norm products
                den(den == 0) = NaN;  % in case vectors have zero magnitude (not decodable)
                cosSigma = BBt ./ den; % Element-wise. Every pair of axes is normalized independently
                similarity(:,:,t) = similarity(:,:,t) + cosSigma;     

            end
            m_f_count = m_f_count + 1;
    end
end
% divide by number of permutations
encoded_similarity = similarity ./ m_f_count; % n_axes, n_axes, times (n_meta, n_perm are averaged inside loop!)

if local==0
    fname = sprintf((['/scratch/gheza/RSA_fits/encoding_fullmodel_similarity_' suffix 'sbj_%d.mat']), sn(si));
else
    fname = sprintf((['/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/preprocessed_data/rsa_data/encoding_fullmodel_similarity_' suffix 'sbj_%d.mat']), sn(si));
end

save(fname, 'encoded_similarity'); 
save(fname, 'sing_counter', '-append');


% example design matrix
% domains  = {'cue','cuedLoc','conflict','prevConflict'};
% features = {'color','shape','edge','pattern'};
% headers = [];
% for d = 1:numel(domains)
%     for f = 1:numel(features)
%         headers{end+1} = sprintf('%s_%s', domains{d}, features{f});
%     end
% end
% 
% xmat_tbl = array2table(xmat, 'VariableNames', headers);
% writetable(xmat_tbl, 'design_matrix.csv');
