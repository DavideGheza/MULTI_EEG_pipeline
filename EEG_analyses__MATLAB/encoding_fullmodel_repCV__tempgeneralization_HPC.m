%% Fit the full encoding model with temporal generalization
% Uses repeated cross-validation to estimate axes at every training time and
% project them onto held-out data at every test time. Dimension-private axes
% are orthogonalized from the other dimensions; dimension-shared axes use
% cross-condition projection. Two-stage whitening is applied to raw power and
% time-resolved encoding residuals. Outputs support Figure 7B.
% Author: Davide Gheza

function [result] = encoding_fullmodel_repCV__tempgeneralization_HPC(si)

local=0;

%% 1. Load settings

raw_whitening_on = 1; % if 0, avoid whitening the raw data (default for encoding models: 1) (if toggle off, then ideally first z-score raw data before encoding)
%Select onset or response events
onset = true;
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
all_prev_conflict_regs = behav.prev_congruency;
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
% 3 Encode cued dimension-specific side with training data (EEG features)
% 4 Regress out dimension-general side axis (i.e. QR(side,dim-spec-side-axes)) and dimension-general conflict
% 5 Project residuals on test data (segregated by PC)

% preallocate
m_f_count = 0; % counter of repeated-CV iterations

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

    % Preallocate single-trial results.
    encoded_cue(s).s             = zeros(length(data.labels_ses1), n_time, n_time, 4, 'single'); % trials, train times, test times, and 4 axes; repeated-CV iterations are averaged in the loop
    encoded_dim_spec_SR(s).s     = zeros(length(data.labels_ses1), n_time, n_time, 4, 'single');
    encoded_dim_gen_SI(s).s      = zeros(length(data.labels_ses1), n_time, n_time, 4, 'single');
    encoded_dim_spec_conf(s).s   = zeros(length(data.labels_ses1), n_time, n_time, 4, 'single');
    encoded_dim_gen_conf(s).s    = zeros(length(data.labels_ses1), n_time, n_time, 4, 'single');
    encoded_prev_conf(s).s       = zeros(length(data.labels_ses1), n_time, n_time, 4, 'single');
    n_tested(s).c                = zeros(length(data.labels_ses1), 1); % FIXED: number of times each trial is tested
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

        if raw_whitening_on
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
                % [U,S,V] = svd(xmat,'econ');
                % singVals = diag(S);  
                % warning('Design matrix is rank-deficient!');
                sing_counter = sing_counter+1;
            end
            % GLM: fit no-intercept encoding model on train data - multivariate linear regression using the normal equation: B = (X' * X) \ (X' * Y)
            B = (xmat' * xmat) \ (xmat' * data_train);
            data1_fitted = xmat * B;
            res = data_train-data1_fitted;
            % remove mean from cue regressors (one-hot coded) 
            B(1:4,:) = B(1:4,:) - mean(B(1:4,:), 1); % Remove shared mean across cue regressors to eliminate intercept-like (common-mode) variance that would artificially inflate similarity
            
            % Estimation of dimension-specific and dimension-general axes
            
            % side axes (spatial attention)
            % preallocate
            B_side = B(5:8,:);  % 4 x n_features
            B_resid_s = zeros(size(B_side));
            B_composite_heldout = zeros(size(B_side));    % for general. Store held-out 3-d composite
            for i = 1:4 % iterate over dims
                rest_idx = setdiff(1:4, i);
                rest = B_side(rest_idx, :);  % 3 x n_features
                B_composite_heldout(i,:) = mean(B_side(rest_idx,:), 1); % for general: 3-d composite
                
                [U, ~, ~] = svd(rest', "econ");  % U: n_features x 3, orthonormal basis of rest subspace
                B_resid_s(i,:) = B_side(i,:) - B_side(i,:) * (U * U'); % remove rest subspace from axis i via hat matrix
                % or equivalent with regression: (rest is the new 3-d composite vector, replacing B_bar (4-d) in previous code)
                % beta = (rest * rest') \ (rest * B_side(i,:)'); % fit B_side(i,:)' as linear combination of rest axes -> 3 scalar betas
                % B_side_specific2(i,:) = (B_side(i,:)' - rest' * beta); % to get orthogonal, subtract the fitted (rest' * beta) that spans only the rest subspace
            end
            B(5:8,:) = B_resid_s;% replace original with residualized axes in full B matrix 
            B(17:20,:) = B_composite_heldout; % add composite heldout for each dimension to estimate 'dimension-general attention axis'

            % conflict axes (apply same residualization)
            % preallocate
            B_conf = B(9:12,:);  % 4 x n_features
            B_resid_c = zeros(size(B_conf));
            B_composite_heldout = zeros(size(B_conf));    % for general. Store held-out 3-d composite
            for i = 1:4 % iterate over dims
                rest_idx = setdiff(1:4, i);
                rest = B_conf(rest_idx, :);  % 3 x n_features
                B_composite_heldout(i,:) = mean(B_conf(rest_idx,:), 1); % for general: 3-d composite
                
                [U, ~, ~] = svd(rest', "econ");  % U: n_features x 3, orthonormal basis of rest subspace
                B_resid_c(i,:) = B_conf(i,:) - B_conf(i,:) * (U * U'); % remove rest subspace from axis i via hat matrix
            end
            B(9:12,:) = B_resid_c;% replace original with residualized axes in full B matrix 
            B(21:24,:) = B_composite_heldout; % add composite heldout for each dimension to estimate 'dimension-general conflict axis'
            
            % prewhiten the axes using the residuals of the initial GLM
            % W = inv(cov1para(res));  % Ledoit Wolf to shrink/regularize covariance, in case matrix is rank deficient
            % B_whitened  = W * B';
            B_whitened  = cov1para(res) \ B'; % A\b faster than inv(A)*b




            % project held-out data into axes
            for s=1:size(Xs,2)
                % 4. Compute and store classification accuracies
                y_test = y_tesT(s).y([test_indices(s).idx{:}]);
                y_test = reshape(y_test,[],1);
                test_order = ord_test(s).ord([test_indices(s).idx{:}]); % recover trial indices to pull all features (cued and non-cued) from behavioral file
                if t == 1 % FIXED: count each test once per (repetition, fold)
                    n_tested(s).c(test_order) = n_tested(s).c(test_order) + 1; % FIXED
                end % FIXED
                % X matrix of cued AND non-cued side regressors, to rectify projections based on expected sign per dimension
                dim_lbls = nan(length(y_test),4);
                for dim = 1:4
                    dim_lbls(:,dim) =  dim + (behav.(fns{dim})(test_order,1)-1) * 4;
                end
                dim_specific_side_x = (2 * (dim_lbls>4) - 1); % note: it overlaps with the non-zero elements of side_regs computed on test labels (which would be the cued-only dim.specific side design matrix)
                % for conflict: rectify projections based on expected sign per dimension (conflict coding applies only to non-cued dims - cued have 0s)
                dim_specific_conflict_x = all_conflict_regs(test_order,:); % current conflict
                dim_specific_p_conflict_x = all_prev_conflict_regs(test_order,:); % previous conflict
                % for cue: rectify cued only
                dim_specific_cue_x = all_cue_regs(test_order,:);

                % Temporal genealization: test for same time t as training, but also all the other times t2 
                for t2=1:n_time % TEST time
                    data_test = X_test_permuted(s).X(:, :, t2);        % CORRECT time tested will be on diagonal of t,t2 slice

                    % projection & rectify (all of the following are of size: nr trials x 4 axes)
                    proj = @(Bcols, gate) (data_test * B_whitened(:,Bcols)) .* gate; % helper function

                    projection_cue              = proj(1:4,   dim_specific_cue_x); % cue is one-hot coded.
                    projection_dim_spec_SR      = proj(5:8,   dim_specific_side_x);% dimension-specific sr 
                    projection_dim_spec_conflict= proj(9:12,  dim_specific_conflict_x);
                    projection_prev_conflict    = proj(13:16, dim_specific_p_conflict_x);
                                            
                    projection_dim_gen_side     = proj(17:20, dim_specific_side_x);    % composite heldout for each dimension, to estimate dimension-general attention axis ..contrast-gated into dimension-specific readouts
                    projection_dim_gen_conflict = proj(21:24, dim_specific_conflict_x);% composite heldout for each dimension, to estimate dimension-general conflict axis ..contrast-gated into dimension-specific readouts

                    % collect single trials, place in order as from original data (segregated by session); 
                    % Repeated-CV and fold dimensions are averaged on the fly.
                    test_trial_n = size(projection_cue,1);
                    axes_n =  size(projection_cue,2); % WARNING: this assumes all projection sets have same n of axes (4)
                    encoded_cue(s).s(test_order,t,t2,:) = encoded_cue(s).s(test_order,t,t2,:) + reshape(projection_cue, [test_trial_n,1,1,axes_n]); 
                    encoded_dim_spec_SR(s).s(test_order,t,t2,:) = encoded_dim_spec_SR(s).s(test_order,t,t2,:) + reshape(projection_dim_spec_SR, [test_trial_n,1,1,axes_n]); 
                    encoded_dim_gen_SI(s).s(test_order,t,t2,:) = encoded_dim_gen_SI(s).s(test_order,t,t2,:) + reshape(projection_dim_gen_side, [test_trial_n,1,1,axes_n]); 
                    encoded_dim_spec_conf(s).s(test_order,t,t2,:) = encoded_dim_spec_conf(s).s(test_order,t,t2,:) + reshape(projection_dim_spec_conflict, [test_trial_n,1,1,axes_n]);
                    encoded_dim_gen_conf(s).s(test_order,t,t2,:) = encoded_dim_gen_conf(s).s(test_order,t,t2,:) + reshape(projection_dim_gen_conflict, [test_trial_n,1,1,axes_n]); 
                    encoded_prev_conf(s).s(test_order,t,t2,:) = encoded_prev_conf(s).s(test_order,t,t2,:) + reshape(projection_prev_conflict, [test_trial_n,1,1,axes_n]); 

                    % indexing explainer:   PC 'sessions' populate the structure (different trial n / session) - meta-perm - perm (kfold) - trial - time - probability of each class
                    % keep meta m and kfold f indices (redundancy of trials across kfolds - test includes voidPC + others beyond PCmin). Then average over f and m (omit nan)
                    % RN keeping sessions apart, each with full experiment trials (NAN out of session), but could drop .s and save all in single result matrix (in principle).
                end

      
            end

        end
        m_f_count = m_f_count + 1;
    end
end
% divide by number of permutations
for s=1:size(Xs,2)
    % FIXED: per-trial mean over the iterations in which the trial was tested
    % (untested trials, including all trials of other PCs, become NaN)
    encoded_cue(s).s = encoded_cue(s).s ./ n_tested(s).c; % FIXED
    encoded_dim_spec_SR(s).s = encoded_dim_spec_SR(s).s ./ n_tested(s).c; % FIXED: trials, train times, test times, and 4 axes
    encoded_dim_gen_SI(s).s = encoded_dim_gen_SI(s).s ./ n_tested(s).c; % FIXED
    encoded_dim_spec_conf(s).s = encoded_dim_spec_conf(s).s ./ n_tested(s).c; % FIXED
    encoded_dim_gen_conf(s).s = encoded_dim_gen_conf(s).s ./ n_tested(s).c; % FIXED
    encoded_prev_conf(s).s = encoded_prev_conf(s).s ./ n_tested(s).c; % FIXED
end

if local==0
    fname = sprintf((['/scratch/gheza/RSA_fits/encoding_fullmodel_' suffix 'sbj_%d.mat']), sn(si));
else
    fname = sprintf((['/Users/davidegheza/Library/CloudStorage/OneDrive-Personal/Documenti/MATLAB/dg_workspace22/dimensions/analyses/preprocessed_data/rsa_data/encoding_fullmodel_' suffix 'sbj_%d.mat']), sn(si));
end

% save(fname, 'encoded_dim_spec_SR', '-v7.3'); # very disk-heavy (~3gb total)
% save(fname, 'encoded_dim_gen_SI', '-append');
% save(fname, 'encoded_dim_spec_conf', '-append');
% save(fname, 'encoded_dim_gen_conf', '-append');


%% extract and save trial averages

save(fname, 'behav') % to initialize the save

ax_name = {'encoded_cue', 'encoded_dim_spec_SR', 'encoded_dim_gen_SI', 'encoded_dim_spec_conf', 'encoded_dim_gen_conf','encoded_prev_conf'};

% Store previous-congruency stratification for each exported axis.
% Each field below will contain prevI, prevC, their contrast, and trial counts.
prevc_by_axis = struct();

for pro_id = 1:length(ax_name)

    projections = eval(ax_name{pro_id});
    projections = cat(5,projections(1).s,... % cat the 4 PCs 
                                   projections(2).s,...
                                   projections(3).s,...
                                   projections(4).s);
    projections = permute(projections,[1,5,2,3,4]); % trial, pc, train time, test time, dimension
    
    % These exports keep PC conditions separate; reported encoding results use
    % the neutral-PC columns (selected in the plotting scripts), where bootstrapping
    % is unnecessary. To explore PC-biased (MI/MC) conditions, uncomment the
    % bootstrap branch below (requires `bst`, `bootstrap_n`, and
    % cong_control_ST_singlecued_function.m).
    % if ~bst
    idxs = 1:length(behav.task); % index all the trials without bootstrap
    % else    
    %     % bootstrap results to match congruency within each cued task and PC
    %     idxs_bypc = cong_control_ST_singlecued_function(behav,bootstrap_n); 
    %     idxs = vertcat(idxs_bypc(:).bst)'; % catenate indexes from all PCs
    % end
    
    % task_count = behav.taskcount(idxs);
    % congruency = behav.congruency(idxs,:);
    % prev_congruency = behav.prev_congruency(idxs,:);
    % post_correct = behav.post_correct(idxs);
    % prev_sum_con = behav.prev_sum_con(idxs);
    task = behav.task(idxs);
    congruency = behav.congruency(idxs,:);
    prev_congruency = behav.prev_congruency(idxs,:);
    task_count = behav.taskcount(idxs);
    post_correct = behav.post_correct(idxs);
    
    % parse cued from non-cued projections
    proj_sel = projections(idxs,:,:,:,:);   % [nTtr × pc × train × test × dist]
    nTtr = size(proj_sel,1); % if bootstrapping: trial n from bootstrapping selection, else simply number of trials

    cued_proj    = nan(nTtr, size(proj_sel,2), size(proj_sel,3), size(proj_sel,4));
    noncued_proj = proj_sel;
    
    for t = 1:nTtr
        d = task(t);
        cued_proj(t,:,:,:)        = proj_sel(t,:,:,:, d); % trial, pc, train time, test time
        noncued_proj(t,:,:,:,d)  = NaN;                  % trial, pc, train time, test time, dist
    end

    % To test whether previous trial state modulates current conflict-axis magnitude:
    % Additional export for each axis type:
    % stratify current projections by previous congruency on each distractor
    % (prevI minus prevC), using basic CSE filters.
    n_pc = size(noncued_proj,2);
    n_train = size(noncued_proj,3);
    n_test = size(noncued_proj,4);
    n_dims = size(noncued_proj,5);

    prevI = nan(n_pc, n_train, n_test, n_dims);
    prevC = nan(n_pc, n_train, n_test, n_dims);
    trial_counts = nan(2, n_dims); % row 1: prevI count, row 2: prevC count

    for dist = 1:n_dims
        idx_prevI = prev_congruency(:,dist)==0 & congruency(:,dist)~=-1 & task_count>1 & post_correct==1;
        idx_prevC = prev_congruency(:,dist)==1 & congruency(:,dist)~=-1 & task_count>1 & post_correct==1;

        prevI(:,:,:,dist) = squeeze(mean(noncued_proj(idx_prevI,:,:,:,dist),1,'omitnan'));
        prevC(:,:,:,dist) = squeeze(mean(noncued_proj(idx_prevC,:,:,:,dist),1,'omitnan'));

        trial_counts(1,dist) = sum(idx_prevI);
        trial_counts(2,dist) = sum(idx_prevC);
    end

    prevc_by_axis.(ax_name{pro_id}).prevI = prevI;
    prevc_by_axis.(ax_name{pro_id}).prevC = prevC;
    prevc_by_axis.(ax_name{pro_id}).contrast_prevI_minus_prevC = prevI - prevC;
    prevc_by_axis.(ax_name{pro_id}).trial_counts = trial_counts;



    cued_proj = squeeze(mean(cued_proj,1,'omitnan'));          % trial average
    noncued_proj = squeeze(mean(noncued_proj,1,'omitnan'));
        
    cued.(ax_name{pro_id}) = cued_proj;
    save(fname, 'cued', '-append');
    noncued.(ax_name{pro_id}) = noncued_proj;
    save(fname, 'noncued', '-append');
    

end

% Save stratified summaries for all axis types in one struct.
save(fname, 'prevc_by_axis', '-append');


