% Run congruency_control decimation on RSA single-trial results.
% 07/02/2024 controlling for CNG within each distractor and PC
% first parse each distractor by congruent vs incongruent 
% based off 03/28/2023 cong_control_function 
% 07/25/2024 Added control for CNG of EACH distractor, within each dist.
% Therefore, integrated with simpler cong_control_ST_function.
% Summary: 1) balance cng by PC session. 2) take bootstrapped
% indices by s, and subselect to ensure exact balance of cng within distractor.

function [di,diSD] = cong_control_ST_singledstr_function(behav,bootstrap_n)

% vectorize
task = behav.task;
congruency = behav.congruency;
blkwise_PCs = behav.blkwise_PCs;

% define vectors for indexing trials by PC
Neu_1_2 = nan(length(blkwise_PCs),2);     %preallocate - two neutral dimensions per trial
MI = nan(length(blkwise_PCs),1);   
MC = nan(length(blkwise_PCs),1);
% find task based on PC
for t = 1:length(blkwise_PCs)      
    Neu_1_2(t,:) = find(blkwise_PCs(t,:) == 1/2);
    MI(t,:) = find(blkwise_PCs(t,:) == 1/3);
    MC(t,:) = find(blkwise_PCs(t,:) == 2/3);
    if Neu_1_2(t,1) == 1 && Neu_1_2(t,2) == 4               % Neu_1 and Neu_2 are artificially segregated as neutral dimension on left vs right column (of blkwise_PCs). Hence, if Neu1 happens to be contiguous to Neu2, Neu1 never contains cond 4, and Neu2 never contains cond 1.
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
PC_matrix = [Neu_1,MI,MC,Neu_2]; %nb: at difference than elsewhere, this maps onto PC order of sessions
% single trial fits are organized by session (PC) on the 3rd dimension. 
% here, to define decimated/bootstrapped indexes within PC, first identify the PC of the cued dim
cued_PC_matrix = PC_matrix==task; % identify PC (column of PC_matrix) of the cued dim

% Preallocate structure arrays
di = struct('bst', cell(1, 4));
diSD = struct('bst', cell(4, 4));

sum_cong = sum(congruency,2)+1; 
cngidx = sum_cong==(0:3);

for bst = 1:bootstrap_n

    for s=1:4 % select within sessions (PC)


        sum_cong_PC = sum_cong(cued_PC_matrix(:,s)); % cued_PC_matrix(:,s) indexes a given cued task when assigned to the PC(s)
        hstcong = hist(sum_cong_PC,(0:3));
        mincong = min(hstcong);
        
        % define indexes based on CNG and PC
        idx = arrayfun(@(x) randsample(find(cngidx(:, x)==1 & cued_PC_matrix(:,s)), mincong), 1:4, 'UniformOutput', false);  % give indices over full array of all trials
        decimated_indices_byses = vertcat(idx{:});        % will be less than original trial # - some are dropped to match congruency
        di(s).bst = [di(s).bst; decimated_indices_byses]; % catenate the new bst iteration - store for simpler bst by PC


        for dist=1:4
            
            dec_indices = decimated_indices_byses(congruency(decimated_indices_byses, dist) == 0 | congruency(decimated_indices_byses, dist) == 1);
            dist_x_cng = congruency(dec_indices, dist) == 1;

            hstcong = hist(dist_x_cng,(0:1));
            mincong = min(hstcong);
            
            % define indexes based on CNG and PC
            cng0id = dec_indices(congruency(dec_indices,dist)==0);
            cng1id = dec_indices(congruency(dec_indices,dist)==1);
            cng0id = cng0id(1:mincong);
            cng1id = cng1id(1:mincong);
        
            decimated_indices_bst = [cng0id;cng1id]; % will be less than original trial # - some are dropped to match congruency
            diSD(s,dist).bst = [diSD(s, dist).bst; decimated_indices_bst]; % catenate the new bst iteration
            
            if bst==bootstrap_n % if last bst
                diSD(s,dist).dist_inc = diSD(s,dist).bst(congruency(diSD(s,dist).bst,dist)==0);
                diSD(s,dist).dist_con = diSD(s,dist).bst(congruency(diSD(s,dist).bst,dist)==1);
            end
        end

     end

end

% mean(sum_cong(di(1).bst)) % must equal 1.5 (mean(0:3))
% mean(sum_cong(di(2).bst)) % must equal 1.5 (mean(0:3))
% mean(sum_cong(di(3).bst)) % must equal 1.5 (mean(0:3))
% mean(sum_cong(di(4).bst)) % must equal 1.5 (mean(0:3))

end