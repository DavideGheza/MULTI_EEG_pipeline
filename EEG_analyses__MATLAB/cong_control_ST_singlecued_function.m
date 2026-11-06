% Run congruency_control decimation on RSA single-trial results.
% 07/02/2024 controlling for CNG within each cued dimension and PC (ses).
% based off 03/28/2023 cong_control_function 

function di = cong_control_ST_singlecued_function(behav,bootstrap_n)


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
PC_matrix = [Neu_1,MI,MC,Neu_2]; %nb: at difference than elsewhere, this maps onto PC order of sessions
% single trial fits are organized by session (PC) on the 3rd dimension. 
% here, to define decimated/bootstrapped indexes within PC, first identify the PC of the cued dim
cued_PC_matrix = PC_matrix==behav.task; % identify PC (column of PC_matrix) of the cued dim

for s=1:4
    for task=1:4
        di(s,task).bst = [];
    end
end
sum_cong = sum(behav.congruency,2)+1; 
cngidx = sum_cong==(0:3);


for s=1:4 % select within sessions (PC)

    for task=1:4
        
        for bst = 1:bootstrap_n

            sum_cong_PC = sum_cong(cued_PC_matrix(:,s) & behav.task==task); % cued_PC_matrix(:,s) indexes a given PC.
            hstcong = hist(sum_cong_PC,(0:3));
            mincong = min(hstcong);
            if mincong == 0
                fprintf("mincong is zero")
                di(s,task).bst = [];
            else
            % alternative: defining 2 mins (for mid and extreme cng) allows more samples, but oversampling on cng 1 and 2
    %         mincong_03 = min(hstcong([1,4])); % defining min for extreme cng (cng0idx and cng3idx)
    %         mincong_12 = min(hstcong([2,3])); % defining min for mid cng (cng1idx and cng2idx)
            
            % define indexes based on CNG and PC
            cng0idx = randsample(find(cngidx(:,1)==1 & cued_PC_matrix(:,s) & behav.task==task), mincong); % give indices over full array of all trials
            cng1idx = randsample(find(cngidx(:,2)==1 & cued_PC_matrix(:,s) & behav.task==task), mincong);
            cng2idx = randsample(find(cngidx(:,3)==1 & cued_PC_matrix(:,s) & behav.task==task), mincong);
            cng3idx = randsample(find(cngidx(:,4)==1 & cued_PC_matrix(:,s) & behav.task==task), mincong);
        
            decimated_indices_bst = [cng0idx;cng1idx;cng2idx;cng3idx]; % will be less than original trial # - some are dropped to match congruency
            di(s,task).bst = cat(1,di(s,task).bst,decimated_indices_bst); % catenate the new bst iteration
            end
        end

    end

end



end