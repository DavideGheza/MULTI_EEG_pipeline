% Generate matched bootstrap indices for single-trial CSE analyses.
% Applies random decimation within proportion-congruency condition while matching current
% and previous congruency. The 16 previous-by-current combinations are pooled
% into four groups to retain adequate trial counts.
% Author: Davide Gheza

function di = cong_control_ST_precur_function(behav,bootstrap_n)


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

di = struct('bst', cell(1, 4));

sum_cong = sum(behav.congruency,2)+1; 
prev_sum_con = behav.prev_sum_con;
cngidx = sum_cong==(0:3);
p_cngidx = prev_sum_con==(0:3);
exclude_trials = (behav.post_correct==0 | behav.taskcount==1);


for bst = 1:bootstrap_n

    for s=1:4 % select within sessions (PC)
        

        sum_cong_PC = sum_cong(cued_PC_matrix(:,s) & ~exclude_trials); % cued_PC_matrix(:,s) indexes a given PC.        % trim excluded trials
        prev_sum_con_PC = prev_sum_con(cued_PC_matrix(:,s) & ~exclude_trials);
        
        hstcong = zeros(4, 4);
        for i = 1:4
            for e = 1:4
                hstcong(i,e) = sum(sum_cong_PC==i-1 & prev_sum_con_PC==e-1);
            end
        end
        % mincong = min(hstcong,[],'all');
        % % define indexes based on CNG and PC
        % for cur_cng = 1:4
        %     for pre_cng = 1:4
        %         cng_cp_idx(cur_cng,pre_cng,:) = randsample(find(cngidx(:,cur_cng)==1 & p_cngidx(:,pre_cng)==1 & cued_PC_matrix(:,s) & ~exclude_trials), mincong); % give indices over full array of all trials
        %     end
        % end
        
        % selecting the min across the 16 combinations often gives 0.
        % Therefore, we instead balance previous and post by selecting the
        % min across the sums within each quadrant of the 4x4 matrix
        q1 = sum(hstcong([1,2,5,6]));
        q2 = sum(hstcong([3,4,7,8]));
        q3 = sum(hstcong([9,10,13,14]));
        q4 = sum(hstcong([11,12,15,16]));
        mincong = min([q1,q2,q3,q4]);
        
        % define indexes based on CNG and PC
        cng_cp_idx = zeros(2,2,mincong);
        for cur_cng = 1:2
            for pre_cng = 1:2
                q_x_index = (cur_cng*(cur_cng==1) + (cur_cng+1)*(cur_cng>1)); % 1 or 3
                q_y_index = (pre_cng*(pre_cng==1) + (pre_cng+1)*(pre_cng>1)); % 1 or 3
                cng_cp_idx(cur_cng,pre_cng,:) = randsample(find(...
                                            (cngidx(:,q_x_index)==1 | cngidx(:,q_x_index+1)==1 ) & ...
                                            (p_cngidx(:,q_y_index)==1 | p_cngidx(:,q_y_index+1)==1 ) & ...
                                            cued_PC_matrix(:,s) & ~exclude_trials), ...
                                            mincong); % give indices over full array of all trials
            end
        end
    
        decimated_indices_bst =  reshape(cng_cp_idx,1,[])'; % will be less than original trial # - some are dropped to match congruency
        di(s).bst = cat(1,di(s).bst,decimated_indices_bst); % catenate the new bst iteration
        
     end

end

% mean(sum_cong(di(1).bst)) % must equal 1.5 (mean(0:3))
% mean(sum_cong(di(2).bst)) % must equal 1.5 (mean(0:3))
% mean(sum_cong(di(3).bst)) % must equal 1.5 (mean(0:3))
% mean(sum_cong(di(4).bst)) % must equal 1.5 (mean(0:3))

end
