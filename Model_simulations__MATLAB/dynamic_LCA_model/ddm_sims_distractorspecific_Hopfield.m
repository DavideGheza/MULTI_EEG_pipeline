function output = ddm_sims_distractorspecific_Hopfield(nrtrials,bias_weight,t2hw)
% Legacy DDM-like implementation from the NHB 2025 model, retained for
% reference and comparison with the dynamic LCA model (dynamic_LCA_sims.m). 
% Uses dimension-specific Hopfield-energy conflict signals and control updates.
% Example: output = ddm_sims_distractorspecific_Hopfield(160000,-5,4)
% Author: Davide Gheza and Wouter Kool
nrtasks = 4;

output = generate_task(nrtasks,nrtrials);
model = init_model_Hopfield(nrtasks,bias_weight,t2hw);

alpha = 0.02; % evidence_based: set to 0.1;
sigma = 0.1;

lambda = 0;
scaling_constant = 3.5;
beta = 0;

evidence = zeros(nrtrials,1);
rt = zeros(nrtrials,1);
response = zeros(nrtrials,1);
accuracy = zeros(nrtrials,1);

output.conflict = zeros(nrtasks,nrtrials);
output.control = zeros(nrtasks,nrtrials);

counter = 0;
for i = 1:nrtrials
    
    if mod(i,nrtrials/10) == 0
        counter = counter + 1;
        disp(counter)
    end    
    
    net_stim = model.stim2hidden_weights*output.stimulus_inputs(:,i);
    
    if i > 1  && output.task(i)==output.task(i-1) % only adapt on stay trials
       
        adaptation = output.control(:,i);
        adaptation(output.task(i)) = 0;
        
        net_task = model.task2hidden_weights*(output.task_inputs(:,i)-adaptation);
        
    else
        net_task = model.task2hidden_weights*output.task_inputs(:,i);
    end
    
    net_hidden = net_stim+net_task+model.hiddenbiases_weights;
    
    hiddenActs = logistic(net_hidden);
    
    net_output = model.hidden2output_weights*hiddenActs;
    
    outputActs = logistic(net_output);
    
    % DDM
    
    evidence(i,1) = outputActs(1) - outputActs(2);
    
    t = 0;
    V = 0;
    
    while abs(V) < 1
        t = t + 1;
        V = V + alpha*evidence(i,1)+normrnd(0,sigma);
    end
    
    rt(i,1) = t;
    response(i,1) = sign(V);
    
    if output.correct_response(i) == 1
        accuracy(i,1) = response(i,1)==1;
    else
        accuracy(i,1) = response(i,1)==-1;
    end
    
    % conflict detection                                                   % assuming energy within subset of pathway connections (hidden neurons - response units)
    partial_e = (hiddenActs .* outputActs') .* model.hidden2output_weights'; 
    for j = 1:nrtasks
        indices = 2*j-1:2*j;
        path_energy(j) = -sum(partial_e(indices,:),'all'); % sum(partial_e) is high when path is congruent with response. Is low and negative when incongruent. Path energy is the negative of this sum.
    end
    
    output.conflict(:,i) = path_energy;

    if i > 1
        output.control(:,i+1) = lambda*output.control(:,i) + (1-lambda)*(scaling_constant*output.conflict(:,i)+beta);
    end
        
end

output.evidence = evidence;
output.rt = rt;
output.response = response;
output.accuracy = accuracy;

output.outputActs = outputActs;
output.evidence = evidence;

% analysis

for c = 1:nrtasks
    output.acc_sum_con(c,1) =  mean(output.accuracy(output.sum_congruency==c-1));
    output.rt_sum_con(c,1) =  mean(output.rt(output.sum_congruency==c-1 & output.accuracy == 1));
end
% CSE spaghetti
for i_current = 0:nrtasks-1
    for i_previous = 0:nrtasks-1
        output.prev_sum_rt(i_previous+1,i_current+1) = mean(output.rt(output.accuracy==1&output.task_count>1&output.previous_sum_congruency==i_previous&output.sum_congruency==i_current));
        output.prev_sum_acc(i_previous+1,i_current+1) = mean(output.accuracy(output.task_count>1&output.previous_sum_congruency==i_previous&output.sum_congruency==i_current));
    end
end
% CSE heatmap
for t_previous = 1:nrtasks
    for t_current = 1:nrtasks
        
        output.previnc_task_con_rt(t_previous,t_current) = mean(output.rt(output.accuracy==1&output.task_count>1&output.prev_congruency(:,t_previous)==0&output.congruency(:,t_current)==1));
        output.previnc_task_inc_rt(t_previous,t_current) = mean(output.rt(output.accuracy==1&output.task_count>1&output.prev_congruency(:,t_previous)==0&output.congruency(:,t_current)==0));
        
        output.prevcon_task_con_rt(t_previous,t_current) = mean(output.rt(output.accuracy==1&output.task_count>1&output.prev_congruency(:,t_previous)==1&output.congruency(:,t_current)==1));
        output.prevcon_task_inc_rt(t_previous,t_current) = mean(output.rt(output.accuracy==1&output.task_count>1&output.prev_congruency(:,t_previous)==1&output.congruency(:,t_current)==0));
        
        output.previnc_task_int_rt(t_previous,t_current) = mean(output.rt(output.accuracy==1&output.task_count>1&output.prev_congruency(:,t_previous)==0&output.congruency(:,t_current)==0))-mean(output.rt(output.accuracy==1&output.task_count>1&output.prev_congruency(:,t_previous)==0&output.congruency(:,t_current)==1));
        output.prevcon_task_int_rt(t_previous,t_current) = mean(output.rt(output.accuracy==1&output.task_count>1&output.prev_congruency(:,t_previous)==1&output.congruency(:,t_current)==0))-mean(output.rt(output.accuracy==1&output.task_count>1&output.prev_congruency(:,t_previous)==1&output.congruency(:,t_current)==1));
        
    end
end

% CSE spaghetti 2-lag
for i_current = 0:nrtasks-1
    for i_previous = 0:nrtasks-1
        output.prev_2_sum_rt(i_previous+1,i_current+1) = mean(output.rt(output.accuracy==1&output.task_count>1&output.previous_2_sum_congruency==i_previous&output.sum_congruency==i_current));
        output.prev_2_sum_acc(i_previous+1,i_current+1) = mean(output.accuracy(output.task_count>1&output.previous_2_sum_congruency==i_previous&output.sum_congruency==i_current));
    end
end
% CSE spaghetti 3-lag
for i_current = 0:nrtasks-1
    for i_previous = 0:nrtasks-1
        output.prev_3_sum_rt(i_previous+1,i_current+1) = mean(output.rt(output.accuracy==1&output.task_count>1&output.previous_3_sum_congruency==i_previous&output.sum_congruency==i_current));
        output.prev_3_sum_acc(i_previous+1,i_current+1) = mean(output.accuracy(output.task_count>1&output.previous_3_sum_congruency==i_previous&output.sum_congruency==i_current));
    end
end

end
