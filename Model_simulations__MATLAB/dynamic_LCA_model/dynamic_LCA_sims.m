function output = dynamic_LCA_sims(nrtrials,bias_weight,t2hw)
% Simulate the dynamic dimension-specific-timescale accumulator model.
% Replaces the earlier DDM-like decision stage with two response accumulators
% and integrates hidden-layer evidence with dimension-specific time constants.
% Example: output = dynamic_LCA_sims(160000,-5,4)
% Author: Davide Gheza
nrtasks = 4;

output = generate_task(nrtasks,nrtrials);
model = init_model_Hopfield(nrtasks,bias_weight,t2hw);

tau = 0.4; % response layer
tau_hidden = [0.1,0.1, 0.025,0.025, 0.015,0.015, 0.005,0.005]'; % dimension-specific rate constant (hidden layer). Weights past activations. Affects the running average of activation over time (large = fast update)
thres = 1; % response threshold

alpha = 0.05; % scale the response accumulator update (large = faster; together with sigma, it sets the SNR)
% Note: alpha scales the response accumulator update --> larger values lead to faster integration at response level, and to increase of interference of a fast dimension (which develops both an early and STRONG interference). Lower this value and an early interference is weaker than later ones.
sigma = 0.15;

lambda = 0;
scaling_constant = 0.2; % for AUC energy conflict detection
% scaling_constant = 3.5; % for peak-energy conflict detection (alternative metric)
beta = 0;

evidence = zeros(nrtrials,1);
rt = zeros(nrtrials,1);
response = zeros(nrtrials,1);
accuracy = zeros(nrtrials,1);

output.conflict = zeros(nrtasks,nrtrials);
output.conflict_peak = zeros(nrtasks,nrtrials); % for testing correlation - for plotting only
output.control = zeros(nrtasks,nrtrials);
output.path_energy = nan(nrtasks,nrtrials,100); % allocate 1000 time points based on max cycles expected
output.hidden_evidence = nan(nrtasks,nrtrials,100); % allocate 1000 time points based on max cycles expected - for plotting only (reduced to 100 for RAM)

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
    
    
    t = 0;
    outputACCs = zeros(2,1); % initialize response accumulators
    net_hidden_timeunit = net_stim+net_task+model.hiddenbiases_weights; % net input at (any) time t. Time-increment from input is constant
    net_hidden = zeros(8,1); % initialize net-hidden at t=1
    net_output = zeros(2,1); % initialize response units at t=1
    path_energy = zeros(4,1);% initialize path_energy at t=1 
    hidden_evidence = zeros(4,1);            % for plotting only

    while all(outputACCs(:,t+1) < thres)    % check if either accumulator crosses threshold
        t = t + 1;
        net_hidden(:,t+1) = tau_hidden.*net_hidden_timeunit + (1-tau_hidden).*net_hidden(:,t);
        
        hiddenActs = logistic(net_hidden(:,t+1));
        
        net_output(:,t+1) = tau*model.hidden2output_weights*hiddenActs + (1-tau)*net_output(:,t);
        
        outputActs = logistic(net_output(:,t+1)); % to force positive both output units
        
        % response accumulators
        delta_output = diff(outputActs);
        mu = alpha*delta_output;
        acc_update = normrnd(mu,sigma);

        outputACCs(:,t+1) = outputACCs(:,t) + [-acc_update; acc_update];

        % dmension-specific conflict detection: measures energy within subset of pathway connections (hidden neurons - response units)
        partial_e = (hiddenActs .* outputActs') .* model.hidden2output_weights'; 
        for j = 1:nrtasks
            indices = 2*j-1:2*j;
            path_energy(j,t+1) = -sum(partial_e(indices,:),'all'); % sum(partial_e) is high when path is congruent with response. Is low and negative when incongruent. Path energy is the negative of this sum.
        end
    
        % for plotting only
        hidden_contrasts = diff(hiddenActs);
        hidden_evidence_update = hidden_contrasts(1:2:end); % 
        hidden_evidence(:,t+1) = hidden_evidence(:,t) + hidden_evidence_update; % dimension-specific evidence (r - l), rectified to reflect congruency
    end
    

    rt(i,1) = t+1;
    response(i,1) = find(outputACCs(:,t+1)>=thres);

    if output.correct_response(i) == 1
        accuracy(i,1) = response(i,1)==1;
    else
        accuracy(i,1) = response(i,1)==2;
    end
    

    % plot(path_energy','LineWidth',2)
    % xlabel('Time (cycles)')
    % ylabel('Energy')
    % legend
    
    output.conflict(:,i) = trapz(path_energy'); % AUC energy, integrated across all processing cycles (t) for a trial (i)
    output.conflict_peak(:,i) = max(abs(path_energy),[],2) .* sign(mean(path_energy,2)); % (NOT USED) Alternative metric: peak energy across all processing cycles (t) for a trial (i). Warning: energy can flip sign over time! AUC is more robust than peak-detection

    if i > 1
        output.control(:,i+1) = lambda*output.control(:,i) + (1-lambda)*(scaling_constant*output.conflict(:,i)+beta);
    end
    % if i<101     % (reduced to 100 for RAM)
    %     output.path_energy(:,i,1:t+1) = path_energy;                    
    %     output.hidden_evidence(:,i,1:t+1) = hidden_evidence; % for plotting only
    % end
end

output.rt = rt;
output.response = response;
output.accuracy = accuracy;
output.outputActs = outputActs;




%% analysis

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
        
        
        output.previnc_task_int_acc(t_previous,t_current) = mean(output.accuracy(output.task_count>1&output.prev_congruency(:,t_previous)==0&output.congruency(:,t_current)==0))-mean(output.accuracy(output.task_count>1&output.prev_congruency(:,t_previous)==0&output.congruency(:,t_current)==1));
        output.prevcon_task_int_acc(t_previous,t_current) = mean(output.accuracy(output.task_count>1&output.prev_congruency(:,t_previous)==1&output.congruency(:,t_current)==0))-mean(output.accuracy(output.task_count>1&output.prev_congruency(:,t_previous)==1&output.congruency(:,t_current)==1));

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

% conflict dynamics. Plot energy by pathway congruency (trial average)
for j = 1:nrtasks
    for cong = 0:1
        output.path_en_con(:,j,cong+1) = mean(output.hidden_evidence(j,output.path_energy(:,j)==cong,:),2,'omitnan'); % omitnan: cycles=RT varies by trial. NAN were used as fillers.
        % output.path_en_con(:,j,cong+1) = mean(output.hidden_evidence(j,output.congruency(:,j)==cong,:),2,'omitnan'); % omitnan: cycles=RT varies by trial. NAN were used as fillers.
    end
end

% delta-plots. Plot interference effects by dimension and RT quantile
cuts = 19;
quants = nan(nrtasks,2,cuts);
for j = 1:nrtasks
    for cong = 0:1
        quants(j,cong+1,:) = quantile(output.rt(output.congruency(:,j)==cong),cuts);
        quantileAssignment = discretize(output.rt, [-Inf, squeeze(quants(j,cong+1,:))', Inf]);
        for bin = 1:cuts+1
            output.acc_bins(bin,j,cong+1) = mean(output.accuracy(output.congruency(:,j)==cong & quantileAssignment==bin));
        end
    end
end
output.int_bins = output.acc_bins(:,:,1) - output.acc_bins(:,:,2); % incongruent - congruent
%%%%%%%%
figure
fullmap = viridis(256);           
idx = round(linspace(1,256,4));  % pick evenly spaced points, for 4 lines
cmap = fullmap(idx,:);
hold on
for i = 1:4
    plot(output.int_bins(:,i)*100, 'Color', cmap(i,:), 'LineWidth', 6);
end
xlabel('RT quantile', 'FontSize', 22)
ylabel('Delta accuracy inc-con (%)', 'FontSize', 22)
set(gca, 'FontSize', 22);
% legend({'Color','Shape','Edge','Motion'})
hold off
% save out figure
aspectRatio = 1.5;
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
newWidth = screenposition(4) * aspectRatio;  % width and height based on aspect ratio
newHeight = screenposition(4);  % Keep the height the same
set(gcf, ...
    'PaperPosition', [0 0 newWidth newHeight], ...
    'PaperSize', [newWidth newHeight]);

filename = 'acc_sims_delta_plot_accuracy';
print('-dpdf', '-vector', filename);
%%%%%%%%

% testing correlation bw conflict detection methods (AUC vs. peak)
figure
plot(output.conflict(1,:), output.conflict_peak(1,:)) % e.g., dim 1

end

