function output = generate_task(nrtasks,nrtrials)
% Generate MULTI task sequences, task inputs, and binary feature inputs.
% Author: Davide Gheza and Wouter Kool

available_tasks = 1:nrtasks;

switch_trial = 1;
counter = 0;

stimulus_inputs = nan(nrtasks*2,nrtrials);
task_inputs = nan(nrtasks,nrtrials);

for t = 1:nrtrials
    counter = counter + 1;
    output.task_count(t,1) = counter;
    
    if switch_trial
        task(t) = randsample(available_tasks,1);
        if t > 1
            while task(t)==task(t-1)
                task(t) = randsample(available_tasks,1);
            end
        end
        available_tasks(available_tasks==task(t)) = [];
        if isempty(available_tasks)
            available_tasks = 1:nrtasks;
        end
        switch_trial = 0;
    else
        
        task(t) = task(t-1);
        
        p = 0;
        if counter == 3
            p = 1/3;
        end
        if counter == 4
            p = 1/2;
        end
        if counter == 5
            p = 1;
        end
        
        if rand < p
            switch_trial = 1;
            counter = 0;
        end
        
    end
    
    for i = 1:nrtasks
        
        
        stimulus_inputs((1:2)+(i-1)*2,t) = randperm(2)-1;
        
        if i == task(t)
            task_inputs(i,t) = 1;
            correct_response(t,1) = find(stimulus_inputs((1:2)+(i-1)*2,t));
        else
            task_inputs(i,t) = 0;
        end
        
    end
    
    for i = 1:nrtasks
        
        if i == task(t)
            congruency(i,t) = nan;
        else
            congruency(i,t) = mean(stimulus_inputs((1:2)+(i-1)*2,t) == stimulus_inputs((1:2)+(task(t)-1)*2,t));
        end
        
    end
    
end

output.task = task';
output.stimulus_inputs = stimulus_inputs;
output.task_inputs = task_inputs;
output.correct_response = correct_response;
output.congruency = congruency';
output.sum_congruency = nansum(output.congruency,2);
output.prev_congruency = [nan(1,nrtasks); output.congruency(1:(end-1),:)];
output.previous_sum_congruency = [nan;output.sum_congruency(1:(end-1))];
output.previous_2_sum_congruency = [nan;output.previous_sum_congruency(1:(end-1))];
output.previous_3_sum_congruency = [nan;output.previous_2_sum_congruency(1:(end-1))];

end
