function model = init_model_Hopfield(nrtasks,bias_weight,t2hw)
% Initialize task, stimulus, hidden-bias, and response-layer weights.
% Author: Davide Gheza

model.task2hidden_weights = zeros(nrtasks,nrtasks*2);

hidden2output_weights = [];

for t = 1:nrtasks
    model.task2hidden_weights(t,((1:2)+2*(t-1))) = t2hw; 
    
    hidden2output_weights = [hidden2output_weights; 6 -6; -6 6];
end

model.task2hidden_weights = model.task2hidden_weights';
model.stim2hidden_weights = eye(nrtasks*2);
model.hiddenbiases_weights = bias_weight*ones(nrtasks*2,1);
model.hidden2output_weights = hidden2output_weights';

% edit to simulate graded dimensions' processing strength 
% model.stim2hidden_weights(:,[end-5,end-4]) = model.stim2hidden_weights(:,[end-5,end-4]) * 0.95;
% model.stim2hidden_weights(:,[end-3,end-2]) = model.stim2hidden_weights(:,[end-3,end-2]) * 0.85;
% model.stim2hidden_weights(:,[end-1,end]) = model.stim2hidden_weights(:,[end-1,end]) * 0.5; % reduces weight on last dim. stim-hidden

end
