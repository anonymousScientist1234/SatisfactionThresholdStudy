function [taskSwitches, stateSwitches, loading, DOL, timeToEq] = makeMeasurements(stateMat, sMat, T, N, lambda, K, timesteps, stimType) 

[a, b] = size(stateMat); 

%calculate the number of state switches, which will be the number of times
%the current state is different from the previous state, summed across all
%individuals in the colony 
stateSwitches = 0; 
for i = 1:b
    states = stateMat(:, i); 
    for j = 2:a
        if states(j) ~= states(j-1)
            stateSwitches = stateSwitches+1; 
        end
    end
end

%calculate number of task switches, which ignores resting state
taskSwitches = 0; 
for i = 1:b
    states = stateMat(:, i);
    states = states(states>0);
    for j = 2:length(states)
        if states(j) ~= states(j-1)
            taskSwitches = taskSwitches+1; 
        end
    end
end

%for loading, we calculate n.. (optimal n), equation S.11
optimalN = (N+N*lambda*timesteps)/K; 
%loading is equation S.14
loading = sum(sum(stateMat>0))-optimalN; 

%here we calculate division of labor according to the procedure outlined in
%Gorelick et al., 2004
for i = 1:N
    behaviorVec = stateMat(:,i); 
    for j = 1:T
        task_mat(i,j) = sum(behaviorVec==j); 
    end
end

task_mat=task_mat'; 

probability_mat = task_mat ./ (sum(sum(task_mat)));
p_ind = repelem(1/N, N);
p_task = sum(probability_mat')./sum(sum(probability_mat)); 
h_ind = (-1)*sum(p_ind.*log(p_ind));
h_task = (-1)*sum(p_task.*log(p_task));

for q = 1:T
  for w = 1:N
    I_mat(q, w) = probability_mat(q, w) * (log(probability_mat(q, w) / (p_ind(w)*p_task(q))));
  end
end

I_mat(isnan(I_mat))=0;
I = sum(sum(I_mat));
%dol_task_into_ind = I/h_ind;
%dol_ind_into_task = I/h_task;
%symmetricDOL = I/sqrt(h_ind*h_task);
DOL = I/sqrt(h_ind*h_task);

%finally, we find the point at which the cues reach equilibrium within the
%colony. This will be the first timestep when the average cue across tasks 
% within the precision limit of its final value (averaged across the last 100 timesteps 
precision = .005; 

if T == 1
    averageCue = sMat';
else
    averageCue = mean(sMat'); 
end

meanS = mean(averageCue(timesteps-100:end));
index = find(averageCue>(meanS-precision) & averageCue < (meanS + precision)); 
timeToEq = min(index); 
