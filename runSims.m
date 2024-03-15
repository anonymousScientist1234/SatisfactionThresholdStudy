function [stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
    thresholds2, M, threshType, stimType, T, N, s, l) 

%initialize initial demand (or completion) of simulation. This is the
%initial 'shock' of stimulus
if stimType == "Task Demand"
    sMat = zeros(M+1, T); 
    sMat(1,:) = 1; 
else
    sMat = zeros(M+1, T)+1; 
    sMat(1,:) = 0; 
end

%this matrix contains the states of all workers for all timesteps. They are
%initilized to 0, as all ants start in the resting state 
stateMat = zeros(M+1, N); 

%iterate through timesteps
for t = 1:M

    %iterate through ants 
    for n = 1:N

        %calculate bias for each task
        if stimType == "Task Demand"
            bias = sMat(t,:)./sum(sMat(t,:)); 
        else
            bias = (1-sMat(t,:))./sum(1-sMat(t,:)); 
        end

        %if bias is not a number, that just means there is one task, so the
        %bias for that task relative to the zero other tasks is one 
        if isnan(bias)
            bias = 1; 
        end

        %calculate the probability that an ant will start and stop a task
        %given which decision rule is being used and given the current value of the cue (Table 1) 
        if threshType == "Response Threshold" && stimType == "Task Demand"
            pStart = bias.*(sMat(t,:).^beta)./(sMat(t,:).^beta+thresholds(:,n)'.^beta);
            pStop = repmat(P, 1, T); 
        elseif threshType == "Response Threshold" && stimType == "Task Completion"
            pStart = bias.*(thresholds(:,n)'.^beta)./(sMat(t,:).^beta+thresholds(:,n)'.^beta);
            pStop = repmat(P, 1, T); 
        elseif threshType == "Satisfaction Threshold" && stimType == "Task Demand"
            pStart = (1-P).*bias; 
            pStop = (thresholds(:,n)'.^beta)./(sMat(t,:).^beta+thresholds(:,n)'.^beta);
        elseif threshType == "Satisfaction Threshold" && stimType == "Task Completion"
            pStart = (1-P).*bias; 
            pStop = (sMat(t,:).^beta)./(sMat(t,:).^beta+thresholds(:,n)'.^beta);
        elseif threshType == "Composite Threshold" && stimType == "Task Demand"
            pStart = bias.*(sMat(t,:).^beta)./(sMat(t,:).^beta+thresholds(:,n)'.^beta);
            pStop = (thresholds2(:,n)'.^beta)./(sMat(t,:).^beta+thresholds2(:,n)'.^beta);
        elseif threshType == "Composite Threshold" && stimType == "Task Completion" 
            pStart = bias.*(thresholds(:,n)'.^beta)./(sMat(t,:).^beta+thresholds(:,n)'.^beta);
            pStop = (sMat(t,:).^beta)./(sMat(t,:).^beta+thresholds2(:,n)'.^beta);
        else
            pStart = repmat(s/T, 1, T); 
            pStop = repmat(l, 1, T); 
        end

        %if the ant is currently resting, we simulate the choice of either
        %continuing to rest or to start any task j
        if stateMat(t, n) == 0
         
           stateMat(t+1, n) = T-sum(rand<cumsum(pStart))+1; 

           if stateMat(t+1, n) == T+1
                stateMat(t+1, n) = 0; 
           end

        %if the ant is not resting, we determine whether or not she stops
        %her current task
        else 

            if rand < pStop(stateMat(t, n))
                stateMat(t+1, n) = 0; 
            else
                stateMat(t+1, n) = stateMat(t, n); 
            end

        end

    end

    %to calculate hoe much the cue changes this timestep, we find how many
    %workers are working on each; task for this timestep 
    states = stateMat(t+1,:);

    %we update cues according to equations 9 and 10
    if stimType == "Task Demand"
        for u = 1:T
            sMat(t+1, u) = sMat(t, u) + lambda - T*(K*sum(states==u))/N;
            if sMat(t+1, u) < 0
                sMat(t+1, u) = 0;
            end
            if sMat(t+1, u) > 1
                sMat(t+1, u) = 1; 
            end
        end
    else 
        for u = 1:T
            sMat(t+1, u) = sMat(t, u) - lambda + T*(K*sum(states==u))/N;
            if sMat(t+1, u) < 0
                sMat(t+1, u) = 0;
            end
            if sMat(t+1, u) > 1
                sMat(t+1, u) = 1; 
            end

        end

    end

end