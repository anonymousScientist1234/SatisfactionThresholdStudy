clear

%% Set constants (including ranges for free parameters) and initialize vectors which will contain data 
sims = 100; 
Nmin = 10;
Nmax = 1000; 
betaMin = 1;
betaMax = 10;
lambda = .001;
muMin = 0; 
muMax = 1;
sigmaMin = 0;
sigmaMax = 1; 
pMin = .5;
pMax = 1; 
TValues = [1 2 4 8]; 
tag = "sim";
M = 5000;
kMin = .001;
kMax = .1;
tolerance = .005; 

taskSwitchVec = zeros(sims*4*2*4, 1);
loadingVec = zeros(sims*4*2*4, 1);
DOLVec = zeros(sims*4*2*4, 1);
timeToEqVec = zeros(sims*4*2*4, 1);

stimTypeVec = strings(sims*4*2*4, 1); 
threshTypeVec = strings(sims*4*2*4, 1); 
TVec = zeros(sims*4*2*4, 1);

NVec = zeros(sims*4*2*4, 1);
betaVec = zeros(sims*4*2*4, 1);
PVec = zeros(sims*4*2*4, 1);
KVec = zeros(sims*4*2*4, 1);
muVec = zeros(sims*4*2*4, 1);
sigmaVec = zeros(sims*4*2*4, 1);

sVec = zeros(sims*4*2*4, 1);
lVec = zeros(sims*4*2*4, 1);

%% Run sims

counter = 0; 

%iterate through each simulation
for i = 1:sims

    %iterate through task numbers 
    for j = 1:length(TValues)

        T = TValues(j);
        %choose parameter values for this simulation. as we later run
        %through all combinations of thresholds and cues, this amounts to a
        %full factorial design where every parameter combination is tested
        %within each threshold type 
        [N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
        muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

        %iterate through cue types
        for k = 1:2

            if k == 1
                stimType = "Task Demand";
            else
                stimType = "Task Completion";
            end

            %iterate through threshold types 
            for u = 1:4

                if u == 1
                    threshType = "Response Threshold";
                    s = NaN; 
                    l = NaN; 
                elseif u == 2
                    threshType = "Satisfaction Threshold";
                    s = NaN; 
                    l = NaN; 
                elseif u == 3
                    threshType = "Composite Threshold";
                    s = NaN; 
                    l = NaN; 
                else
                    threshType = "Random Choice";
                    ratio = (lambda+1/M)/K;
                    s = rand;
                    l = rand;
                    counter2 = 0; 
                    flag = 0;
                    %check if ratio for random choice is in the tolerance
                    %limit
                    while flag == 0
                        counter2 = counter2 + 1;
                        s = rand;
                        l = rand;
                        if s/(s+l)>(ratio-tolerance) && s/(s+l)<(ratio+tolerance) 
                            flag = 1;
                        end
                    
                        if counter2 > 2500
                            flag = 1;
                        end

                    end
                end

                counter = counter+1; 
    
                %run the simulation 
                [stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
                    thresholds2, M, threshType, stimType, T, N, s, l); 
    
                %measure the outcomes of the simulation 
                [taskSwitch, loading, DOL, timeToEq] = makeMeasurements(stateMat, sMat, T, N, lambda, K, M, stimType);  

                %collect data on performance and parameter values 
                taskSwitchVec(counter) = taskSwitch;
                loadingVec(counter) = loading;
                DOLVec(counter) = DOL;
                timeToEqVec(counter) = timeToEq; 
    
                stimTypeVec(counter) = stimType;
                threshTypeVec(counter) = threshType;
                TVec(counter) = T;

                NVec(counter) = N;
                betaVec(counter) = beta;
                PVec(counter) = P;
                KVec(counter) = K;
                muVec(counter) = mu;
                sigmaVec(counter) = sigma;

                sVec(counter) = s;
                lVec(counter) = l;

                counter/(sims*4*2*4)

            end

        end

    end

end

%output data as a csv file 
data = table(threshTypeVec, stimTypeVec, TVec, NVec, KVec, PVec, betaVec, muVec, sigmaVec, sVec, lVec, DOLVec, taskSwitchVec,loadingVec, timeToEqVec); 
data.Properties.VariableNames = ["ThresholdType", "CueType", "T", "N", "K", "P", "Beta", "Mu", "Sigma", "s", "l", "DOL", "TaskSwitches", "Loading", "TimetoEq"];
writetable(data,'thresholdComparisons.csv');

%system('shutdown -s')