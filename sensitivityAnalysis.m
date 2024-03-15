clear

data = readtable('thresholdComparisons.csv');
thresholds = ["Response Threshold", "Satisfaction Threshold", "Composite Threshold"];
cues = ["Task Demand", "Task Completion"]; 
counter = 0; 

%for every threshold and cue type, we calculate the partial correlations of
%the 6 free parameters on each of the 3 performance metrics 
for i = 1:3
    for j = 1:2
        threshTemp = thresholds(i);
        cueTemp = cues(j); 
        dataTemp = data(data.ThresholdType == threshTemp & data.CueType ==cueTemp ,:);

        for k = 1:3
            Response = k+12; %this is done to ignore irrelevant columns in the dataset
            dataCorr = [table2array(dataTemp(:, Response)) dataTemp.N dataTemp.K dataTemp.P dataTemp.Beta dataTemp.Mu dataTemp.Sigma];
            [rho, pVal] = partialcorr(dataCorr);

            counter = counter+1; 
            rhoMat(counter, :) = rho(1, 2:7); 
            pMat(counter, :) = pVal(1, 2:7); 
    

        end

    end
end

%all of the partial correlations for table S1
rhoMat = (round(rhoMat, 4));

%these are the bolded terms for table S1
pMat<.05
