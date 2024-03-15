clear

%% set consants 
sims = 100;
timesteps = 2500; 
Nmin = 10;
Nmax = 1000; 
betaMin = 1;
betaMax = 10; 
lambda = .001;
kMin = .001;
kMax = .05;
muMin = 0; 
muMax = 1;
sigmaMin = 0;
sigmaMax = 1; 
pMin = .85; 
pMax = 1;
TValues = [1 2 4 8]; 
%alpha1 = 1;
%alpha2 = .85;
%alpha3 = .5;
%alpha4 = .25;
alpha1 = 1;
alpha2 = 1;
alpha3 = 1;
alpha4 = 1;
T = 1; 

%% P1 
stimType = "Task Demand";
threshType = "Response Threshold";
tag = "p1"; 

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

tolerance = .005;
flag = 0;
ratio = (lambda+1/5000)/K;
s = rand;
l = rand;
counter2 = 0; 
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

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatRTp1 = sMat(:,1); 
actionsTensRTp1 = sum(stateMat'>0)./N; 

threshType = "Satisfaction Threshold";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatSTp1 = sMat(:,1); 
actionsTensSTp1 = sum(stateMat'>0)./N; 

threshType = "Composite Threshold";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatCTp1 = sMat(:,1); 
actionsTensCTp1 = sum(stateMat'>0)./N; 

threshType = "Random Choice";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N, s, l); 

oneSMatRCp1 = sMat(:,1); 
actionsTensRCp1 = sum(stateMat'>0)./N; 

%% P2 
threshType = "Response Threshold";
tag = "p2"; 
T = 4; 

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

tolerance = .005;
flag = 0;
ratio = (lambda+1/5000)/K;
s = rand;
l = rand;
counter2 = 0; 
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


[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatRTp2 = mean(sMat')'; 
actionsTensRTp2 = sum(stateMat'>0)./N; 

threshType = "Satisfaction Threshold";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatSTp2 = mean(sMat')'; 
actionsTensSTp2 = sum(stateMat'>0)./N; 

threshType = "Composite Threshold";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatCTp2 = mean(sMat')'; 
actionsTensCTp2 = sum(stateMat'>0)./N; 

threshType = "Random Choice";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N, s, l); 

oneSMatRCp2 = sMat(:,1); 
actionsTensRCp2 = sum(stateMat'>0)./N; 

%% P3 
threshType = "Response Threshold";
tag = "p3"; 
T = 1; 

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

tolerance = .005;
flag = 0;
ratio = (lambda+1/5000)/K;
s = rand;
l = rand;
counter2 = 0; 
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


[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatRTp3 = sMat(:,1); 
actionsTensRTp3 = sum(stateMat'>0)./N; 

threshType = "Satisfaction Threshold";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatSTp3 = sMat(:,1); 
actionsTensSTp3 = sum(stateMat'>0)./N; 

threshType = "Composite Threshold";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatCTp3 = sMat(:,1); 
actionsTensCTp3 = sum(stateMat'>0)./N; 

threshType = "Random Choice";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N, s, l); 

oneSMatRCp3 = sMat(:,1); 
actionsTensRCp3 = sum(stateMat'>0)./N; 

%% P4 
threshType = "Response Threshold";
tag = "p4"; 

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

tolerance = .005;
flag = 0;
ratio = (lambda+1/5000)/K;
s = rand;
l = rand;
counter2 = 0; 
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


[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatRTp4 = sMat(:,1); 
actionsTensRTp4 = sum(stateMat'>0)./N; 

threshType = "Satisfaction Threshold";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatSTp4 = sMat(:,1); 
actionsTensSTp4 = sum(stateMat'>0)./N; 

threshType = "Composite Threshold";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatCTp4 = sMat(:,1); 
actionsTensCTp4 = sum(stateMat'>0)./N; 

threshType = "Random Choice";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N, s, l); 

oneSMatRCp4 = sMat(:,1); 
actionsTensRCp4 = sum(stateMat'>0)./N; 

%% Plots

%colors in hex: https://www.mathworks.com/help/matlab/creating_plots/specify-plot-colors.html
x = 0:timesteps;

figure(1)
tiledlayout(2,4)

%cue 
nexttile

plot(x, oneSMatRTp1, "Color", [.4 .8 .93333 alpha1],'LineWidth',2)
hold on
plot(x, oneSMatSTp1, "Color", [0 .6 .53333 alpha2],'LineWidth',2)
plot(x, oneSMatCTp1, "Color", [.8 .2 .06666 alpha3],'LineWidth',2)
plot(x, oneSMatRCp1, "Color", [.73333 .73333 .73333 alpha4],'LineWidth',2)
xlim([0, timesteps+1])  
ylim([0, 1])
ttl = title('A');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
ttl.FontSize = 16; 

%legend('RT','ST', 'CT')
xlabel("Timestep")
ylabel("Task Demand Cue")

nexttile

plot(x, oneSMatRTp2, "Color", [.4 .8 .93333 alpha1],'LineWidth',2, 'LineStyle', '-')
hold on
plot(x, oneSMatSTp2, "Color", [0 .6 .53333 alpha2],'LineWidth',2, 'LineStyle', '--')
plot(x, oneSMatCTp2, "Color", [.8 .2 .06666 alpha3],'LineWidth',2, 'LineStyle', '-.')
plot(x, oneSMatRCp2, "Color", [.73333 .73333 .73333 alpha4],'LineWidth',2, 'LineStyle', ':')
xlim([0, timesteps+1])  
ylim([0, 1])
ttl = title('B');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
ttl.FontSize = 16; 
xlabel("Timestep")
ylabel("Task Demand Cue")

nexttile

plot(x, oneSMatRTp3, "Color", [.4 .8 .93333 alpha1],'LineWidth',2, 'LineStyle', '-')
hold on
plot(x, oneSMatSTp3, "Color", [0 .6 .53333 alpha2],'LineWidth',2, 'LineStyle', '--')
plot(x, oneSMatCTp3, "Color", [.8 .2 .06666 alpha3],'LineWidth',2, 'LineStyle', '-.')
plot(x, oneSMatRCp3, "Color", [.73333 .73333 .73333 alpha4],'LineWidth',2,'LineStyle',  ':')
xlim([0, timesteps+1])  
ylim([0, 1])
ttl = title('C');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
ttl.FontSize = 16; 
xlabel("Timestep")
ylabel("Task Demand Cue")

nexttile

plot(x, oneSMatRTp4, "Color", [.4 .8 .93333 alpha1],'LineWidth',2,'LineStyle',  '-')
hold on
plot(x, oneSMatSTp4, "Color", [0 .6 .53333 alpha2],'LineWidth',2, 'LineStyle', '--')
plot(x, oneSMatCTp4, "Color", [.8 .2 .06666 alpha3],'LineWidth',2, 'LineStyle', '-.')
plot(x, oneSMatRCp4, "Color", [.73333 .73333 .73333 alpha4],'LineWidth',2, 'LineStyle', ':')
xlim([0, timesteps+1])  
ylim([0, 1])
ttl = title('D');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
ttl.FontSize = 16; 
xlabel("Timestep")
ylabel("Task Demand Cue")
legend('Response Threshold','Satisfaction Threshold', 'Composite Threshold', 'Random Choice', 'Location','southeast')

%n 
nexttile

plot(x, actionsTensRTp1,"Color", [.4 .8 .93333 alpha1],'LineWidth',2, 'LineStyle', '-')
hold on
plot(x, actionsTensSTp1,"Color", [0 .6 .53333 alpha2],'LineWidth',2, 'LineStyle', '--')
plot(x, actionsTensCTp1,"Color", [.8 .2 .06666 alpha3],'LineWidth',2,'LineStyle',  '-.')
plot(x, actionsTensRCp1, "Color", [.73333 .73333 .73333 alpha4],'LineWidth', 2, 'LineStyle', ':')
xlim([0, timesteps+1])  
ylim([0, 1])
ttl = title('E');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
ttl.FontSize = 16; 

%legend('RT','ST', 'CT')
xlabel("Timestep")
ylabel("Prop. Active Workers")

nexttile

plot(x, actionsTensRTp2,"Color", [.4 .8 .93333 alpha1],'LineWidth',2,'LineStyle',  '-')
hold on
plot(x, actionsTensSTp2,"Color", [0 .6 .53333 alpha2],'LineWidth',2,'LineStyle',  '--')
plot(x, actionsTensCTp2,"Color", [.8 .2 .06666 alpha3],'LineWidth',2, 'LineStyle', '-.')
plot(x, actionsTensRCp2, "Color", [.73333 .73333 .73333 alpha4],'LineWidth',2, 'LineStyle', ':')
xlim([0, timesteps+1])  
ylim([0, 1])
ttl = title('F');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
ttl.FontSize = 16; 
xlabel("Timestep")
ylabel("Prop. Active Workers")

nexttile

plot(x, actionsTensRTp3,"Color", [.4 .8 .93333 alpha1],'LineWidth',2,'LineStyle',  '-')
hold on
plot(x, actionsTensSTp3,"Color", [0 .6 .53333 alpha2],'LineWidth',2, 'LineStyle', '--')
plot(x, actionsTensCTp3,"Color", [.8 .2 .06666 alpha3],'LineWidth',2, 'LineStyle', '-.')
plot(x, actionsTensRCp3, "Color", [.73333 .73333 .73333 alpha4],'LineWidth',2, 'LineStyle', ':')
xlim([0, timesteps+1])  
ylim([0, 1])
ttl = title('G');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
ttl.FontSize = 16; 
xlabel("Timestep")
ylabel("Prop. Active Workers")

nexttile

plot(x, actionsTensRTp4,"Color", [.4 .8 .93333 alpha1],'LineWidth',2,'LineStyle',  '-')
hold on
plot(x, actionsTensSTp4,"Color", [0 .6 .53333 alpha2],'LineWidth',2,'LineStyle',  '--')
plot(x, actionsTensCTp4,"Color", [.8 .2 .06666 alpha3],'LineWidth',2, 'LineStyle', '-.')
plot(x, actionsTensRCp4, "Color", [.73333 .73333 .73333 alpha4],'LineWidth',2, 'LineStyle', ':')
xlim([0, timesteps+1])  
ylim([0, 1])
ttl = title('H');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
ttl.FontSize = 16; 
xlabel("Timestep")
ylabel("Prop. Active Workers")

%% Task Completion 

%% P1 
stimType = "Task Completion";
threshType = "Response Threshold";
tag = "p1"; 

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

tolerance = .005;
flag = 0;
ratio = (lambda+1/5000)/K;
s = rand;
l = rand;
counter2 = 0; 
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


[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatRTp1 = sMat(:,1); 
actionsTensRTp1 = sum(stateMat'>0)./N; 

threshType = "Satisfaction Threshold";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatSTp1 = sMat(:,1); 
actionsTensSTp1 = sum(stateMat'>0)./N; 

threshType = "Composite Threshold";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatCTp1 = sMat(:,1); 
actionsTensCTp1 = sum(stateMat'>0)./N; 

threshType = "Random Choice";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N, s, l); 

oneSMatRCp1 = sMat(:,1); 
actionsTensRCp1 = sum(stateMat'>0)./N; 

%% P2 
threshType = "Response Threshold";
tag = "p2"; 
T = 4; 

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

tolerance = .005;
flag = 0;
ratio = (lambda+1/5000)/K;
s = rand;
l = rand;
counter2 = 0; 
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


[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatRTp2 = mean(sMat')'; 
actionsTensRTp2 = sum(stateMat'>0)./N; 

threshType = "Satisfaction Threshold";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatSTp2 = mean(sMat')'; 
actionsTensSTp2 = sum(stateMat'>0)./N; 

threshType = "Composite Threshold";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatCTp2 = mean(sMat')'; 
actionsTensCTp2 = sum(stateMat'>0)./N; 

threshType = "Random Choice";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N, s, l); 

oneSMatRCp2 = sMat(:,1); 
actionsTensRCp2 = sum(stateMat'>0)./N; 

%% P3 
threshType = "Response Threshold";
tag = "p3"; 
T = 1; 

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

tolerance = .005;
flag = 0;
ratio = (lambda+1/5000)/K;
s = rand;
l = rand;
counter2 = 0; 
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


[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatRTp3 = sMat(:,1); 
actionsTensRTp3 = sum(stateMat'>0)./N; 

threshType = "Satisfaction Threshold";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatSTp3 = sMat(:,1); 
actionsTensSTp3 = sum(stateMat'>0)./N; 

threshType = "Composite Threshold";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatCTp3 = sMat(:,1); 
actionsTensCTp3 = sum(stateMat'>0)./N; 

threshType = "Random Choice";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N, s, l); 

oneSMatRCp3 = sMat(:,1); 
actionsTensRCp3 = sum(stateMat'>0)./N; 

%% P4 
threshType = "Response Threshold";
tag = "p4"; 

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

tolerance = .005;
flag = 0;
ratio = (lambda+1/5000)/K;
s = rand;
l = rand;
counter2 = 0; 
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

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatRTp4 = sMat(:,1); 
actionsTensRTp4 = sum(stateMat'>0)./N; 

threshType = "Satisfaction Threshold";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatSTp4 = sMat(:,1); 
actionsTensSTp4 = sum(stateMat'>0)./N; 

threshType = "Composite Threshold";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N); 

oneSMatCTp4 = sMat(:,1); 
actionsTensCTp4 = sum(stateMat'>0)./N; 

threshType = "Random Choice";

[N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag);

[stateMat, sMat] = runSims(P, K, beta, lambda, thresholds, ...
thresholds2, timesteps, threshType, stimType, T, N, s, l); 

oneSMatRCp4 = sMat(:,1); 
actionsTensRCp4 = sum(stateMat'>0)./N; 

%% Plots

x = 0:timesteps;

figure(2)
tiledlayout(2,4)

%cue 
nexttile

plot(x, oneSMatRTp1, "Color", [.4 .8 .93333 alpha1],'LineWidth',2, 'LineStyle', '-')
hold on
plot(x, oneSMatSTp1, "Color", [0 .6 .53333 alpha2],'LineWidth',2, 'LineStyle', '--')
plot(x, oneSMatCTp1, "Color", [.8 .2 .06666 alpha3],'LineWidth',2, 'LineStyle', '-.')
plot(x, oneSMatRCp1, "Color", [.73333 .73333 .73333 alpha4],'LineWidth',2, 'LineStyle', ':')
xlim([0, timesteps+1])  
ylim([0, 1])
ttl = title('A');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
ttl.FontSize = 16; 

%legend('RT','ST', 'CT')
xlabel("Timestep")
ylabel("Task Completion Cue")

nexttile

plot(x, oneSMatRTp2, "Color", [.4 .8 .93333 alpha1],'LineWidth',2, 'LineStyle', '-')
hold on
plot(x, oneSMatSTp2, "Color", [0 .6 .53333 alpha2],'LineWidth',2, 'LineStyle', '--')
plot(x, oneSMatCTp2, "Color", [.8 .2 .06666 alpha3],'LineWidth',2, 'LineStyle', '-.')
plot(x, oneSMatRCp2, "Color", [.73333 .73333 .73333 alpha4],'LineWidth',2, 'LineStyle', ':')
xlim([0, timesteps+1])  
ylim([0, 1])
ttl = title('B');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
ttl.FontSize = 16; 
xlabel("Timestep")
ylabel("Task Completion Cue")

nexttile

plot(x, oneSMatRTp3, "Color", [.4 .8 .93333 alpha1],'LineWidth',2, 'LineStyle', '-')
hold on
plot(x, oneSMatSTp3, "Color", [0 .6 .53333 alpha2],'LineWidth',2, 'LineStyle', '--')
plot(x, oneSMatCTp3, "Color", [.8 .2 .06666 alpha3],'LineWidth',2, 'LineStyle', '-.')
plot(x, oneSMatRCp3, "Color", [.73333 .73333 .73333 alpha4],'LineWidth',2, 'LineStyle', ':')
xlim([0, timesteps+1])  
ylim([0, 1])
ttl = title('C');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
ttl.FontSize = 16; 
xlabel("Timestep")
ylabel("Task Completion Cue")

nexttile

plot(x, oneSMatRTp4, "Color", [.4 .8 .93333 alpha1],'LineWidth',2, 'LineStyle', '-')
hold on
plot(x, oneSMatSTp4, "Color", [0 .6 .53333 alpha2],'LineWidth',2, 'LineStyle', '--')
plot(x, oneSMatCTp4, "Color", [.8 .2 .06666 alpha3],'LineWidth',2, 'LineStyle', '-.')
plot(x, oneSMatRCp4, "Color", [.73333 .73333 .73333 alpha4],'LineWidth',2, 'LineStyle', ':')
xlim([0, timesteps+1])  
ylim([0, 1])
ttl = title('D');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
ttl.FontSize = 16; 
xlabel("Timestep")
ylabel("Task Completion Cue")
legend('Response Threshold','Satisfaction Threshold', 'Composite Threshold', 'Random Choice', 'Location','northeast')

%n 
nexttile

plot(x, actionsTensRTp1,"Color", [.4 .8 .93333 alpha1],'LineWidth',2, 'LineStyle', '-')
hold on
plot(x, actionsTensSTp1,"Color", [0 .6 .53333 alpha2],'LineWidth',2, 'LineStyle', '--')
plot(x, actionsTensCTp1,"Color", [.8 .2 .06666 alpha3],'LineWidth',2, 'LineStyle', '-.')
plot(x, actionsTensRCp1, "Color", [.73333 .73333 .73333 alpha4],'LineWidth',2, 'LineStyle', ':')
xlim([0, timesteps+1])  
ylim([0, 1])
ttl = title('E');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
ttl.FontSize = 16; 

%legend('RT','ST', 'CT')
xlabel("Timestep")
ylabel("Prop. Active Workers")

nexttile

plot(x, actionsTensRTp2,"Color", [.4 .8 .93333 alpha1],'LineWidth',2, 'LineStyle', '-')
hold on
plot(x, actionsTensSTp2,"Color", [0 .6 .53333 alpha2],'LineWidth',2, 'LineStyle', '--')
plot(x, actionsTensCTp2,"Color", [.8 .2 .06666 alpha3],'LineWidth',2, 'LineStyle', '-.')
plot(x, actionsTensRCp2, "Color", [.73333 .73333 .73333 alpha4],'LineWidth',2, 'LineStyle', ':')
xlim([0, timesteps+1])  
ylim([0, 1])
ttl = title('F');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
ttl.FontSize = 16; 
xlabel("Timestep")
ylabel("Prop. Active Workers")

nexttile

plot(x, actionsTensRTp3,"Color", [.4 .8 .93333 alpha1],'LineWidth',2, 'LineStyle', '-')
hold on
plot(x, actionsTensSTp3,"Color", [0 .6 .53333 alpha2],'LineWidth',2, 'LineStyle', '--')
plot(x, actionsTensCTp3,"Color", [.8 .2 .06666 alpha3],'LineWidth',2, 'LineStyle', '-.')
plot(x, actionsTensRCp3, "Color", [.73333 .73333 .73333 alpha4],'LineWidth',2, 'LineStyle', ':')
xlim([0, timesteps+1])  
ylim([0, 1])
ttl = title('G');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
ttl.FontSize = 16; 
xlabel("Timestep")
ylabel("Prop. Active Workers")

nexttile

plot(x, actionsTensRTp4,"Color", [.4 .8 .93333 alpha1],'LineWidth',2, 'LineStyle', '-')
hold on
plot(x, actionsTensSTp4,"Color", [0 .6 .53333 alpha2],'LineWidth',2, 'LineStyle', '--')
plot(x, actionsTensCTp4,"Color", [.8 .2 .06666 alpha3],'LineWidth',2, 'LineStyle', '-.')
plot(x, actionsTensRCp4, "Color", [.73333 .73333 .73333 alpha4],'LineWidth',2, 'LineStyle', ':')
xlim([0, timesteps+1])  
ylim([0, 1])
ttl = title('H');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
ttl.FontSize = 16; 
xlabel("Timestep")
ylabel("Prop. Active Workers")
