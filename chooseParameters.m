function [N, beta, P, K, mu, sigma, thresholds, thresholds2] = chooseParameters(kMin, kMax,...
                    muMin, muMax, sigmaMin, sigmaMax, pMin, pMax, Nmin, Nmax, betaMin, betaMax, T, tag)

%Set constants within sim. Rules for sampling are set in table 2. 
%If/else statements are for the figs7andS2 script, as different plots
%require different and specific parameter values. Main sims just use the
%'else' condition at the end 

if tag == "p1"
    T = 1; 
    N = 1000;
    beta = 6; 
    P = .9; 
    K = .01; 
    mu = .5; 
    sigma = .5; 
    pd = makedist('Normal','mu',mu,'sigma',sigma); 
    pdt = truncate(pd,0,1);
    thresholds = random(pdt, T, N); 
    thresholds2 = random(pdt, T, N); 
elseif tag == "p2"
    T = 4; 
    N = 1000;
    beta = 6; 
    P =  .9; 
    K = .01; 
    mu = .5; 
    sigma = .5; 
    pd = makedist('Normal','mu',mu,'sigma',sigma); 
    pdt = truncate(pd,0,1);
    thresholds = random(pdt, T, N); 
    thresholds2 = random(pdt, T, N); 
elseif tag == "p3"
    T = 1; 
    N = 50;
    beta = 6; 
    P = .9; 
    K = .01; 
    mu = .5; 
    sigma = .5; 
    pd = makedist('Normal','mu',mu,'sigma',sigma); 
    pdt = truncate(pd,0,1);
    thresholds = random(pdt, T, N); 
    thresholds2 = random(pdt, T, N); 
elseif tag == "p4"
    T = 1; 
    N = 1000;
    beta = 6; 
    P = .9; 
    K = .0015; 
    mu = .5; 
    sigma = .5; 
    pd = makedist('Normal','mu',mu,'sigma',sigma); 
    pdt = truncate(pd,0,1);
    thresholds = random(pdt, T, N); 
    thresholds2 = random(pdt, T, N); 
else
    N = randi([Nmin Nmax],1, 1); 
    beta = unifrnd(betaMin,betaMax); 
    P = unifrnd(pMin,pMax); 
    %P = 1./unifrnd(1/pMax,1/pMin); 
    K = 1./unifrnd(1/kMax,1/kMin);
    mu = unifrnd(muMin,muMax); 
    sigma = unifrnd(sigmaMin,sigmaMax); 
    pd = makedist('Normal','mu',mu,'sigma',sigma); 
    pdt = truncate(pd,0,1);
    thresholds = random(pdt, T, N); 
    thresholds2 = random(pdt, T, N); 
end