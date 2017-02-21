clear all;
close all;
clc;

len = @(y) 1:length(y);

%% Main simulation loop
simType = 0;    %0 for single KF run
                %1 for Monte Carlo simulation
                %2 for several MC simulations (used for colour maps)

simLen = 500;
mcRuns = 1;

adaSwitch = 1;
alphaStep = 7;
sampleSize = 5;

step = 1;

x0 = [0;0];

%Signal Model Settings
sig.Qpre = 0.1;
sig.Qpost = 0.1;

sig.Rpre = 1;
sig.Rpost = 16;

T = 10000;

%sig.phi = [1 0; 0 1-step/T];
sig.phi = [1 0; 0 1-step/T];
sig.H = [1 1];

% KF Settings
kfMod = sig; % kf model
%kfMod.Rpre = 1;
%kfMod.Qpre = .5;
ada = tanisEst(1, alphaStep, sampleSize);

p0 = eye(2);

%%

switch simType
    case 0
        z = genSig(sig, simLen, x0);
        [xest, p, s, alpha, wincovar, inno] = kf(z, x0, p0, kfMod, ada, adaSwitch);
        [xestRef, pRef, sRef, alphaRef, wincovarRef, innoRef] = kf(z, x0, p0, kfMod, ada, 0);
        subplot(2,1,1)
        plot(len(xest), xest(1,:), len(xestRef), xestRef(1,:), len(z),z);
        
        subplot(2,1,2)
        plot(len(xest), z-xest(1,:), len(xestRef), z-xestRef(1,:), len(z),z);
        
    case 1
        [xEst, S, alpha, inno] = mc(mcRuns, simLen, sig, kfMod, ada, adaSwitch, x0, p0);
        %[xEstRef, Sref, alpharef, innoref] = mc(mcRuns, simLen, sig, kfMod, ada, 0, x0, p0);
        
    case 2
        alphaSteps = [1.1 1.25 1.5 1.75 2 2.5 3 4 5 7 10 15 25 50];
        sampleSizes = [2 3 4 5 7 10 15 20 25 40 50];
        
        xCol.mean = zeros(length(alphaSteps),length(sampleSizes));
        xCol.std = zeros(length(alphaSteps),length(sampleSizes));
        
        sCol.mean = zeros(length(alphaSteps),length(sampleSizes));
        sCol.std = zeros(length(alphaSteps),length(sampleSizes));
        
        alphaCol.mean = zeros(length(alphaSteps),length(sampleSizes));
        alphaCol.std = zeros(length(alphaSteps),length(sampleSizes));
        
        innoCol.sum = zeros(length(alphaSteps),length(sampleSizes));
        innoCol.std = zeros(length(alphaSteps),length(sampleSizes));
        
        rng(3141592)
        [xtemp, Stemp, alphatemp, innoTemp] = mc(mcRuns, simLen, sig, kfMod, ada, 0, x0, p0);
        
        ref.x.mean = mean(xtemp.val(:).^2);
        ref.x.std = std(xtemp.val(:).^2);
        
        ref.inno.sum = sum(innoTemp.val(:));
        
        for m = 1:length(alphaSteps)
            for n = 1:length(sampleSizes)
                [alphaSteps(m) sampleSizes(n)]
                rng(3141592)
                ada = tanisEst(1,alphaSteps(m), sampleSizes(n));
                
                [xEst, S, alpha] = mc(mcRuns, simLen, sig, kfMod, ada, adaSwitch, x0, p0);
                
                xCol.mean(m,n) = mean(xEst.val(:).^2) / ref.x.mean;
                xCol.std(m,n) = std(xEst.val(:).^2) / ref.x.std;
                
                sCol.mean(m,n) = mean((S.S.val(:) - S.wincovar(:)).^2);
                sCol.std(m,n) = std((S.S.val(:) - S.wincovar(:)).^2);
                
                alphaCol.mean(m,n) = alpha.mean;
                alphaCol.std(m,n) = alpha.std;
            end
        end
        
    otherwise
        
end
