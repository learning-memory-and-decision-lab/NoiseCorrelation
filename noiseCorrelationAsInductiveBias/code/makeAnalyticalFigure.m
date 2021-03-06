% Code for figure 4 of:
% Noise Correlations for Faster and More Robust Learning
% Matthew R. Nassar, Daniel Scott, Apoorva Bhandari
% Journal of Neuroscience 4 August 2021, 41 (31) 6740-6752; DOI: 10.1523/JNEUROSCI.3045-20.2021 

% Code for analytical simulations of noise correlation learning advantage


% Instructions: Navigate to noiseCorrelationAsInductiveBias directory
% run script. Note neural populations are smaller than paper for faster
% simulation. 
clear

addpath(genpath('.'))
%% Get vertical and horizontal gradients: 


numReps               = 20;
testTrials            = 81:100;
%input.targetSEM       = 10;    % standard error on mean population representation
input.sigToNoise      = 2;
input.nPools          = 2;
input.nNeuronsPerPool = 100;
input.targFR          = 1;
input.nonTargFR       = -1;
input.nTrials         = 100;
input.LR              = .0001;
input.doRL            = true;
input.invT            = 10000;
input.weightScale     = 1;
input.useRPE          = false; % this produces weird results -- leave off for now.
input.doOpt           = true
input.all_fracCorrNoise=0:.01:.2;


nCorrs=length(input.all_fracCorrNoise);

allOptAcc=[]; allAcc=[]; allFinalWeights={};
   

for i = 1:numReps
disp(sprintf('running rep %g', i));
output1=runCorrLearningSim_noNorm_dan(input);
allWsn(:,:, i)= output1.wsn;
allWpn(:,:, i)= output1.wpn;
end    
    
% close all
% figure; hold on
 exNC=2 % choose noise correlations
% % Analytical weight components
% plot(input.LR*0.5*(1:input.nTrials)*norm(output1.Mu), '-', 'Color', 'b');
% plot(input.LR*0.5*sqrt(1:input.nTrials)*sqrt(trace(output1.allCovMat(:,:,exNC))-100), '-', 'Color', 'r');
% % Observed weight components
% hold on
% plot(mean(allWsn(exNC,:,:), 3), 'o', 'Color', 'b');
% plot(mean(allWpn(exNC,:,:),3), 'o', 'Color', 'r');
% title('Weight Components')
% xlabel('Trial')
% ylabel('Component Norm')
% ff=legend({'Ex. simulated \Delta w_s','Ex. simulated \Delta w_{\perp}', 'Analytic \Delta w_s','Analytic \Delta w_{\perp}'}, 'Location', 'NorthWest')
% set(ff, 'box', 'off')



%% Figure out how much of an advantage you get for noise correlations at
% different levels of signal and different numbers of neurons. 
allMu=.2:.2:5;

allN =2.^(2:13);

for i = 1:length(allMu)
    for j = 1:length(allN)
        input.mu = allMu(i);
        input.n = allN(j);
        [noiseCorrAdvantage(i,j), SNR(i,j)] = get_noiseCorrAdvantage(input);
    end
end

output=get_learning_speed
output3 = get_learning_speed(10000)
 
 
 
addpath(genpath('~/Dropbox/sharedMatlabUtilities'))

makeMathFig