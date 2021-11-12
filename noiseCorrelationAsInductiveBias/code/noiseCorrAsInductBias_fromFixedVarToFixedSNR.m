% NOISE CORRELATIONS AS INDUCTIVE BIASES for learning?


% This code was spawned from noiseCorrAsInductBias.m
% the goal here is to run 2 noise correlation conditions:
% 0, 0.2 
% within each condition, we will parameterize our assumption on fixed SNR
% So -- if assumpParam=1, we use our previous methods.
%       if assumpParam=0, we assume an equal variance across all NC
%                       conditions
%       for intermediate assumpParams, we set sigma equal to a weighted average of these methods



% The math describing how this will work is in a document called:
% multSignalAndNoiseMethods.docx

% cd ~/Dropbox (Brown)/noiseCorrelationAsInductiveBias/


%% RUN THIS SECTION FOR EITHER SIMULATION?

targetSEM       =  1;    % standard error on mean population representation
nPools          =  2;
nNeuronsPerPool = 100;
targFR          =  1; % Now these will be manipulated. 
nonTargFR       = -1; % 
nTrials         =100;
LR              =.0001;
doRL            =true;
invT            =10000;
weightScale     =1;
useRPE          =false; % online predictions
doOpt           =true;
initWtVar       =0.00001;
whichComp       =1;

% Sort out path:
switch whichComp
    case 1
        baseDir='~/Dropbox';
    case 2
        baseDir='~/Dropbox (Brown)' % Put your dropbox path here
end
addpath(genpath(fullfile(baseDir, 'sharedMatlabUtilities')));
cd(fullfile(baseDir, 'noiseCorrelationAsInductiveBias'));
addpath(genpath('./code'));


if ~doRL
    numOutNeurons = 1;
else
    numOutNeurons = nPools;
end

totVar = nNeuronsPerPool.*targetSEM.^2; % FIXED to a target level!
all_fracCorrNoise= 0:.1:1;
allStims=randi(nPools, nTrials, 1); % pick one set of stims for all mods

% Target values currently only make sense for 2 pools... 
targetOutputs=[100; -100;]; % choose values for decision neuron output

% create a list of neuron pool assignments:
poolID=[];
for i = 1:nPools
    poolID(end+1:end+nNeuronsPerPool,1)=ones(nNeuronsPerPool,1).*i;
end
meanFR=(targFR+nonTargFR.*(nPools-1))./nPools;

if ~doRL
optWeights=(poolID==1).*targetOutputs(1)./sum((poolID==1)) + ...
    (poolID==2).*targetOutputs(2)./sum((poolID==2));
else
    % Setup optimal weights here. 
    optWeights=[]
    for pp=1:max(poolID)
    optWeights(poolID==pp,1:(pp-1))=nonTargFR;
    optWeights(poolID==pp,pp)=targFR;
    optWeights(poolID==pp,(pp+1):max(poolID))=nonTargFR;
    end
    
    % normalize optimal weights using the normalization scheme implemented:
    
    optWeights(:,1)=  optWeights(:,1)./(sum(optWeights(:,1).^2)./weightScale);
    optWeights(:,2)=  optWeights(:,2)./(sum(optWeights(:,2).^2)./weightScale);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Network #1  -- two layers, simplest possible proof of concept:
%          Run multiple sims to get "averaged" performance:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE -- initial weights matter quite a bit for the "weight profile" plot
     

numReps               = 10000;
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
input.all_fracCorrNoise=[0, .2];
input.fixNoise        = true;
input.assumpParam     = 0:.1:1;


nCorrs=length(input.all_fracCorrNoise);





allOptAcc=[]; allAcc=[]; allFinalWeights={};
for i = 1:numReps
    disp(sprintf('running rep %g', i));
    output=runCorrLearningSim_altAssumptions(input);
    allOptAcc=cat(4, allOptAcc, output.optAccuracy);
    allAcc=cat(4, allAcc, output.accuracy); 
    allFinalWeights{i}=output.finalWeights
end


% optAcc=mean(squeeze(nanmean(allOptAcc(:,1,:,:), 3)), 2);

zeroNC_acc=mean(squeeze(nanmean(allAcc(:,1,:,:), 3)), 2);
someNC_acc=nanmean(squeeze(allAcc(:,2,:,:)), 3);





hold on
a=plot(zeroNC_acc, 'r', 'lineWidth', 2)
for j = 1:size(someNC_acc, 2)
    frac=j./11;
    z(j)=plot(someNC_acc(:,j), '-', 'color', [.9, .9, .9].*frac, 'lineWidth', 2) 

end




ff=legend([a, z(1), z(end)], 'No NC', '0.2 NC -- fixed unit signal & variance',  '0.2 NC -- fixed SNR')
set(ff, 'box', 'off', 'location', 'south')
ylabel('Accuracy')
xlabel('Trials')
set(gca, 'box', 'off')
saveas(gcf, 'suppFig_assumptions.eps', 'epsc2')
close all


save diffAssumptionsWorkspace_3-28-21.mat













