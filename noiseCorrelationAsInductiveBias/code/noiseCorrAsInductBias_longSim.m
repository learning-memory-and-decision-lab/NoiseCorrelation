% NOISE CORRELATIONS AS INDUCTIVE BIASES for learning?

% Run one long simulation to show that low noise correlation conditions do
% eventually learn near-optimal readout. 



%% RUN THIS SECTION FOR EITHER SIMULATION?

whichComp       =  1;
targetSEM       =  1;    % standard error on mean population representation
nPools          =  2;
nNeuronsPerPool = 100;
targFR          =  1;
nonTargFR       = -1; % 
LR              =.0001;
doRL            =true;
weightScale     =1;
useRPE          =false; % online predictions
doOpt           =true;

% Sort out path:
switch whichComp
    case 1
        baseDir='~/Dropbox';
    case 2
        baseDir='~/Dropbox (Brown)' % Put your dropbox path here
end
addpath(genpath(fullfile(baseDir, 'sharedMatlabUtilities')));
addpath(genpath(fullfile(baseDir, 'noiseCorrelationAsInductiveBias')))
cd(fullfile(baseDir, 'noiseCorrelationAsInductiveBias'));


totVar = nNeuronsPerPool.*targetSEM.^2; % FIXED to a target level!
all_fracCorrNoise= 0:.1:1;

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
     

numReps               = 1000;
testTrials            = 4981:5000;
%input.targetSEM       = 10;    % standard error on mean population representation
input.sigToNoise      = 2;
input.nPools          = 2;
input.nNeuronsPerPool = 100;
input.targFR          = targFR;
input.nonTargFR       = nonTargFR;
input.nTrials         = 5000;
input.LR              = LR;
input.doRL            = doRL;
input.invT            = 10000;
input.weightScale     = weightScale;
input.useRPE          = useRPE; % this produces weird results -- leave off for now.
input.doOpt           = doOpt
input.all_fracCorrNoise=0:.01:.2;


nCorrs=length(input.all_fracCorrNoise);

allOptAcc=[]; allAcc=[]; allFinalWeights={};
for i = 1:numReps
    disp(sprintf('running rep %g', i));
    output=runCorrLearningSim(input);
    allOptAcc=cat(3, allOptAcc, output.optAccuracy);
    allAcc=cat(3, allAcc, output.accuracy); 
    allFinalWeights{i}=output.finalWeights

end

if doRL
weightDiff=[squeeze(output.finalWeights(:,1,:))-squeeze(output.finalWeights(:,2,:))] ;
normWeightDiff=weightDiff./repmat(std(weightDiff), size(weightDiff, 1), 1) ;
optWeightDiff=[(optWeights(:,1,:))-(optWeights(:,2,:))] ;
normWeightDiff=optWeightDiff./repmat(std(optWeightDiff), size(optWeightDiff, 1), 1) ;
else
weightDiff=finalWeights;   
end


% Make plots of the "final" weights from 3 different noise correlations:

toShow=[1, 11, 21]
input.all_fracCorrNoise(toShow)
% 
% figure
% hold on
% subplot(3, 1, 1)
% hist([normWeightDiff(1:100,toShow(1)), normWeightDiff(101:200,toShow(1))], -3:.3:3)
% xlim([-4, 4])
% set(gca, 'box', 'off')
% 
% ylabel(num2str(input.all_fracCorrNoise(toShow(1))))
% 
% subplot(3, 1, 2)
% hist([normWeightDiff(1:100,toShow(2)), normWeightDiff(101:200,toShow(2))], -3:.3:3)
% xlim([-4, 4])
% ylabel(num2str(input.all_fracCorrNoise(toShow(2))))
% set(gca, 'box', 'off')
% 
% subplot(3, 1, 3)
% hist([normWeightDiff(1:100,toShow(3)), normWeightDiff(101:200,toShow(3))], -3:.3:3)
% xlim([-4, 4])
% set(gca, 'box', 'off')
% ylabel(num2str(input.all_fracCorrNoise(toShow(3))))
% 

% MRN STOPPED HERE!!!!
% Now go through and grab final weights for ALL simulations:
eucDistToBound=[]; idealIsAccurate=[]; optDistToBound=[];
for rep=1:1:numReps
    fw=allFinalWeights{rep};
    
    if doRL
        weightDiff=[squeeze(fw(:,1,:))-squeeze(fw(:,2,:))] ;
        normWeightDiff=weightDiff./repmat(std(weightDiff), size(weightDiff,1), 1) ;
        optWeightDiff=[(optWeights(:,1,:))-(optWeights(:,2,:))] ;
        normOptDiff=optWeightDiff./repmat(std(optWeightDiff), size(optWeightDiff, 1), 1) ;  
    else
        weightDiff=finalWeights;
    end


  % 
    C = 1;
    
    Mu(poolID==C)=targFR;
    Mu(poolID~=C)=nonTargFR;
    FR = Mu; % Start at a "perfect representation of the stimulus)
    
    
    for k = 1:nCorrs
        
        kWeights=weightDiff(:,k);
        normPred=FR*kWeights;

        idealIsAccurate(k,rep)=sign(normPred)==sign(targetOutputs(C));
        
        % how much would you need to change each neurons output to get to the
        % boundary?
        
        distToBoundary_perDim=abs(normPred./kWeights);
        % in which direction:
        dirToBoundary=sign(normPred./kWeights); % only important if we want to create adversarial examples.
        
        minDistToBound=distToBoundary_perDim(1);
        for dim=2:length(kWeights);
            % Consider a new dimension:
            newDimDistToBound=distToBoundary_perDim(dim);
            distRatio=newDimDistToBound./minDistToBound;
            minDistToBound=minDistToBound.*distRatio./sqrt(distRatio^2+1);
        end
        
        % If we're already past boundary -- make distance negative.
        if idealIsAccurate(k,rep)
            eucDistToBound(k,rep)= minDistToBound;
        else
            eucDistToBound(k,rep)= -minDistToBound;
        end
                
    end
    
    % repeat same procedure, but for optimal weights:
    
    
    
    kWeights=optWeightDiff;
    normPred=FR*kWeights;

    % how much would you need to change each neurons output to get to the
    % boundary?
    
    distToBoundary_perDim=abs(normPred./kWeights);
    % in which direction:
    dirToBoundary=sign(normPred./kWeights); % only important if we want to create adversarial examples.
    
    minDistToBound=distToBoundary_perDim(1);
    for dim=2:length(kWeights);
        % Consider a new dimension:
        newDimDistToBound=distToBoundary_perDim(dim);
        distRatio=newDimDistToBound./minDistToBound;
        minDistToBound=minDistToBound.*distRatio./sqrt(distRatio^2+1);
    end
    
    optEucDistToBound = minDistToBound;

end


hold on
STD=nanstd(eucDistToBound')
H=shadedErrorBar(input.all_fracCorrNoise, nanmean(eucDistToBound'), STD, {'color', [.8, 0, .3]});
plot([input.all_fracCorrNoise(1), input.all_fracCorrNoise(end)], [optEucDistToBound, optEucDistToBound], 'g')


%plot(input.all_fracCorrNoise, eucDistToBound(:,1))
ylabel('Distance to category boundary')
xlabel('Noise correlations')
set(gca, 'box', 'off')

% 
% 
% close all
% subplot(3, 1, 1)
% hold on
% %a=plot(optWeights, '--k');
% for i = 1:(nCorrs)
%     plot(weightDiff(:,i), 'color', ones(3,1).*i./(nCorrs+1));
% end
% ylabel('weight')
% xlabel('neuron')
% 
% %ff=legend([a, b], 'optimal', 'learned');
% 
% subplot(3, 1, 2)
% 



% 
% hold on
% c=plot(all_fracCorrNoise, nanmean(allOptAcc), 'b')
% d=plot(all_fracCorrNoise, nanmean(accuracy), 'r')
% aa=legend([c, d], 'Optimal readout', 'Learned readout')
% set(aa, 'location' ,'southeast', 'box', 'off')
% ylabel('Accuracy')
% 
% 
% subplot(3, 1, 3)
% 
% hold on
% e=plot(all_fracCorrNoise, nanmean(eucDistToBound,2), 'r')
% ylabel('Distance to bound')
% xlabel('In pool correlation')
% % aa=legend([c, d], 'Optimal readout', 'Learned readout')
% % set(aa, 'location' ,'southeast', 'box', 'off')
% % 
% saveas(gcf, 'noiseCorrsAsInductBias.eps', 'epsc2')

% optimal accuracy:
meanOptAcc=nanmean(allOptAcc, 3);
semOptAcc=nanstd(allOptAcc, [], 3)./sqrt(size(allOptAcc, 3));

% actual accuracy:
meanAcc=nanmean(allAcc, 3);
semAcc=nanstd(allAcc, [], 3)./sqrt(size(allAcc, 3));


% Do they look same?
firstHalfMean=nanmean(allOptAcc(:,:,1:500), 3);
secondHalfMean=nanmean(allOptAcc(:,:,501:1000), 3);

hold on
plot(firstHalfMean)
plot(secondHalfMean,'r')
corr(firstHalfMean(900:1000,1), secondHalfMean(900:1000,1))




a=normpdf(-2:.5:2,0, 1 );
a=a./sum(a);

clear smMeanAcc
for i = 1:size(meanAcc, 2)

    
    normTerm=conv(ones(length(meanAcc(:,i)), 1), a', 'same');

    smMeanAcc(:,i)=conv(meanAcc(:,i), a', 'same')./normTerm;

end





nCorrs=length(input.all_fracCorrNoise);
hold on
clear testAcc testOptAcc testAccMat
%a=plot(optWeights, '--k');
for i = 1:(nCorrs)
    plot(meanAcc(:,i), 'color', ones(3,1).*i./(nCorrs+1));
    testAcc(i)=nanmean(meanAcc(testTrials,i));
    testOptAcc(i)=nanmean(meanOptAcc(testTrials,i));
    testAccMat(:,i)=squeeze(mean(allAcc(testTrials,i,:)));
    noiseCorrMat(:,i)=ones(size(squeeze(mean(allAcc(testTrials,i,:))))).*input.all_fracCorrNoise(i);
end

% MRN added 7/8 to get stats. 
[RHO,PVAL] =corr(noiseCorrMat(:), testAccMat(:));
flNoiseCorrMat=noiseCorrMat';
[RHO2,PVAL2] =corr(flNoiseCorrMat(:), eucDistToBound(:));







hold on
c=plot(input.all_fracCorrNoise, testOptAcc, 'b')
d=plot(input.all_fracCorrNoise, testAcc, 'r')


aa=legend([c, d], 'Optimal readout', 'Learned readout')
set(aa, 'location' ,'southeast', 'box', 'off')
ylabel('Accuracy')
xlabel('Fraction correlated noise')

% 
% subplot(3, 1, 1)
% hold on
% %a=plot(optWeights, '--k');
% for i = 1:(nCorrs)
%     plot(normWeights(:,i), 'color', ones(3,1).*i./(nCorrs+1));
% end
% ylabel('normalized weight')
% xlabel('neuron')
% 


%make_CCN_fig
%make_perceptLearnFig


close all
nCorrs=length(input.all_fracCorrNoise);

hold on
plot([0, 8.5], [nanmean(testOptAcc(:)), nanmean(testOptAcc(:))], '-', 'color', cbColors(7,:))

clear testAcc testOptAcc
%a=plot(optWeights, '--k');
for i = (nCorrs):-1:1
    plot(log(1:length(meanAcc(:,i))), meanAcc(:,i), 'color', ones(3,1).*i./(nCorrs+3));
    testAcc(i)=nanmean(meanAcc(testTrials,i));
    testOptAcc(i)=nanmean(meanOptAcc(testTrials,i));
end

xlim([0, 8.5])

ylabel('Accuracy')
xlabel('Time (log scale)')
set(gca, 'box', 'off')
saveas(gcf, 'longSimSuppFig.eps', 'epsc2')
close all

        



make_perceptLearnFig_longSim



% CREATE Figure 2 based on perceptual learning simulations:


save longSimResults_3-27-21.mat



