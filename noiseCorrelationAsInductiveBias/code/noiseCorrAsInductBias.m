% NOISE CORRELATIONS AS INDUCTIVE BIASES for learning.

% Code for figures 3&6 of:
% Noise Correlations for Faster and More Robust Learning
% Matthew R. Nassar, Daniel Scott, Apoorva Bhandari
% Journal of Neuroscience 4 August 2021, 41 (31) 6740-6752; DOI: 10.1523/JNEUROSCI.3045-20.2021 

% Set location of folder on local machine in "baseDir" on line 30
% Then code should run without additional changes. 
%% RUN THIS SECTION FOR EITHER SIMULATION?

clear

targetSEM       =  1;    % standard error on mean population representation
nPools          =  2;
nNeuronsPerPool = 100;
targFR          =  1;
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
        baseDir='~/Downloads';
    case 2
        baseDir='~/Dropbox (Brown)' % Put your dropbox path here
end
addpath(genpath(fullfile(baseDir, 'noiseCorrelationAsInductiveBias')));
cd(fullfile(baseDir, 'noiseCorrelationAsInductiveBias'));
%addpath(genpath('.'))


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
     

numReps               = 100;
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
    %output=runCorrLearningSim_noNorm(input);
    output=runCorrLearningSim(input); % MRN switched to non-normalized version
    
    allOptAcc=cat(3, allOptAcc, output.optAccuracy);
    allAcc=cat(3, allAcc, output.accuracy); 
    allFinalWeights{i}=output.finalWeights;

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



a=normpdf(-2:.5:2,0, 1 );
a=a./sum(a);

clear smMeanAcc
for i = 1:size(meanAcc, 2)

    
    normTerm=conv(ones(length(meanAcc(:,i)), 1), a', 'same');

    smMeanAcc(:,i)=conv(meanAcc(:,i), a', 'same')./normTerm;

end



nCorrs=length(input.all_fracCorrNoise);
hold on
clear testAcc testOptAcc testAccMat noiseCorrMat
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
make_perceptLearnFig




% CREATE Figure 2 based on perceptual learning simulations:








%% Network #2  -- three layers, hebbian learning reduces dimensionality and 
%                 produces noise correlations

% OK, now -- what if we create a slightly different network in which:

% 1) Layer one has uncorrelated noise as described above
% 2) Layer two receives inputs from layer one that are tuned with Hebbian learning
% 3) Layer two projects to outputs and weights between layers 2/3 are
%    learned using supervised or gated hebbian learning. 

% HEBBIAN LEARNING IS UNSTABLE -- need to make it more robust 

nTrials=200;
hebbLRs=[0, .0005];
rndDecWeightStd= .0000001;
rndInputWeightStd =.01;
noiseCorrWindSize=100;


preInitSelfWeights=false; % if false, will use random normal weights. 



% more trials here -- so recreate "stims"
allStims=randi(nPools, nTrials, 1); % pick one set of stims for all mods

% poolID

storeWeights=[1, 100];
noiseCorrsInPool=nan(nTrials, length(hebbLRs));
noiseCorrsOutPool=nan(nTrials, length(hebbLRs));

clear accuracy finalWeights2   firingRate firingRateLayer2 storeNoiseCorrMat finalWeights1
for k = 1:length(hebbLRs)
    
    
    hebbLR=hebbLRs(k)
    % Inputs will be truly uncorrelated:
    fracCorrNoise=0; %all_fracCorrNoise(k);
    
    % Choose an overall level of neural variance based on equation for fixed
    % lambda and tot variance:
    frVar = totVar ./ (nNeuronsPerPool + nNeuronsPerPool.*(nNeuronsPerPool-1).*fracCorrNoise);
    inPoolCov       = frVar.*fracCorrNoise; % some correlations within pools.
    outPoolCov      = 0;  % no correlations across pools.
    
    
    % create a matrix stipulating which neurons are in same pool:
    samePool=repmat(poolID, 1, length(poolID))==repmat(poolID, 1, length(poolID))';
    sameNeuron=logical(eye(length(poolID)));
    
    
    % create a covariance matrix across neural population:
    covMat=nan(size(sameNeuron));
    covMat(sameNeuron)=frVar;
    covMat(samePool&~sameNeuron)=inPoolCov;
    covMat(~samePool&~sameNeuron)=outPoolCov;
    
    %% GET READY:
    
    % preallocate space for firing rates:
    firingRate=nan(nNeuronsPerPool.*nPools, nTrials);
    firingRateLayer2=nan(nNeuronsPerPool.*nPools, nTrials);

    % create a random weight matrix that projects onto a "decision" neuron:
    wtMatrix2=normrnd(0, rndDecWeightStd, length(poolID),numOutNeurons);
   
    % And an N X N matrix connecting input to output
    % initialize to identity matrix plus some normal noise:
    if preInitSelfWeights
        wtMatrix1=normrnd(0, rndInputWeightStd, length(poolID),  length(poolID))+eye(length(poolID));
    else
        % random normal initialization:
         wtMatrix1=normrnd(0, rndInputWeightStd, length(poolID));        
    end
    
    
    
    stored_FR1=nan(nTrials, length(poolID));
    stored_FR2=nan(nTrials, length(poolID));
    
    
    % Loop through trials and simulate firing rates:
    Mu=nan(size(poolID));
    for i = 1:nTrials
        Mu(poolID==allStims(i))=targFR;
        Mu(poolID~=allStims(i))=nonTargFR;
        % compute and store firing rate for input layer neurons:
        FR1 = mvnrnd(Mu,covMat);
        firingRate(:,i)=FR1;
        
        % compute and store firing rate in layer 2:
        FR2=FR1*wtMatrix1; % each column reflects the weights IN to a layer 2 neuron
        firingRateLayer2(:,i)=FR2; % store firing rates... 

%       compute activation of downstream neuron:
        decNeuronOutput=FR2*wtMatrix2;
        
        if ~isfinite(decNeuronOutput)
            disp('houston, we have a problem!')
            keyboard
        end
        
        % make a choice:
        if ~doRL
            [pChoice choice(i)]=softMax_epsilon([decNeuronOutput, 0], invT, 0);
            expRew=pChoice(choice(i));
        else
            [pChoice choice(i)]=softMax_epsilon([decNeuronOutput], invT, 0);
            if useRPE
            expRew=pChoice(choice(i));
            else
            expRew=.5;
            end
        end

   
        % store running average of in pool and out of pool noise
        % correlations:
        if i>noiseCorrWindSize
            sel=(i-noiseCorrWindSize+1):i;
            winFR=firingRateLayer2(:,sel)';
            xMat=[ones(length(sel), 1), allStims(sel)];
            
            % do regression to get betas:
            betas=(xMat'* xMat)\(xMat'*winFR); 
            residual=   winFR-(xMat*betas);
            
            noiseCorrMat=corr(residual);
            noiseCorrsInPool(i, k)=nanmean(noiseCorrMat(samePool&~sameNeuron));
            noiseCorrsOutPool(i, k)=nanmean(noiseCorrMat(~samePool&~sameNeuron));
        else
            clear noiseCorrMat
        end

     % If you already knew the best decision rule:
        optOutput=FR1*optWeights;

        meanFR=(targFR+nonTargFR)./2;
        % UPDATE layer 1 to layer 2 connections with hebbian learning:
        deltaWt=hebbLR.*(FR1-meanFR)'*(FR2-meanFR); % get 

        if isfinite(wtMatrix1+deltaWt)
            wtMatrix1=wtMatrix1+deltaWt;   % Do it. 
        else 
           disp('here comes trouble' )
           keyboard
        end

        % Normalize weight matrix (eg. homeostatic plasticity)                
        wtMatrix1=wtMatrix1 ./repmat(sqrt(sum(wtMatrix1.^2)),  length(poolID),1);
         

        if ~doRL
            % compute error:
            PE=targetOutputs(allStims(i))-decNeuronOutput;
            % update weights:
            wtMatrix2=wtMatrix2+ LR.*FR2'.*PE;
            storePE(i)=PE;
        else
            % deltaW(C) = alpha * RPE * [X-E(X)]
            RPE=(choice(i)==allStims(i)) - expRew; % compute RPE
            deltaWt=LR.*RPE.*(FR2-meanFR); % compute weight update
            wtMatrix2(:,choice(i))= wtMatrix2(:,choice(i))+deltaWt'; %implement weight update
            storePE(i)=RPE;  % store RPE
          %  wtMatrix2(:,choice(i))=wtMatrix2(:,choice(i))./(sum(wtMatrix2(:,choice(i)).^2)./weightScale); % normalize weights
        end
        
        
         
        
        
        %
        %         % compute error:
        %         PE=targetOutputs(allStims(i))-decNeuronOutput;
        %
        %         % update weights:
        %         wtMatrix2=wtMatrix2+ LR.*FR2'.*PE;
        %
        %
        %
        accuracy(i,k)=choice(i)==allStims(i);
        %         optAccuracy(i,k)=sign(optOutput)==sign(targetOutputs(allStims(i)));
        %         storePE(i)=PE;
        %
        
        if any(i==storeWeights)
            ind=find(i==storeWeights);
            
            finalWeights1(:,:,ind,k)=wtMatrix1;
            
            if exist('noiseCorrMat')
                storeNoiseCorrMat(:,:,ind,k)=noiseCorrMat;
            else
                storeNoiseCorrMat(:,:,ind,k)=nan(length(poolID));
            end
            

        end
    end
    if ~doRL
        finalWeights2(:, k)=wtMatrix2;
    else
        finalWeights2(:, :, k)=wtMatrix2;
    end
end


selStart=[100];
selEnd=[200];



allWeights=wtMatrix2(:,1)-wtMatrix2(:,2);
% ok, what if we resort the neurons according to the connections to output
% neurons:
 [B,I] = sort(allWeights) ;



for i = 1:length(selEnd)
    selDat=selStart(i):selEnd(i);
    subplot(length(selStart), 2, i)
    clear res_FR2 yHat2
    for k = 1:size(firingRateLayer2, 1)
        xes=[ones(length(selDat), 1), allStims(selDat)];
        [B,BINT,res_FR2(:,k)] = regress(firingRateLayer2(k,selDat)',xes);
        yHat2(:,k)=xes*B;
    end
    noiseCorrMat=corr((res_FR2));
    
    noiseCorrMatSort=corr((res_FR2(:,I)));
    imagesc(noiseCorrMat, [-1, 1])
    colorbar

    selDat=selStart(i):selEnd(i);
    subplot(length(selStart), 2, i+1)
    clear res_FR1 yHat1
    for k = 1:size(firingRate, 1)
        xes=[ones(length(selDat), 1), allStims(selDat)];
        [B,BINT,res_FR1(:,k)] = regress(firingRate(k,selDat)',xes);
        yHat1(:,k)=xes*B;
    end
    noiseCorrMat=corr((res_FR1));
    imagesc(noiseCorrMat, [-1, 1])
    colorbar

end


inPool=[true(nNeuronsPerPool,nNeuronsPerPool,1), false(nNeuronsPerPool,nNeuronsPerPool,1); ...
    false(nNeuronsPerPool,nNeuronsPerPool,1), true(nNeuronsPerPool,nNeuronsPerPool,1)];

sameNeuron=eye(nNeuronsPerPool.*2);

inPoolCorrs=nanmean(noiseCorrMat(inPool&~sameNeuron))
outPoolCorrs=nanmean(noiseCorrMat(~inPool&~sameNeuron))



% what does read out of layer 2 look like:

subplot(2, 2, 3)
title('Noise correlation')
imagesc((corr((res_FR1))), [-1, 1])
colorbar
subplot(2, 2, 1)
title('Signal correlation')
imagesc((corr((yHat1))), [-1, 1])
colorbar    


subplot(2, 2, 4)
title('Noise correlation')
imagesc((corr((res_FR2(:,I)))), [-1, 1])
colorbar
subplot(2, 2, 2)
title('Signal correlation')
imagesc((corr((yHat2(:,I)))), [-1, 1])
colorbar 

close all


% Make figure

if preInitSelfWeights   
make_hebbNoiseCorrFig
else
make_hebbNoiseCorrFig_randomWeights
end


