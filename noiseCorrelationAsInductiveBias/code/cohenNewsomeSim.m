% Code for figure 7 of:
% Noise Correlations for Faster and More Robust Learning
% Matthew R. Nassar, Daniel Scott, Apoorva Bhandari
% Journal of Neuroscience 4 August 2021, 41 (31) 6740-6752; DOI: 10.1523/JNEUROSCI.3045-20.2021 

% Model more complex task (abstract version of Cohen/newsome task):
% step 1 = do basic x,y,+,- task encoding with 2 layers... show that
% correlations ala cohen and newsome produce faster learning than fixed
% correlations across trials

% Instructions for reproducing results:
% 1) navigate to noiseCorrelationAsInductiveBias folder on local machine. 
% 2) run script. Note that for paper, nReps (line 49) was set to 10,000,
% but using this value requires substantial time for simulation. 

clear

addpath(genpath('.'))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              Set parameterization for model:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup to run both relevant and irrelevant correlations 
% -- use corrType =1 for relevant dimension correlations and 2 for
% irrelevant

corrTypes=[1,2]

for corrType=1:2




targetSEM       = 20;    % standard error on mean population representation
nPools          = 4;    % 4 =  x+, x-, y+, y-
nNeuronsPerPool = 100;  % plenty neurons per pool
targFR          = 1;    % now targFR controls the scale of tuning for xy and +-
nTrials         = 100;
LR              =.005;
doRL            =true;  % RIGHT NOW THIS CODE ONLY USES RL type LR... learn with RL rather than with supervised signal
invT            =10000; % pretty deterministic...
weightScale     =1;     % normalize weights to some level...
useRPE          =false; % this produces weird results -- leave off for now.
doOpt           =true;
initWtVar       =0.00001;
whichComp       =1;     % who is running simulation (matt =1 , apoorva = 2)
meanFR          =0;     % currently not doing much...
taskNeuronWt    =1000;    % strength of hard coded inputs from "task neurons" 
nReps           =100;  % 10000   % ran with 1000, should ramp up to 10,000 or so overnight for paper
testTrials      =80:100;   
normType        =3;     % 1 = divide by squared weight, 2 = divide by standard deviation of weights (eg. sqrt(sum(wt^2)./n), 3 = no normalization. 
ms              =8;
%corrType        =1;     % 1 = relevant pool, 2 = irrelevant pool

chooseMax       =true;  % Use greedy action selection, with perfect task bias. 
makeSchem       =false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Sort out path:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch whichComp
    case 1
        baseDir='~/Downloads';
    case 2
        baseDir='~/Dropbox (Brown)' % Put your dropbox path here
end
addpath(genpath(fullfile(baseDir, 'sharedMatlabUtilities')));
cd(fullfile(baseDir, 'noiseCorrelationAsInductiveBias'));

if ~doRL
    numOutNeurons = 1;
else
    numOutNeurons = nPools;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                         Setup the correlation structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a population of tuned neurons:
poolLabels = {'x+', 'x-', 'y+', 'y-'};
poolXTune     = [1, 1, -1, -1].*targFR;
poolPlusTune  = [1, -1, 1, -1].*targFR;

nTrialTypes = 2; % 1 = xy, 2 = +-
% output neurons report x/y for task 1, or +/- for task 2:
outputNeuronLabels={'1x', '1y', '2+', '2-'};
outputTask1 = [1,1, -1, -1];


% what do we want correlations to look like?
% We'll control 3 types of correlations:

corrStruct.inPool  = 0:.04:.2;  % level of correlations between neurons in same pool
if corrType ==1;     % 1 = relevant pool, 2 = irrelevant pool
    corrStruct.inRel   =0:.04:.2;  % level of correlations between neurons in diff pools with similar encoding of relevant stimulus dimension
    corrStruct.inIrrel = 0;  % level of correlations between neurons in diff pools with similar encoding of relevant stimulus dimension
elseif corrType == 2;
    corrStruct.inRel   = 0;  % level of correlations between neurons in diff pools with similar encoding of relevant stimulus dimension
    corrStruct.inIrrel = 0:.04:.2;  % level of correlations between neurons in diff pools with similar encoding of relevant stimulus dimension
end


% Make all possible combinations of correlation structure:
allCorrCombos=makeAllCombos(corrStruct);

% get rid of correlation structures that aren't possible:
allCorrCombos=selBehav(allCorrCombos, allCorrCombos.inPool>=allCorrCombos.inIrrel & allCorrCombos.inPool>=allCorrCombos.inRel)


% How much variance should be in summed output?
totVar = nNeuronsPerPool.*2.*targetSEM.^2; % There are always two pools contributing to each trial



% create a list of neuron pool assignments:
poolID=[];
for i = 1:nPools
    poolID(end+1:end+nNeuronsPerPool,1)=ones(nNeuronsPerPool,1).*i;
end

% Make maps of each type of covariance:
% poolLabels = {'x+', 'x-', 'y+', 'y-'}; % JUST FOR REFERENCE:
isX=false(size(poolID)); isPlus=false(size(poolID));
isX(poolID==1|poolID==2)=true;
isPlus(poolID==1|poolID==3)=true;

% OK -- what should weights look like?
optWeights=[]

%poolLabels = {'x+', 'x-', 'y+', 'y-'};

% trial types:
% 1x, 1y, 2+, 2-   % Outputs 1-4

% Hand code for now:
% 1x
optWeights(poolID==1|poolID==2,1)=1;
optWeights(poolID==3|poolID==4,1)=-1;

% 1y
optWeights(poolID==1|poolID==2,2)=-1;
optWeights(poolID==3|poolID==4,2)=1;
optWeights(:,2)=optWeights(:,2)./(sum(optWeights(:,2).^2)./weightScale);

% 2+
optWeights(poolID==1|poolID==3,3)=1;
optWeights(poolID==2|poolID==4,3)=-1;
optWeights(:,3)=optWeights(:,3)./(sum(optWeights(:,3).^2)./weightScale);


% 2-
optWeights(poolID==1|poolID==3,4)=-1;
optWeights(poolID==2|poolID==4,4)=1;
optWeights(:,4)=optWeights(:,4)./(sum(optWeights(:,4).^2)./weightScale);


for i = 1:4
    if normType==1
    optWeights(:,i)=optWeights(:,i)./(sum(optWeights(:,i).^2)./weightScale);
    elseif normType==2
    optWeights(:,i)=optWeights(:,i)./  sqrt(nanmean(optWeights(:,i).^2));
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               LOOP through trials, do task, and learn!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear finalWeights


for rep = 1:nReps
    disp(sprintf('Repeat number %g', rep));
    
    % initialize task variables:
    % MRN moved this into loop 4-29-20
    
    % Choose X/Y for each trial: 1 = x, -1 = y
    allStimsX=(randi(nPools./2, nTrials, 1)-1.5).*2; % pick one set of stims for all mods
    % Choose +/- for each trial: 1 = +, -1 = -
    allStimsPlus=(randi(nPools./2, nTrials, 1)-1.5).*2; % pick one set of stims for all mods
    % Choose trial type (XY versus +-) for each trial:
    allTrialTypes=randi(nTrialTypes, nTrials, 1);
    
    % Figure out correct output for each trial:
    corrOutput=zeros(size(allTrialTypes));
    corrOutput(allTrialTypes==1&allStimsX==1)=1;
    corrOutput(allTrialTypes==1&allStimsX==-1)=2;
    corrOutput(allTrialTypes==2&allStimsPlus==1)=3;
    corrOutput(allTrialTypes==2&allStimsPlus==-1)=4;
    
    for k = 1:length(allCorrCombos.inPool)
        
        
        
        % Pick noise correlation structure:
        simCorrs.inPool=allCorrCombos.inPool(k);
        simCorrs.inRel=allCorrCombos.inRel(k);
        simCorrs.inIrrel=allCorrCombos.inIrrel(k);
        
        
        % Choose an overall level of neural variance based on equation for fixed
        % lambda and tot variance:
        
        
        
        %Mean variance = tot variance / [number neurons + num covariance elements * correlations between those elements]
        %  frVar1 = totVar ./ (nNeuronsPerPool + nNeuronsPerPool.*(nNeuronsPerPool-1).* simCorrs.inPool);
        
        %  frVar2 = totVar ./ (nNeuronsPerPool + nNeuronsPerPool.*(nNeuronsPerPool-1).* simCorrs.inPool);
        
        % total variance for a pool / [num neurons + numWithinPoolPairs*inPoolCorrs
        % number of crossPoolPairs *
        frVar_2Pools=totVar ./ (nNeuronsPerPool + nNeuronsPerPool.*(nNeuronsPerPool-1).*simCorrs.inPool + nNeuronsPerPool.*(nNeuronsPerPool).*simCorrs.inRel - nNeuronsPerPool.*(nNeuronsPerPool).* simCorrs.inIrrel);
        %keyboard
                                                                                                      
        

        % Now this is the sum of variance across the 2 pools encoding the
        % stim...
        
        % OK -- if we wanted to do this but include all covariance elements --
        % what would it look like?
        
       
        inPoolCov        = frVar_2Pools.*simCorrs.inPool;  % some correlations within pools.
        relPoolCov       = frVar_2Pools.*simCorrs.inRel;   % some correlations in pools that share tuning for relevant dimension
        irrelPoolCov     = frVar_2Pools.*simCorrs.inIrrel; % for now lets just make this zero.
        outPoolCov       = 0;                              % everything thats not being set...

        % create a matrix stipulating which neurons are in same pool:
        samePool=repmat(poolID, 1, length(poolID))==repmat(poolID, 1, length(poolID))';
        sameNeuron=logical(eye(length(poolID)));
        sameRel1=repmat(isX, 1, length(poolID))==repmat(isX, 1, length(poolID))';
        sameRel2=repmat(isPlus, 1, length(poolID))==repmat(isPlus, 1, length(poolID))';
        
        
        if corrType ==1;     % 1 = relevant pool, 2 = irrelevant pool

            % create a covariance matrix for task 1 (x/y):
            covMat1=nan(size(sameNeuron));
            covMat1(sameNeuron)=frVar_2Pools;
            covMat1(samePool&~sameNeuron)=inPoolCov;
            covMat1(sameRel1&~samePool)=relPoolCov;
            covMat1(~samePool&~sameNeuron&~sameRel1)=outPoolCov;
            
            covMat2=nan(size(sameNeuron));
            covMat2(sameNeuron)=frVar_2Pools;
            covMat2(samePool&~sameNeuron)=inPoolCov;
            covMat2(sameRel2&~samePool)=relPoolCov;
            covMat2(~samePool&~sameNeuron&~sameRel2)=outPoolCov;
            
            
            
            
            if makeSchem 
                keyboard
                figure
                imagesc(covMat1)
                saveas(gcf, 'covMat1.eps', 'epsc2')
                close all
                
            end
            
            
            
            
        elseif corrType==2 % If we're making correlations to irrelevant pool. 
            % create a covariance matrix for task 1 (x/y):
            covMat1=nan(size(sameNeuron));
            covMat1(sameNeuron)=frVar_2Pools;
            covMat1(samePool&~sameNeuron)=inPoolCov;
            covMat1(sameRel2&~samePool)=irrelPoolCov;
            covMat1(~samePool&~sameNeuron&~sameRel2)=outPoolCov;
            
            covMat2=nan(size(sameNeuron));
            covMat2(sameNeuron)=frVar_2Pools;
            covMat2(samePool&~sameNeuron)=inPoolCov;
            covMat2(sameRel1&~samePool)=irrelPoolCov;
            covMat2(~samePool&~sameNeuron&~sameRel1)=outPoolCov;
        end
        
            
        %% GET READY:
        
        % preallocate space for firing rates:
        firingRate=nan(nNeuronsPerPool.*nPools, nTrials);
        % create a random weight matrix that projects onto a single "decision" neuron:
        wtMatrix=normrnd(0, initWtVar, length(poolID),numOutNeurons);
        
        
        % Loop through trials and simulate firing rates:
        Mu=nan(size(poolID));
        for i = 1:nTrials
             
            Mu=nan(size(poolID));
            poolMu=allStimsX(i)*poolXTune + allStimsPlus(i)*poolPlusTune;
            for pp = 1:max(poolID)
                Mu(poolID==pp)=poolMu(pp);
            end

            
            % Select covariance matrix appropriate to trial type:
            if allTrialTypes(i)==1
                FR = mvnrnd(Mu,covMat1);
            else
                FR = mvnrnd(Mu,covMat2);
            end
            firingRate(:,i)=FR;

            
            % compute activation of downstream neuron(s):
            decNeuronOutput=FR*wtMatrix;
            
            
            if ~chooseMax
            % Bias decision neurons toward relevant task:
            decNeuronOutput=decNeuronOutput+ outputTask1*2.*((allTrialTypes(i)==1)-.5).*taskNeuronWt;
            [pChoice choice(i)]=softMax_epsilon([decNeuronOutput], invT, 0);
            else
                
            decNeuronOutput(outputTask1~=2.*((allTrialTypes(i)==1)-.5))=-inf;
            [Y, choice(i)] = max(decNeuronOutput);
            end
           
            if useRPE
                expRew=pChoice(choice(i));
            else
                expRew=.5;
            end

            if ~chooseMax&any(~isfinite(decNeuronOutput))
                disp('houston, we have a problem!')
                keyboard
            end
            
            % If you already knew the best decision rule:
            if doOpt
                if ~chooseMax
                    optOutput=FR*optWeights + outputTask1*2.*((allTrialTypes(i)==1)-.5).*taskNeuronWt;
                else
                    optOutput=FR*optWeights;
                    optOutput(outputTask1~=2.*((allTrialTypes(i)==1)-.5))=-inf;
                end
            end
            
            
            Rew=choice(i)==corrOutput(i);
            
            % deltaW(C) = alpha * RPE * [X-E(X)]
            RPE=Rew - expRew; % compute RPE
            deltaWt=LR.*RPE.*(FR'-meanFR); % compute weight update
            wtMatrix(:,choice(i))= wtMatrix(:,choice(i))+deltaWt; %implement weight update
            storePE(i)=RPE;  % store RPE
            
            if normType ==1
            wtMatrix(:,choice(i))=wtMatrix(:,choice(i))./(sum(wtMatrix(:,choice(i)).^2)./weightScale); % normalize weights
            elseif normType ==2
            wtMatrix(:,choice(i))=wtMatrix(:,choice(i))./ sqrt(nanmean(wtMatrix(:,choice(i)).^2)); % normalize weights
            end    
            
            accuracy(i,k)=Rew;
            
            if doOpt
                [val,ind]=max(optOutput);
                optAccuracy(i,k)=corrOutput(i)==ind;
            end

            storedWeightUpdates(:,i)=deltaWt;

        end
        
         allStoredWeightUpdates(:,:,k)=storedWeightUpdates;
         
         if normType<=2
         finalWeights(:,:,k)=wtMatrix;
         finalWeightDist(k)=mean(mean((finalWeights(:,:,k)-optWeights).^2));
         else
             
             % Normalize both optimal and learned weights such that the
             % expected euclidean norm is equal to one:
             normWeights=wtMatrix./repmat(sqrt(sum(wtMatrix.^2)), length(wtMatrix),1);
             finalWeights(:,:,k)=normWeights;
             normOptWeights=   optWeights./repmat(sqrt(sum(optWeights.^2)), length(optWeights),1);
             finalWeightDist(k)=mean(sqrt(sum((finalWeights(:,:,k)-normOptWeights).^2)));
         end
         
         
         
    end
    % Can we compute how: it would be good to estimate


    % Store mean accuracy results in a matrix:
    if corrType==1
        outCorrs=allCorrCombos.inRel;
    else
        outCorrs=allCorrCombos.inIrrel;
    end

    a=unique(allCorrCombos.inPool);
    
    for i = 1:length(a)
        for j = 1:length(a)
            sel= allCorrCombos.inPool==a(i)&outCorrs==a(j);
            if sum(sel)>0
                meanAccMat(i,j,rep)=nanmean(accuracy(testTrials,sel));
                optAccMat(i,j,rep)= nanmean(optAccuracy(testTrials,sel));
                meanDistMat(i,j,rep)=finalWeightDist(sel);
            end
        end
    end
    
    
    allAccuracy(rep,:,:)=accuracy;
    allOptAccuracy(rep,:,:)=optAccuracy;
end


defaultPlotParameters




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            MAKE plots to summarize simulated data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for k =1:21

updateInDims=[ones(1,200), ones(1,200).*-1; ones(1,100), ones(1,100).*-1, ones(1,100), ones(1,100).*-1]*allStoredWeightUpdates(:,:,k)


figure
hold on
Scale=max(abs(updateInDims(:)));
title(sprintf('in pool= %g, rel pool= %g, irrel pool = %g', allCorrCombos.inPool(k), allCorrCombos.inRel(k), allCorrCombos.inIrrel(k)))
plot([-Scale, Scale], [0, 0], '--k')
plot([0, 0], [-Scale, Scale],  '--k')
% a=plot(updateInDims(1,allTrialTypes==1), updateInDims(2,allTrialTypes==1), 'or', 'markerFaceColor','r', 'markerEdgeColor', 'k', 'lineWidth', 1, 'markerSize', 8)
% b=plot(updateInDims(1,allTrialTypes==2), updateInDims(2,allTrialTypes==2), 'ob', 'markerFaceColor','b', 'markerEdgeColor', 'k', 'lineWidth', 1, 'markerSize', 8)
a=quiver(zeros(sum(allTrialTypes==1),1), zeros(sum(allTrialTypes==1),1), updateInDims(1,allTrialTypes==1)', updateInDims(2,allTrialTypes==1)', 0,  'r', 'lineWidth', 2)
b=quiver(zeros(sum(allTrialTypes==2),1), zeros(sum(allTrialTypes==2),1), updateInDims(1,allTrialTypes==2)', updateInDims(2,allTrialTypes==2)', 0,  'b', 'lineWidth', 2)

ff=legend([a, b], 'X or Y', '+ or -')
set(ff, 'box', 'off')
xlabel('Update in XY')
ylabel('Update in +-')
xlim([-Scale, Scale])
ylim([-Scale, Scale])
set(gca, 'box', 'off')
fn=sprintf('WeightUpdateQuiver_%g_rel_%g_irrel_%g.eps', allCorrCombos.inPool(k),allCorrCombos.inRel(k), allCorrCombos.inIrrel(k));
saveas(gcf, fn, 'epsc2')
close all

end


% PLOT RAW WEIGHTS?
% k=21
% 
% hold on
% plot(optWeights(:,1)./std(optWeights(:,1)))
% plot(finalWeights(:,1,k)./std(finalWeights(:,1,k)), 'r')
% 



getCbColors

allCorrs=unique(allCorrCombos.inPool)
hold on
for i =1:length(allCorrs)
    sel=allCorrCombos.inPool==allCorrs(i)
    h(i)=plot( allCorrCombos.inRel(sel),   finalWeightDist(sel), '-o', 'color', cbColors(i+1,:), 'markerSize', ms, 'markerFaceColor', cbColors(i+1,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
    label{i}=num2str(allCorrs(i));
end
ylabel('Distance from optimal readout')
xlabel('Relevant feature correlations')
legend(h, label)
close all



figure
subplot(1, 2, 1)
imagesc(nanmean(meanAccMat, 3),[.5,1])
colorbar
title('Learned readout')
a=xlabel('Irrelevent pool correlations')
set(a, 'HorizontalAlignment', 'left')
ylabel('Same pool correlations')
set(gca, 'xtick', 2:2:6, 'xticklabel', label(2:2:6), 'ytick', 2:2:6, 'yticklabel', label(2:2:6)) 
subplot(1, 2, 2)
imagesc(nanmean(optAccMat, 3),[.5,1])
set(gca, 'xtick', 2:2:6, 'xticklabel', label(2:2:6), 'ytick', 2:2:6, 'yticklabel', []) 
title('Optimal readout')
colorbar
fn=sprintf('accuracyHeatmap_corrType_%s.pdf', num2str(corrType))
print(fn, '-dpdf')



% Run regression to do stats on performance:
% Regress accuracy onto 

accMat=[]; expMat=[];

for k =1:size(allAccuracy,3)
    % Build a performance matrix:
    accMat=[accMat; allAccuracy(:,:,k)];
    
    % And an explanatory matrix:
    if corrType==1
        newExp=[ones(size(allAccuracy(:,:,k),1),1), ones(size(allAccuracy(:,:,k),1),1).*allCorrCombos.inPool(k), ones(size(allAccuracy(:,:,k),1),1).*allCorrCombos.inRel(k)];
        legendText={'no correlation', 'in pool', 'in pool + rel pool'};
    elseif corrType==2
        newExp=[ones(size(allAccuracy(:,:,k),1),1), ones(size(allAccuracy(:,:,k),1),1).*allCorrCombos.inPool(k), ones(size(allAccuracy(:,:,k),1),1).*allCorrCombos.inIrrel(k)];
        legendText={'no correlation', 'in pool', 'in pool + irrel pool'};
    end
    expMat=[expMat; newExp];
end



bins=0:10:100;
for i = 1:(length(bins)-1)
   
     meanPerf(i,1)=nanmean(nanmean(accMat(expMat(:,2)==0&expMat(:,3)==0, (bins(i)+1):bins(i+1) )));
     meanPerf(i,2)=nanmean(nanmean(accMat(expMat(:,2)==.2&expMat(:,3)==0, (bins(i)+1):bins(i+1) )));
     meanPerf(i,3)=nanmean(nanmean(accMat(expMat(:,2)==.2&expMat(:,3)==.2, (bins(i)+1):bins(i+1) )));
     
     
     tmp1=accMat(expMat(:,2)==0&expMat(:,3)==0, (bins(i)+1):bins(i+1));
     tmp2=accMat(expMat(:,2)==.2&expMat(:,3)==0, (bins(i)+1):bins(i+1));
     tmp3=accMat(expMat(:,2)==.2&expMat(:,3)==.2, (bins(i)+1):bins(i+1));
     
     
     semPerf(i,1)=nanstd(tmp1(:))./sqrt(nReps);     
     semPerf(i,2)=nanstd(tmp2(:))./sqrt(nReps);
     semPerf(i,3)=nanstd(tmp3(:))./sqrt(nReps);
end



close all

defaultPlotParameters
hold on
H1=shadedErrorBar(bins(2:end),meanPerf(:,1),semPerf(:,1),{'color', cbColors(2,:)})
H2=shadedErrorBar(bins(2:end),meanPerf(:,2),semPerf(:,2),{'color', cbColors(3,:)})
H3=shadedErrorBar(bins(2:end),meanPerf(:,3),semPerf(:,3),{'color', cbColors(4,:)})

ff=legend([H1.mainLine, H2.mainLine, H3.mainLine], legendText)
set(ff, 'location', 'east',  'box', 'off')
ylabel('Percent Correct')
xlabel('Trials')
xlim([10,100])
set(gca, 'box', 'off')
fn=sprintf('learningTimecourse_%s_%s.pdf', num2str(corrType), date)
print(fn, '-dpdf')




hold on
a=plot(nanmean(accMat(expMat(:,2)==0&expMat(:,3)==0, :)), 'b')
b=plot(nanmean(accMat(expMat(:,2)==.2&expMat(:,3)==0, :)), 'g')
c=plot(nanmean(accMat(expMat(:,2)==.2&expMat(:,3)==.2, :)), 'r')


figure
imagesc(nanmean(meanDistMat, 3),[.5,2])
set(gca, 'xtick', 1:6, 'xticklabel', label, 'ytick', 1:6, 'yticklabel', label) 
title('Distance from optimal readout')
xlabel('irrelevant pool correlations')
ylabel('In pool correlations')
colorbar
fn=sprintf('weightDistHeatmap_corrType_%s.pdf', num2str(corrType));
print(fn, '-dpdf')
close all


fn=sprintf('abstractCorrsWorkspace_corrType%s_%s.mat', num2str(corrType), date);
save(fn)

end



input('pause here')




%% Load both workspaces and make key figures:

% RESET DATE HERE if you rerun simulation!!!
% relCorrWorkspace='abstractCorrsWorkspace_corrType1_27-Mar-2021.mat'
% irrelCorrWorkspace='abstractCorrsWorkspace_corrType2_27-Mar-2021.mat'
fn1=sprintf('abstractCorrsWorkspace_corrType1_%s.mat',  date);
relCorrWorkspace=fn1;
fn2=sprintf('abstractCorrsWorkspace_corrType2_%s.mat',  date);
irrelCorrWorkspace=fn2;

% load both relevant and irrelevant correlation workspaces:
relCorrSim=load(relCorrWorkspace);
irrelCorrSim=load(irrelCorrWorkspace);

% Make figure based on both:
makeTwoDimFig






%% Get statistics for paper:

noCorrSim=relCorrSim.allCorrCombos.inPool==0 &relCorrSim.allCorrCombos.inRel==0;
inPoolCorSim=relCorrSim.allCorrCombos.inPool==.2 &relCorrSim.allCorrCombos.inRel==0;

noCorrAcc=((mean(relCorrSim.allAccuracy(:,:,noCorrSim),2)))
mean(noCorrAcc)
std(noCorrAcc)

inCorrAcc=((mean(relCorrSim.allAccuracy(:,:,inPoolCorSim), 2)))
mean(inCorrAcc)
std(inCorrAcc)

[H,P,CI,STATS]=ttest2(inCorrAcc, noCorrAcc)

% Now get stats for adding relevant pool correlations:
inPoolCorrSim=relCorrSim.allCorrCombos.inPool==.2 &relCorrSim.allCorrCombos.inRel==0;
relPoolCorSim=relCorrSim.allCorrCombos.inPool==.2 &relCorrSim.allCorrCombos.inRel==.2;

inCorrAcc=((mean(relCorrSim.allAccuracy(:,:,inPoolCorrSim),2)))
mean(inCorrAcc)
std(inCorrAcc)

relCorrAcc=((mean(relCorrSim.allAccuracy(:,:,relPoolCorSim), 2)))
mean(relCorrAcc)
std(relCorrAcc)

[H,P,CI,STATS]=ttest2(relCorrAcc, inCorrAcc)

% Now get stats for adding irrelevant pool correlations:
inPoolCorrSim=irrelCorrSim.allCorrCombos.inPool==.2 &irrelCorrSim.allCorrCombos.inIrrel==0;
relPoolCorSim=irrelCorrSim.allCorrCombos.inPool==.2 &irrelCorrSim.allCorrCombos.inIrrel==.2;

inCorrAcc=((mean(irrelCorrSim.allAccuracy(:,:,inPoolCorrSim),2)))
mean(inCorrAcc)
std(inCorrAcc)

irrelCorrAcc=((mean(irrelCorrSim.allAccuracy(:,:,relPoolCorSim), 2)))
mean(irrelCorrAcc)
std(irrelCorrAcc)

[H,P,CI,STATS]=ttest2(irrelCorrAcc, inCorrAcc)


%% MAke supplementary figure showing that optimal readout leads to same performance:



