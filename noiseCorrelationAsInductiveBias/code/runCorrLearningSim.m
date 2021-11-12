function output=runCorrLearningSim(input)


% unpack inputs:
sigToNoise       = input.sigToNoise; 
%targetSEM       = input.targetSEM;    % standard error on mean population representation
nPools           = input.nPools;
nNeuronsPerPool  = input.nNeuronsPerPool;
targFR           = input.targFR;
nonTargFR        = input.nonTargFR;
nTrials          = input.nTrials;
LR               = input.LR;
doRL             = input.doRL;
invT             = input.invT;
weightScale      = input.weightScale;
useRPE           = input.useRPE; % this produces weird results -- leave off for now.
doOpt            = input.doOpt;
all_fracCorrNoise= input.all_fracCorrNoise;


% MRN added 3/1/21 to address reviewer concerns (eg. fix noise to show same
% effects)
if ~isfield(input, 'fixNoise')
fixNoise=false;
else
fixNoise =input.fixNoise;
disp('Fixing noise and manipulating signal according to correlations')
end




% create output structure:
output = struct;


% Set key variables according to inputs:

if ~doRL
numOutNeurons = 1;
else
numOutNeurons = nPools;
end


% SNR = signal / sigma

% fix signal to noise ratio -- note this is just for one pool, population
% SNR is computed below for 2 pools. 
totSig= (input.targFR-input.nonTargFR).*nNeuronsPerPool;
totVar=     (totSig./input.sigToNoise)^2;

% MRN trying to match SNR between my code and dan's code:
% recompute SNR based on true signal and variance:
trueSignal=totSig.*2; % twice the signal (2 pools)
trueSigma =sqrt(totVar+totVar); % twice the variance in the signal direction
trueSNR   =trueSignal./trueSigma; % numerator doubles, denominator increases by root 2
% This is SNR for entire population response (not just 1 pool). 


%targetSEM=
%totVar = nNeuronsPerPool.*targetSEM.^2; % FIXED to a target level!
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
    optWeights=[];
    for pp=1:max(poolID)
    optWeights(poolID==pp,1:(pp-1))=nonTargFR;
    optWeights(poolID==pp,pp)=targFR;
    optWeights(poolID==pp,(pp+1):max(poolID))=nonTargFR;
    end
end

clear finalWeights
for k = 1:length(all_fracCorrNoise)
    fracCorrNoise=all_fracCorrNoise(k);
    
    
    if ~fixNoise
    % Choose an overall level of neural variance based on equation for fixed
    % lambda and tot variance:
    frVar = totVar ./ (nNeuronsPerPool + nNeuronsPerPool.*(nNeuronsPerPool-1).*fracCorrNoise);
    inPoolCov       = frVar.*fracCorrNoise; % some correlations within pools.
    outPoolCov      = 0;  % no correlations across pools.
    else
    frVar           = totVar ./ (nNeuronsPerPool); % Neurons have same variability in firing... 
    inPoolCov       = frVar.*fracCorrNoise; % some correlations within pools.
    outPoolCov      = 0;  % no correlations across pools.
    % But adjust target firing rate according to noise correlation
    % condition -- see multSignalAndNoiseMethods.docx:
    S_neuron_base=sqrt(frVar./nNeuronsPerPool);
    S_neuron_NC  =sqrt(frVar.*(1+(nNeuronsPerPool-1).*fracCorrNoise)./nNeuronsPerPool);
    scaleFactor=S_neuron_NC./S_neuron_base;
    targFR           = input.targFR.*scaleFactor;
    nonTargFR        = input.nonTargFR.*scaleFactor;
    end

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
    % create a random weight matrix that projects onto a single "decision" neuron:
    wtMatrix=normrnd(0, .00001, length(poolID),numOutNeurons);
%     if doRL % normalize squared weighgts for both output neurons:
%       wtMatrix(:,1)=wtMatrix(:,1)./(sum(wtMatrix(:,1).^2)./weightScale);    
%       wtMatrix(:,2)=wtMatrix(:,2)./(sum(wtMatrix(:,2).^2)./weightScale);    
%     end
    
    
    
    % Loop through trials and simulate firing rates:
    Mu=nan(size(poolID));
    for i = 1:nTrials
        Mu(poolID==allStims(i))=targFR;
        Mu(poolID~=allStims(i))=nonTargFR;
        FR = mvnrnd(Mu,covMat);
        firingRate(:,i)=FR;
        
        % compute activation of downstream neuron(s):
        decNeuronOutput=FR*wtMatrix;

        
        % make a choice:
        if ~doRL
            [pChoice, choice(i)]=softMax_epsilon([decNeuronOutput, 0], invT, 0);
            expRew=pChoice(choice(i));
        else
            [pChoice, choice(i)]=softMax_epsilon([decNeuronOutput], invT, 0);
            if useRPE
            expRew=pChoice(choice(i));
            else
            expRew=.5;
            end
            
        end
        
        
        if any(~isfinite(decNeuronOutput))
            disp('houston, we have a problem!')
            keyboard
        end
        
        % If you already knew the best decision rule:
        if doOpt
        optOutput=FR*optWeights;
        end
        
        if ~doRL
        % compute error:
          PE=targetOutputs(allStims(i))-decNeuronOutput;
        % update weights:
          wtMatrix=wtMatrix+ LR.*FR'.*PE;
          storePE(i)=PE;
        else
             % deltaW(C) = alpha * RPE * [X-E(X)]
             RPE=(choice(i)==allStims(i)) - expRew; % compute RPE
             deltaWt=LR.*RPE.*(FR'-meanFR); % compute weight update
             wtMatrix(:,choice(i))= wtMatrix(:,choice(i))+deltaWt; %implement weight update
             storePE(i)=RPE;  % store RPE
             % MRN ditching normalization: 3/27/21
             % wtMatrix(:,choice(i))=wtMatrix(:,choice(i))./(sum(wtMatrix(:,choice(i)).^2 ).^0.5./weightScale); % normalize weights

        end

        accuracy(i,k)=choice(i)==allStims(i);
        
        if doOpt
            if ~doRL
            optAccuracy(i,k)=sign(optOutput)==sign(targetOutputs(allStims(i)));
            else
                [val,ind]=max(optOutput);
                optAccuracy(i,k)=allStims(i)==ind;
            end
        end
    end
    
    
    if ~doRL
        finalWeights(:,k)=wtMatrix;
    else
        finalWeights(:,:,k)=wtMatrix;
    end
    % Can we compute how: it would be good to estimate

end

% spit out the things we care about:
output.finalWeights=finalWeights;
output.optAccuracy =optAccuracy;
output.accuracy    =accuracy; 

