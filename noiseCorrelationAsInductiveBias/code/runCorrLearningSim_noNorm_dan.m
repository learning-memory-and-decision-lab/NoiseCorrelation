function output=runCorrLearningSim_noNorm_dan(input)


% this is script that Dan was using to generate middle panel in his
% analytical figure. 



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



totSig=(input.targFR-input.nonTargFR).*nNeuronsPerPool;
totVar=     (totSig./input.sigToNoise)^2;

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

%    load('1-step-wgts.mat')
%     if doRL % normalize squared weighgts for both output neurons:
%       wtMatrix(:,1)=wtMatrix(:,1)./(sum(wtMatrix(:,1).^2)./weightScale);    
%       wtMatrix(:,2)=wtMatrix(:,2)./(sum(wtMatrix(:,2).^2)./weightScale);    
%     end
    
   %%
   % Pools
   n = 200;
   m1 = [ones(100,1) ; zeros(100,1)];
   m2 = [zeros(100,1);  ones(100,1)];
   m1h = m1./norm(m1);
   m2h = m2./norm(m2);
   
   % Signal dimension
   v_s = (m1 - m2)./norm(m1 - m2);
   a_s = (m1 + m2)./norm(m1 + m2);

   % Get a noise representative
   xi = mvnrnd(zeros(n,1), eye(n))';
   
   % Orthogonalize it to the signal
   %xi = xi - (m1h-m2h)*xi'*(m1h - m2h);
   
   % Since the projection onto ones(200,1) really matters
   % orthogonalize it to this as well.
   %xi = xi - (m1h + m2h)*xi'*(m1h + m2h);
   
   xi = xi - m1h*xi'*m1h - m2h*xi'*m2h;

   
   % Normalize
   xi = xi./norm(xi);
   
   [v,d] = eig(covMat);
   d = diag(d);
   [d,inds] = sort(d,'descend');
   v = v(:,inds);
   %%


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
            % wtMatrix(:,choice(i))=wtMatrix(:,choice(i))./(sum(wtMatrix(:,choice(i)).^(0.5))./weightScale); % normalize weights
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
        
      % 1D weights
      w  = [1,-1]*wtMatrix';
      dw = deltaWt;

      % "Actual" firing rate (sign associated with stimulus removed)
      aFR = (-1)^mod(allStims(i)+1,2)*FR;
      
      % Firing rate component norms
      FRn(k,i)   = norm(aFR);
      FRpn(k,i)  = norm(aFR - aFR*v_s*v_s');
      FRfn(k,i)  = norm(aFR - aFR*v_s*v_s' - aFR*a_s*a_s');

      % Observed weight components over time
      wsn(k,i)  = w*v_s;
      wpn(k,i)  = norm(w - w*v_s*v_s');
      wfn(k,i)  = norm(w - w*v_s*v_s' - w*a_s*a_s');

      % Observed component increments
      dwn(k,i)  = norm(dw);
      dwpn(k,i) = norm(dw' - dw'*v_s*v_s');
      dwfn(k,i) = norm(dw' - dw'*v_s*v_s' - dw'*a_s*a_s');

      % Noise transfer into signal space
      ns(k,i)  = v_s'*aFR';
      nxf(k,i) = (w - w*v_s*v_s')*aFR';
    end
    
    
    
%    subplot(1,3,2)
%    plot(nxf(k,:))
%    
%    title('Noise Transfer')
%    xlabel('Trial')
%    ylabel('')
%    grid on
   
   
    if ~doRL
        finalWeights(:,k)=wtMatrix;
    else
        finalWeights(:,:,k)=wtMatrix;
    end
    % Can we compute how: it would be good to estimate

    allCovMat(:,:,k)=covMat;
end

% spit out the things we care about:
output.finalWeights=finalWeights;
output.optAccuracy =optAccuracy;
output.accuracy    =accuracy; 

% Spit out stuff needed for Dan's plot
output.wsn=wsn;
output.wpn=wpn;
output.allCovMat=allCovMat;
output.Mu = Mu;





 

