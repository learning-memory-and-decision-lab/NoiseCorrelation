
%% RUN THIS SECTION FOR EITHER SIMULATION?

whichComp       =  1;
targetSEM       =  1;    % standard error on mean population representation
nPools          =  2;
nNeuronsPerPool = 100;
targFR          =  1;
nonTargFR       = -1; %
%LR              =.0001;  % Now doing a bunch of LRs
doRL            =true;
weightScale     =1;
useRPE          =false; % online predictions
doOpt           =true;
sigToNoise      = 2;
all_fracCorrNoise=0:.01:.2;
fixNoise        =false;


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

% Target values currently only make sense for 2 pools...
targetOutputs=[100; -100;]; % choose values for decision neuron output

% create a list of neuron pool assignments:
poolID=[];
for i = 1:nPools
    poolID(end+1:end+nNeuronsPerPool,1)=ones(nNeuronsPerPool,1).*i;
end
meanFR=(targFR+nonTargFR.*(nPools-1))./nPools;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Network #1  -- two layers, simplest possible proof of concept:
%          Run multiple sims to get "averaged" performance:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE -- initial weights matter quite a bit for the "weight profile" plot

nCorrs=length(all_fracCorrNoise);


if ~doRL
    numOutNeurons = 1;
else
    numOutNeurons = nPools;
end

% Compute total signal and choose a total amount of "decision relevant noise" accordingly:
totSig=(targFR-nonTargFR).*nNeuronsPerPool;
totVar=(totSig./sigToNoise)^2;


% create a list of neuron pool assignments:
poolID=[];
for i = 1:nPools
    poolID(end+1:end+nNeuronsPerPool,1)=ones(nNeuronsPerPool,1).*i;
end
meanFR=(targFR+nonTargFR.*(nPools-1))./nPools;


clear finalWeights

alLRs=[.001,.002, .005 .01];

for L =1:length(alLRs);
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
            targFR           = targFR.*scaleFactor;
            nonTargFR        = nonTargFR.*scaleFactor;
        end
        
        % create a matrix stipulating which neurons are in same pool:
        samePool=repmat(poolID, 1, length(poolID))==repmat(poolID, 1, length(poolID))';
        sameNeuron=logical(eye(length(poolID)));
        
        
        % create a covariance matrix across neural population:
        covMat=nan(size(sameNeuron));
        covMat(sameNeuron)=frVar;
        covMat(samePool&~sameNeuron)=inPoolCov;
        covMat(~samePool&~sameNeuron)=outPoolCov;
        
        %% compute learning speed analytically based on EXPECTED gradients:
        % ala Jan D code -- this doesn't make sense for us -- since the
        % expected gradients are not very typical in the face of noisy
        % inputs. 
        
        Mu(poolID==1)=targFR;
        Mu(poolID~=1)=nonTargFR;
        wStar=[ones(nNeuronsPerPool, 1); ones(nNeuronsPerPool, 1).*-1]; % scale could be off here -- should do this analytically too
        wN=zeros(nNeuronsPerPool.*2,1);                                 % start from nothing -- could do this from any position though.
        wStarMinusWn=wStar-wN;
        learningSpeed(k)=wStarMinusWn'*covMat*wStarMinusWn + (wStarMinusWn'*Mu')^2;%
        
        
        %     %% compute learning speed for W=zeros, but considering full distribution over X:
        %
        %     hold on
        %     possSpeeds=0:.001:.1
        %      eq6_Mu= 2.*LR.*Mu*wStarMinusWn;
        %     eq6_Cov=(2.*LR.*wStarMinusWn)'*covMat*(2.*LR.*wStarMinusWn);
        %     learnSpeed_fullDist=normpdf(possSpeeds,eq6_Mu,sqrt(eq6_Cov)) + LR^2.*Mu*Mu';
        %     plot(possSpeeds, learnSpeed_fullDist);
        %
        
        %% compute learning speed according to full distribution as described here:
        %  https://math.stackexchange.com/questions/442472/sum-of-squares-of-dependent-gaussian-random-variables/442916#442916
        %  Q(x) = X^{t}AX = \sum_{j=1}^n \lambda_j(U_j+b_j)^2
        %  OK, in order to get P, we need to find eigenvectors of covariance
        %  matrix:
        %  covMat^.5*eye(length(covMat))*covMat^.5 - covMat
 
        
        %testing:
        %covMat=[1, 1; 1, 1];
        
        
        %% Validate eigenvectors/eigenvalues:
        %     a=covMat*v
        %     b=v*d
        % Get variables necessary to compute euclidean norm:
        LR=alLRs(L);
        meanError=(wN+LR*Mu')-wStar;
        covError=LR*covMat*LR'; % covariance of firing rates scaled by learning rate (which is multiplied to control update)
        
        % Let Y be cov^-.5 *X -- so Y is non-central multivariate normal
        %                       with equal identity variance elements
        
        % Let Z = (Y-cov^.5 * /mu) -- so z is a mean zero multivariate
        %                       normal with equal variances, but some
        %                       non-zero covariance
        
        % we will get P and L where P is orthonormal projectsion such that PLP=Cov 

        [P, d]=eig(covError); % get eigenvectors and define as P
        lambda=eig(covError); % get eigenvalues
        
        % covError = P*d*P'   -- so P is transposed with respect to
        %                        stack exchange solution

        b=P'*(covError^-.5)   *(meanError);
        
        % for now, lets just take the sum of the expectations:
        % non-central chi squared distribution has mean:
        % k + lambda, where lambda is the noncentrality parameter
        % (equal to sum of the squared mean vector) equal to
        % the sum of the squared mean difference vector
        
        clear expectNorm_eig
        for i =1:length(lambda)
            
            % OK -- solution is weighted chi squared distribution:
            % Mean of chi square is sum of squared means (b(i)^2) plus the
            % number of degrees of freedom -- which is always 1 in our
            % case. 
 
            expectNorm_eig(i)=lambda(i).* (1+ b(i).^2);

        end
        expectNorm(k)=sum(expectNorm_eig);
        initError=wStar'*wStar;
    end
    
    allExpNorms(:,L)=expectNorm;
    
end

learningSpeed2=initError-allExpNorms;
a=plot(repmat(all_fracCorrNoise, 4, 1)', learningSpeed2)
ff=legend(a, 'LR = .001', 'LR = .002', 'LR = .005', 'LR = .01')
set(ff, 'location', 'northeast', 'box', 'off')
set(gca, 'box', 'off')
ylabel('Learning speed')
xlabel('Noise correlations')
saveas(gcf, 'analyticalLearningSpeedPlot.eps', 'epsc2')
close all



