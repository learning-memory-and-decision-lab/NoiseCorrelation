% testing theoretical derivation of vector summation to estimate
% SEM of correlated variables:


Theta=30;
sigmaObs=1;
nObs = 20;
nSim=1000;


% create covariance matrix:
Sigma=eye(nObs);
Sigma(Sigma==0)=cosd(Theta);
Mu=zeros(nObs, 1);


genData=mvnrnd(Mu,Sigma, nSim)
empiricalSEM=nanstd(mean(genData, 2))



% TRIANGLE METHOD: THIS ONLY WORKS FOR TWO SAMPLES!!!!
% OTHERWISE IT IS ONLY INCLUDING 1 set row of covariance matrix (not all
% entries):
% theoreticalSEM=sigmaObs.*sqrt(cosd(Theta)+1)./sqrt(nObs)



% STATS METHOD from here: 
% https://stats.stackexchange.com/questions/44032/what-is-the-standard-deviation-of-the-sum-of-three-correlated-random-variables

%Var_tot = Var_x1 + Var_x2 + Var_x3 ... + 2Var_cov(x1, x2) +
            %2Var_cov(x1,x3) + 2Var_cov(x2,x3)

            
            
            
covParts=tril(true(nObs),-1);
var_tot=sum(Sigma(logical(eye(nObs))))+2.*sum(Sigma(covParts));


length(Sigma).*(length(Sigma)-1)


% sigma_{pop}=SigmaSq_i*n + CovSq*n*(n-1)


% sqrt of total variance divided by the number of observations:
theoreticalSEM2= sqrt(var_tot)./nObs



% LETS CREATE A BUNCH OF 






a=tril(ones(10), -1)
sum(a(:))







%% OK, what about a situation where we have two types of correlations 
% (analagous to the relavant pool simulation).

type  =repmat(1:4, 1, 4)';
relDim=mod(type, 2);

sameNeuron=logical(eye(16))
samePool  =repmat(type, 1, length(type))== repmat(type', length(type), 1);
sameRelDim=repmat(relDim, 1, length(type))== repmat(relDim', length(type), 1);


inPoolVar=.2
inRelVar =-.2;

Sigma=nan(size(sameNeuron));
Sigma(sameNeuron)=1;
Sigma(samePool&~sameNeuron)=inPoolVar;
Sigma(sameRelDim&~samePool)=inRelVar;
Sigma(~sameRelDim&~samePool)=0;

Mu=zeros(length(Sigma), 1);
genData=mvnrnd(Mu,Sigma, nSim);
empiricalSEM=nanstd(mean(genData, 2))

var_tot=sum(sameNeuron(:))+sum(samePool(:)&~sameNeuron(:)).*inPoolVar + sum(sameRelDim(:)&~samePool(:)).*inRelVar
theoreticalSEM2= sqrt(var_tot)./length(Sigma)


