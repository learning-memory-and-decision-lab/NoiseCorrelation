% Make schematic figure

% This script just creates a basic scheumatic that will be the first figure
% in the noise correlations as inductive biases paper. 

totVar = 1;
nNeuronsPerPool=2;
all_fracCorrNoise=[0, .5]
tags={'nc00', 'nc50'}

nPoints=10;
poolID=ones(nNeuronsPerPool,1)



for k = 1:2




fracCorrNoise=all_fracCorrNoise(k);

% Choose an overall level of neural variance based on equation for fixed
% lambda and tot variance:
frVar = totVar ./ (nNeuronsPerPool + nNeuronsPerPool.*(nNeuronsPerPool-1).*fracCorrNoise);

inPoolCov       = frVar.*fracCorrNoise; % some correlations within pools.
outPoolCov      = 0;  % no correlations across pools.
nPools          = 1;





% create a matrix stipulating which neurons are in same pool:
samePool=repmat(poolID, 1, length(poolID))==repmat(poolID, 1, length(poolID))';
sameNeuron=logical(eye(length(poolID)));

targFR=1;
nonTargFR=0;

% create a covariance matrix across neural population:
covMat=nan(size(sameNeuron));
covMat(sameNeuron)=frVar;
covMat(samePool&~sameNeuron)=inPoolCov;
covMat(~samePool&~sameNeuron)=outPoolCov;


allStims=[ones(nPoints, 1);  zeros(nPoints, 1)]

% preallocate space for firing rates:
%firingRate(:,i)=nan(nNeuronsPerPool.*nPools, 1);


clear firing rate

for i = 1:length(allStims)


Mu(poolID==allStims(i))=targFR;
Mu(poolID~=allStims(i))=nonTargFR;
FR = mvnrnd(Mu,covMat);
firingRate(:,i)=FR;


accuracy(i)=(mean(firingRate(:,i))>.5 )   ==allStims(i);


end

allCovMat(:,:,k)=covMat;
allFR(:,:,k)=firingRate;

end



% This should be panel 2
figure
hold on
%plot([-boundSize, boundSize]+.5, [boundSize, -boundSize]+.5, '-k');

c=error_ellipse(allCovMat(1:2,1:2,1),[targFR, targFR], 'conf', .85)
d=error_ellipse(allCovMat(1:2,1:2,1),[nonTargFR, nonTargFR], 'conf', .85)
set(c, 'color', 'b')
set(d, 'color', 'r')
e=error_ellipse(allCovMat(1:2,1:2,2),[targFR, targFR], 'conf', .85)
f=error_ellipse(allCovMat(1:2,1:2,2),[nonTargFR, nonTargFR], 'conf', .85)
set(e, 'color', 'b', 'LineStyle', '--')
set(f, 'color', 'r', 'LineStyle', '--')
set(gca, 'box', 'off', 'xticklabel', '', 'yticklabel', '');
saveas(gcf,'fixSigToNoise_100Units.eps', 'epsc2')
close all
  


figure
hold on

a=0;
%plot([-boundSize, boundSize]+.5, [boundSize, -boundSize]+.5, '-k'); 
C=error_ellipse([1, a; a, 1] ,[targFR, targFR], 'conf', .85)
D=error_ellipse([1, a; a, 1],[nonTargFR, nonTargFR], 'conf', .85)
set(C, 'color', 'b')
set(D, 'color', 'r')
hold on
a=.5;
E=error_ellipse([1, a; a, 1],[targFR, targFR], 'conf', .85)
F=error_ellipse([1, a; a, 1],[nonTargFR, nonTargFR], 'conf', .85)
set(E, 'color', 'b', 'LineStyle', '--')
set(F, 'color', 'r', 'LineStyle', '--')
set(gca, 'box', 'off', 'xticklabel', '', 'yticklabel', '');
saveas(gcf,'theStandardStory.eps', 'epsc2')



NC_makeFig1;








