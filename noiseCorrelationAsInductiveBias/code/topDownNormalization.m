
%% Set things up:

% create covariance matrix, mu, and weight matrix:
covMat=eye(200);
Mu=[ones(100,1);-ones(100,1)];

weightScale=1;

% Create a weight matrix 
%wtMatrix=[[ones(100,1);-ones(100,1)]   [-ones(100,1);ones(100,1)]]; % ideal weight matrix
%wtMatrix=normrnd(0, 1, 200, 2); % random weight matrix. 

wtMatrix=[[ones(100,1);-ones(100,1)]   [-ones(100,1);ones(100,1)]] + normrnd(0, 2, 200, 2); 




wtDiff=wtMatrix(:,1)-wtMatrix(:,2);
wtDiffNorm=wtDiff./norm(wtDiff);


% Normalize weight matrix:
normWtMatrix=wtMatrix./norm(wtMatrix);


% normWtMatrix(:,1)=wtMatrix(:,1)./(sum(wtMatrix(:,1).^2 ).^0.5)
% normWtMatrix(:,2)=wtMatrix(:,2)./(sum(wtMatrix(:,2).^2 ).^0.5)



% initialize decision neuron
decNeuronOutput=zeros(2,1)';

%% Create signal and take one pass through loop:

% Loop through timesteps within a single trial
FR=mvnrnd(Mu,covMat)+decNeuronOutput*normWtMatrix'; % pass back decision activity to input units
% compute activation of downstream neuron(s):
decNeuronOutput=FR*normWtMatrix;
  

%% take a second pass through the loop -- but now don't add signal
%  The idea is that we should get back the same decision neuron activity we
%  put in... 


FR2 =   decNeuronOutput*normWtMatrix';
decNeuronOutput2=FR2*normWtMatrix;



%% Plot initial versus end:

figure
hold on
minMax=[min([decNeuronOutput, decNeuronOutput2]), max([decNeuronOutput, decNeuronOutput2])];
plot(decNeuronOutput, decNeuronOutput2, 'o')
plot(minMax, minMax, '--r')

ylabel('Decision neuron: round 2')
xlabel('Decision neuron: round 1')


figure
hold on
plot(FR, 'b')
plot(FR2, 'r')




%% can simulation above create noise correlations? 
% Lets loop through trials and find out. 

%% Set things up:

% create covariance matrix, mu, and weight matrix:
covMat=eye(200);
Mu=[ones(100,1);-ones(100,1)];

zeroCoherenceMu=zeros(length(Mu), 1); % mean firing for "zero coherence" trials.


% Create a weight matrix 
%wtMatrix=[[ones(100,1);-ones(100,1)]   [-ones(100,1);ones(100,1)]]; % ideal weight matrix
%wtMatrix=normrnd(0, 1, 200, 2); % random weight matrix. 

% lets make a weight matrix that contains some signal info and some noise:
wtMatrix=[[ones(100,1);-ones(100,1)]   [-ones(100,1);ones(100,1)]] + normrnd(0, 2, 200, 2); 

% Normalize weight matrix:
normWtMatrix=wtMatrix./norm(wtMatrix);


% normWtMatrix(:,1)=wtMatrix(:,1)./(sum(wtMatrix(:,1).^2 ).^0.5)
% normWtMatrix(:,2)=wtMatrix(:,2)./(sum(wtMatrix(:,2).^2 ).^0.5)



%%
% loop through 100 zero coherence trials and see whether noise correlations
% emerge:
for i = 1:100


% initialize decision neuron
decNeuronOutput=zeros(2,1)';

% Create signal and take one pass through loop:

% Loop through timesteps within a single trial -- now only zero coherence
% trials -- so there is no signal difference across units or trials:
FR=mvnrnd(zeroCoherenceMu,covMat)+decNeuronOutput*normWtMatrix'; % pass back decision activity to input units
% compute activation of downstream neuron(s):
decNeuronOutput=FR*normWtMatrix;
  

%% take a second pass through the loop -- but now don't add signal
%  The idea is that we should get back the same decision neuron activity we
%  put in... 

FR2 = decNeuronOutput*normWtMatrix';
decNeuronOutput2=FR2*normWtMatrix;


store_FR1(i,:)=FR;
store_FR2(i,:)=FR2;
end

figure 
imagesc(corr(store_FR1), [-1, 1])
title('noise correlations: initial input')


figure 
imagesc(corr(store_FR2), [-1, 1])
title('noise correlations: through top down')


