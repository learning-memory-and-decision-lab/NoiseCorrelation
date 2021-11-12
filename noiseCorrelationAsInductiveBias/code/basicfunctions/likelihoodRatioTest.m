function [p, devDiff]=likelihoodRatioTest(dev_full, dev_base, degFree)
    % likelihood ratio test for logistic regression models.
    % dev_full is the deviance of the full (larger) model
    % dev_base is the deviance of the restricted (base) model.
    % degFRee is the difference in number of parameters.
    
    
devDiff=  dev_base-dev_full;
prob=chi2cdf(devDiff, degFree);
p=1-prob;