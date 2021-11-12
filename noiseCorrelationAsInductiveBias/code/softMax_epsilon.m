function [pChoice choice]=softMax_epsilon(values, invT, epsilon)

% epilon softmax action selection:

values          =  values-(max(values));
pChoice_softMax =  exp(values.*invT)./nansum(exp(values.*invT));
pChoice         =  pChoice_softMax.* (1-epsilon) + epsilon.*ones(size(values))./length(values);


if nargout>1
    choice=find(cumsum(pChoice)>rand, 1);
end

