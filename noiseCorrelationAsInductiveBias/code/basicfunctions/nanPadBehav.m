function padBehav=nanPadBehav(behav)

%% Pad all fields of an array that are "Too short" so that they contain
%% nans and are the same size as other fields.



allNames=fieldnames(behav);



%eval(sprintf('maxLn=length(behav.%s)', char(allNames(1))))


for i = 1:length(allNames)
        % get data from a field
        ln(i)=eval(sprintf('length(behav.%s)', char(allNames(i))));

end

maxLn=max(ln);

for i = 1:length(allNames)
       
    
    if ln(i)<maxLn
    
    eval(sprintf('behav.%s(end+1:maxLn,:)=nan', char(allNames(i))));
    end
  
    
end


padBehav=behav;
