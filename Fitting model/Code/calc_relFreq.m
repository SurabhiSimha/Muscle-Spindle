% Created 08/01/2016
% Author: KPB
% This function calculates the relative frequencies that a given model was
% the best from the set of models tested given their information criteria
% from the cross-validation. 


function relFreq = calc_relFreq(xvaldata)

IC_array = zeros(6,100);
relFreq = zeros(1,6);

for j = 1:6
    IC_array(j,:) = xvaldata.MLEout(j).test_AICc;
end

[~,i] = min(IC_array);

for k = 1:6
relFreq(k) = sum(i == k)/size(IC_array,2);
end
    
end