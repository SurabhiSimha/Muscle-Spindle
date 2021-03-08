function J = spindle_cost(gains,data,flags)
% Created: 5/2013
% Modified: 5/2016
% Author: Kyle Blum
% Description: This function is the objective function used by 'fmincon.m'
% in 'spindle_findGains.m' to fit kinematics/force data to spindle IFRs.
% Returns cost, J, of current fit which is calculated as the weighted sum
% of squared errors and maximum error between the fitted data and the
% actual IFR.

%%% Choose model to use %%%
switch flags.model
    case 1 
        FR_sim = kinetics(data,gains,flags);
    case 2
        FR_sim = kinematics(data,gains,flags);
%     case 3 
%         FR_sim = kinematics(data,gains,flags);
    case 3
        FR_sim = prochazka(data,gains);
    case 4
        FR_sim = kinematics(data,gains,flags);
    case 5
        FR_sim = kinematics(data,gains,flags);
    case 6
        FR_sim = mixedKin(data,gains);
end

%%% Compute errors between recorded data and fit %%%
FR_recorded = data.firing_rate;
FR_error = FR_recorded - FR_sim;

if isempty(data.firing_rate);
    J = 0;
    return
end

%%% Obective Function %%%
% j1 = flags.weights(1);
% j2 = flags.weights(2);

SSE = sum(FR_error.^2);
SSM = sum((FR_recorded-mean(FR_recorded)).^2);


% L1 = j1*sum((FR_error.*FR_error));
% L2 = j2*max(abs(FR_error));

J = SSE/SSM;