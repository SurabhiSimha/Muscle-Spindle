function J = cost(gains,data,flags)
% Created: 5/2013
% Modified: 5/2016
% Author: Kyle Blum
% Description: This function is the objective function used by 'fmincon.m'
% in 'spindle_findGains.m' to fit kinematics/force data to spindle IFRs.
% Returns the cost, J, of current fit which is calculated as the sum of 
% squared errors fitted data and the measured IFR.

%%% Choose model to use %%%
switch flags.model
    case 1 
        [FR_sim,~,~] = kinetics(data,gains,flags); %Get estimated IFR from force model
    case 2
        FR_sim = kinematics(data,gains,flags);     %Get estimated IFR from length model
end

%%% Compute errors between recorded data and fit %%%
FR_recorded = data.IFR;
FR_error = FR_recorded - FR_sim;

if isempty(data.IFR) %If there are no IFR values (unlikely)
    J = 0;
    return
end

%%% Obective Function %%%

SSE = sum(FR_error.^2); %Calculate sum of squared errors
SSM = sum((FR_recorded-mean(FR_recorded)).^2); %Calculate sum of squares about the mean

J = SSE/SSM; %Calculate cost