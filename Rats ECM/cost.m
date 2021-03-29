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
        [FR_sim,~,~,~,~] = kinetics2(data,gains,flags);     
end

%%% Compute errors between recorded data and fit %%%
FR_recorded = data.IFR;
FR_error = FR_recorded - FR_sim;

if isempty(data.IFR) %If there are no IFR values (unlikely)
    J = 0;
    return
end

%%% Obective Function %%%
switch flags.cost
    case 1
        %R2
        SSE = sum(FR_error.^2); %Calculate sum of squared errors
        SSM = sum((FR_recorded-mean(FR_recorded)).^2); %Calculate sum of squares about the mean
        J = SSE/SSM; %Calculate cost
    case 2
        %VAF uncentered
        SSE = sum(FR_error.^2); %Calculate sum of squared errors
        SSM = sum(FR_recorded.^2); %Calculate sum of squares 
        J = SSE/SSM; %Calculate cost
    case 3
        R = corrcoef(FR_sim,FR_recorded);  % Correlation coefficient for fit
        if numel(R) ==1
            J = R^2;                % Sometimes R only has 1 value...
        else
            J = R(2)^2;             % R-squared for fit
        end
    case 4
        J = max(FR_error.^2);
    case 5
        SSE = sum(FR_error.^2); %Calculate sum of squared errors
        SSM = sum(FR_recorded.^2); %Calculate sum of squares about the mean
        [maxR,imax] = max(FR_recorded);
        J = (SSE/SSM)+(maxR-FR_sim(imax))/50; %Calculate cost
end

