function [opt_gains,terminal_cost,SSE,ybar,SSM,R_squared,exitflag,output] = findGains(gain_limits,data,flags)
% Created: 5/2013
% Modified: 11/2015, 5/2016
% Author: Kyle Blum
% Description: This function performs the fitting of kinematics/force data
% to spindle IFRs. 
% User-defined dependencies: 'cost.m'

%%% Set up fmincon options %%%
optimize_options = optimset('display','off','algorithm',...
    'interior-point','TolX',1e-8,'TolFun',1e-7,'MaxFunEvals',10000);  

%%% Set up parameters %%%
gains_init = gain_limits(1,:);  % Initial guess for gains
lower_bound = gain_limits(2,:); % Lower bound of gains and time delay
upper_bound = gain_limits(3,:); % Upper bound of gains and time delay

%%% Set up optimization constraints %%%
A = [];
B = [];
Aeq = [];
Beq = [];
nonlconstr = [];

%%% RUN OPTIMIZATION WITH PARAMETERS AND CONSTRAINTS FROM ABOVE %%%
switch flags.alg
    case 1
        [opt_gains,test_cost,exitflag,output] = fmincon(@cost, gains_init, ...
            A, B, Aeq, Beq, lower_bound, upper_bound, nonlconstr, optimize_options, ...
            data,flags);
        % IMPORTANT NOTE:
        % Line 26 will only return final values for opt_gains, test_cost, etc.
        % Internally, fmincon.m will search for parameters until conditions for
        % stopping the optimization are met (e.g. tolerances, max function evals,
        % etc.). Stats are calculated based on final values only.
    case 2
        [opt_gains,test_cost,exitflag,output] = bads(@(gains_init)cost(gains_init,data,flags),gains_init,...
            lower_bound, upper_bound,[],[],[],[]);
end

%%% Find optimal gains %%%
switch flags.model
    case 1   % Force-related model optimal fit
        [opt_fit,~,~] = kinetics(data,opt_gains,flags); 
    case 2   % Length-related model optimal fit
        [opt_fit,~,~,~,~] = kinetics2(data,opt_gains,flags); 
end


FR_recorded = data.IFR;             % Recorded firing rate
R = corrcoef(opt_fit,FR_recorded);  % Correlation coefficient for fit
if numel(R) ==1
    R_squared = R^2;                % Sometimes R only has 1 value...
else
    R_squared = R(2)^2;             % R-squared for fit
end

SSE = sum((FR_recorded-opt_fit).^2);% Sum of squared errors for final fit
ybar = mean(FR_recorded);           % Mean value of IFR
SSM = sum((FR_recorded-ybar).^2);   % Sum of squares about the mean of IFR

terminal_cost = test_cost;          % Cost function value for optimal fit
