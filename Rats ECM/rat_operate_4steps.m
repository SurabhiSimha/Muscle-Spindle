%Created 08/01/2017
%Author: KPB
%
%Purpose: The purpose of this script is as a bare-bones driver for
%curve-fitting force- and length-related models to recorded IFRs from
%rat muscle spindles. Models are specified in the 'kinetics.m' and
%'kinematics.m' files. This script requires the user to specify which model
%they wish to use for the curve fitting (default models: 1-kinetics; 
%2-kinematics). A 'proc_data' file containing processed experimental data
%should be loaded into the MATLAB workspace prior to running this script. 

clc, tic
disp('Running Code...')

load aff101_sample.mat

%%% STEPS FOR CURVE FITTING SPINDLE DATA %%%
% Load afferent .mat file into MATLAB workspace
% Figure out which models you want to use for curve-fitting (specify below)
% Figure out which trials you want to (specify below)

model = [1]; %Which model do you want to use for fitting?
pert = [1];  %Which perturbation do you want to fit?

flags.model = model; %Set a flag to signal which model to use later

%% step 1 - search non-contractile force parameters
%%% Initialize optimization parameters %%%

gain_limits(1,:) = [1,1,0,0,0.1,0.1];       %Initial guesses for [k_F,k_dF,o_F,o_dF,A_pas,K_pas]
gain_limits(2,:) = [1,1,0,0,0.0,0.0]; %Lower bounds for [k_F,k_dF,o_F,o_dF,A_pas,K_pas]
gain_limits(3,:) = [1,1,0,0,0.5,1]; %Upper bounds [k_F,k_dF,o_F,o_dF,A_pas,K_pas]

%%% Spike timing (IFR OCCURS AT SECOND SPIKE OF EACH PAIR) %%%
data.IFRtimes = proc_data(pert).spiketimes(2:end) - proc_data(pert).time(1);
fs = 1/(proc_data(1).time(2) - proc_data(1).time(1)); %Sampling frequency

%%% Resample kinematics data at IFR times %%%
data.IFR = proc_data(pert).firing_rate;                            % Firing Rate
data.time_end = proc_data(pert).time(end)-proc_data(pert).time(1); % Trial Duration
data.Length = proc_data(pert).Length(ceil(data.IFRtimes*fs));      % Length
data.Velocity = proc_data(pert).Velocity(ceil(data.IFRtimes*fs));  % Velocity
data.dVdt = proc_data(pert).dVdt(ceil(data.IFRtimes*fs));          % Filtered Accel

%%% Continuous Force & Length Data (resampled after F_contractile is estimated later) %%%
data.L = proc_data(pert).Length;      % Continuous length (for passive tissue force calculation)
data.Force = proc_data(pert).Force;   % Recorded MT Force
data.dFdt = proc_data(pert).dFdt;     % Recorded MT dF/dt

data.time = proc_data(pert).time - proc_data(pert).time(1); %trial time
data.fs = fs;                                               %trial sampling freq
flags.plot = 0;                                             %don't plot data


%%% RUN OPTIMIZATION CODE: RETURNS OPTIMAL GAINS, AND FIT STATS %%%
[opt_gains(pert,:),terminal_cost(pert),SSE(pert),ybar(pert),SSM(pert),R_squared(pert),exitflag(pert),output(pert)] = findGains(gain_limits,data,flags);

%%% STORE DATA THAT IS OVERWRITTEN EVERY PERTURBATION LOOP %%%
IFR(pert).times = data.IFRtimes;                      % Times at which we have a firing rate value
IFR(pert).values = data.IFR;                          % Recorded firing rate
fit_data(flags.model).fitting_data(pert) = data;      % Data used for fitting
fit_data(flags.model).IFR = IFR';                     % Firing Rate data with values and times
fit_data(flags.model).opt_gains = opt_gains;          % Optimal gains from spindle_findGains
fit_data(flags.model).terminal_cost = terminal_cost'; % Final cost value from spindle_findGains
fit_data(flags.model).exitflag = exitflag';           % Exit flag from fmincon
fit_data(flags.model).opt_output = output';           % information about iterations from fmincon
fit_data(flags.model).date = date;                    % date code was run
fit_data(flags.model).R_squared = R_squared';         % R_squared for fit

%% step 2 - fit force gain
%%% Initialize optimization parameters %%%

gain_limits(1,:) = [1,1,0,0,opt_gains(5),opt_gains(6)];       %Initial guesses for [k_F,k_dF,o_F,o_dF,A_pas,K_pas]
gain_limits(2,:) = [0,1,-500,0,opt_gains(5),opt_gains(6)]; %Lower bounds for [k_F,k_dF,o_F,o_dF,A_pas,K_pas]
gain_limits(3,:) = [300,1,500,0,opt_gains(5),opt_gains(6)]; %Upper bounds [k_F,k_dF,o_F,o_dF,A_pas,K_pas]

%%% Spike timing (IFR OCCURS AT SECOND SPIKE OF EACH PAIR) %%%
data.IFRtimes = proc_data(pert).spiketimes(2:end) - proc_data(pert).time(1);
fs = 1/(proc_data(1).time(2) - proc_data(1).time(1)); %Sampling frequency

%%% Resample kinematics data at IFR times %%%
data.IFR = proc_data(pert).firing_rate;                            % Firing Rate
data.time_end = proc_data(pert).time(end)-proc_data(pert).time(1); % Trial Duration
data.Length = proc_data(pert).Length(ceil(data.IFRtimes*fs));      % Length
data.Velocity = proc_data(pert).Velocity(ceil(data.IFRtimes*fs));  % Velocity
data.dVdt = proc_data(pert).dVdt(ceil(data.IFRtimes*fs));          % Filtered Accel

%%% Continuous Force & Length Data (resampled after F_contractile is estimated later) %%%
data.L = proc_data(pert).Length;      % Continuous length (for passive tissue force calculation)
data.Force = proc_data(pert).Force;   % Recorded MT Force
data.dFdt = proc_data(pert).dFdt;     % Recorded MT dF/dt

data.time = proc_data(pert).time - proc_data(pert).time(1); %trial time
data.fs = fs;                                               %trial sampling freq
flags.plot = 0;                                             %don't plot data


%%% RUN OPTIMIZATION CODE: RETURNS OPTIMAL GAINS, AND FIT STATS %%%
[opt_gains(pert,:),terminal_cost(pert),SSE(pert),ybar(pert),SSM(pert),R_squared(pert),exitflag(pert),output(pert)] = findGains(gain_limits,data,flags);

%%% STORE DATA THAT IS OVERWRITTEN EVERY PERTURBATION LOOP %%%
IFR(pert).times = data.IFRtimes;                      % Times at which we have a firing rate value
IFR(pert).values = data.IFR;                          % Recorded firing rate
fit_data(flags.model).fitting_data(pert) = data;      % Data used for fitting
fit_data(flags.model).IFR = IFR';                     % Firing Rate data with values and times
fit_data(flags.model).opt_gains = opt_gains;          % Optimal gains from spindle_findGains
fit_data(flags.model).terminal_cost = terminal_cost'; % Final cost value from spindle_findGains
fit_data(flags.model).exitflag = exitflag';           % Exit flag from fmincon
fit_data(flags.model).opt_output = output';           % information about iterations from fmincon
fit_data(flags.model).date = date;                    % date code was run
fit_data(flags.model).R_squared = R_squared';         % R_squared for fit

%% step 3 - lag and dFt
%%% Initialize optimization parameters %%%
oldcost = 1;
tempoptgain = opt_gains;
for lag = 1:36
    lagtime = 0:0.001:0.035;
    
    gain_limits(1,:) = [tempoptgain(1),10000,tempoptgain(3),0,tempoptgain(5),tempoptgain(6)];       %Initial guesses for [k_F,k_dF,o_F,o_dF,A_pas,K_pas]
    gain_limits(2,:) = [tempoptgain(1),10000,tempoptgain(3),0,tempoptgain(5),tempoptgain(6)]; %Lower bounds for [k_F,k_dF,o_F,o_dF,A_pas,K_pas]
    gain_limits(3,:) = [tempoptgain(1),1000000,tempoptgain(3),0,tempoptgain(5),tempoptgain(6)]; %Upper bounds [k_F,k_dF,o_F,o_dF,A_pas,K_pas]
    
    %%% Spike timing (IFR OCCURS AT SECOND SPIKE OF EACH PAIR) %%%
    data.IFRtimes = proc_data(pert).spiketimes(2:end) - proc_data(pert).time(1) - lagtime(lag);
    fs = 1/(proc_data(1).time(2) - proc_data(1).time(1)); %Sampling frequency
    
    %%% Resample kinematics data at IFR times %%%
    data.IFR = proc_data(pert).firing_rate;                            % Firing Rate
    data.time_end = proc_data(pert).time(end)-proc_data(pert).time(1); % Trial Duration
    data.Length = proc_data(pert).Length(ceil(data.IFRtimes*fs));      % Length
    data.Velocity = proc_data(pert).Velocity(ceil(data.IFRtimes*fs));  % Velocity
    data.dVdt = proc_data(pert).dVdt(ceil(data.IFRtimes*fs));          % Filtered Accel
    
    %%% Continuous Force & Length Data (resampled after F_contractile is estimated later) %%%
    data.L = proc_data(pert).Length;      % Continuous length (for passive tissue force calculation)
    data.Force = proc_data(pert).Force;   % Recorded MT Force
    data.dFdt = proc_data(pert).dFdt;     % Recorded MT dF/dt
    
    data.time = proc_data(pert).time - proc_data(pert).time(1); %trial time
    data.fs = fs;                                               %trial sampling freq
    flags.plot = 0;                                             %don't plot data
    
    
    %%% RUN OPTIMIZATION CODE: RETURNS OPTIMAL GAINS, AND FIT STATS %%%
    [opt_gains(pert,:),terminal_cost(pert),SSE(pert),ybar(pert),SSM(pert),R_squared(pert),exitflag(pert),output(pert)] = findGains(gain_limits,data,flags);
    
    if R_squared<oldcost        
        %%% STORE DATA THAT IS OVERWRITTEN EVERY PERTURBATION LOOP %%%
        IFR(pert).times = data.IFRtimes;                      % Times at which we have a firing rate value
        IFR(pert).values = data.IFR;                          % Recorded firing rate
        fit_data(flags.model).fitting_data(pert) = data;      % Data used for fitting
        fit_data(flags.model).IFR = IFR';                     % Firing Rate data with values and times
        fit_data(flags.model).opt_gains = opt_gains;          % Optimal gains from spindle_findGains
        fit_data(flags.model).terminal_cost = terminal_cost'; % Final cost value from spindle_findGains
        fit_data(flags.model).exitflag = exitflag';           % Exit flag from fmincon
        fit_data(flags.model).opt_output = output';           % information about iterations from fmincon
        fit_data(flags.model).date = date;                    % date code was run
        fit_data(flags.model).R_squared = R_squared';         % R_squared for fit
        fit_data(flags.model).lag = lagtime(lag);
    end
end

%% step 4 - fit +-10%
%%% Initialize optimization parameters %%%

gain_limits(1,:) = [fit_data(flags.model).opt_gains(1),fit_data(flags.model).opt_gains(2),...
    fit_data(flags.model).opt_gains(3),fit_data(flags.model).opt_gains(4),...
    fit_data(flags.model).opt_gains(5),fit_data(flags.model).opt_gains(6)];       %Initial guesses for [k_F,k_dF,o_F,o_dF,A_pas,K_pas]
gain_limits(2,:) = [fit_data(flags.model).opt_gains(1)-fit_data(flags.model).opt_gains(1)*0.1,...
    fit_data(flags.model).opt_gains(2)-fit_data(flags.model).opt_gains(2)*0.1,...
    fit_data(flags.model).opt_gains(3)-fit_data(flags.model).opt_gains(3)*0.1,...
    fit_data(flags.model).opt_gains(4)-fit_data(flags.model).opt_gains(4)*0.1,...
    fit_data(flags.model).opt_gains(5)-fit_data(flags.model).opt_gains(5)*0.1,...
    fit_data(flags.model).opt_gains(6)-fit_data(flags.model).opt_gains(6)*0.1]; %Lower bounds for [k_F,k_dF,o_F,o_dF,A_pas,K_pas]
gain_limits(3,:) = [fit_data(flags.model).opt_gains(1)+fit_data(flags.model).opt_gains(1)*0.1,...
    fit_data(flags.model).opt_gains(2)+fit_data(flags.model).opt_gains(2)*0.1,...
    fit_data(flags.model).opt_gains(3)+fit_data(flags.model).opt_gains(3)*0.1,...
    fit_data(flags.model).opt_gains(4)+fit_data(flags.model).opt_gains(4)*0.1,...
    fit_data(flags.model).opt_gains(5)+fit_data(flags.model).opt_gains(5)*0.1,...
    fit_data(flags.model).opt_gains(6)+fit_data(flags.model).opt_gains(6)*0.1]; %Upper bounds for [k_F,k_dF,o_F,o_dF,A_pas,K_pas]

%%% Spike timing (IFR OCCURS AT SECOND SPIKE OF EACH PAIR) %%%
data.IFRtimes = proc_data(pert).spiketimes(2:end) - proc_data(pert).time(1) - fit_data(flags.model).lag;
fs = 1/(proc_data(1).time(2) - proc_data(1).time(1)); %Sampling frequency

%%% Resample kinematics data at IFR times %%%
data.IFR = proc_data(pert).firing_rate;                            % Firing Rate
data.time_end = proc_data(pert).time(end)-proc_data(pert).time(1); % Trial Duration
data.Length = proc_data(pert).Length(ceil(data.IFRtimes*fs));      % Length
data.Velocity = proc_data(pert).Velocity(ceil(data.IFRtimes*fs));  % Velocity
data.dVdt = proc_data(pert).dVdt(ceil(data.IFRtimes*fs));          % Filtered Accel

%%% Continuous Force & Length Data (resampled after F_contractile is estimated later) %%%
data.L = proc_data(pert).Length;      % Continuous length (for passive tissue force calculation)
data.Force = proc_data(pert).Force;   % Recorded MT Force
data.dFdt = proc_data(pert).dFdt;     % Recorded MT dF/dt

data.time = proc_data(pert).time - proc_data(pert).time(1); %trial time
data.fs = fs;                                               %trial sampling freq
flags.plot = 0;                                             %don't plot data


%%% RUN OPTIMIZATION CODE: RETURNS OPTIMAL GAINS, AND FIT STATS %%%
[opt_gains(pert,:),terminal_cost(pert),SSE(pert),ybar(pert),SSM(pert),R_squared(pert),exitflag(pert),output(pert)] = findGains(gain_limits,data,flags);

%%% STORE DATA THAT IS OVERWRITTEN EVERY PERTURBATION LOOP %%%
IFR(pert).times = data.IFRtimes;                      % Times at which we have a firing rate value
IFR(pert).values = data.IFR;                          % Recorded firing rate
fit_data(flags.model).fitting_data(pert) = data;      % Data used for fitting
fit_data(flags.model).IFR = IFR';                     % Firing Rate data with values and times
fit_data(flags.model).opt_gains = opt_gains;          % Optimal gains from spindle_findGains
fit_data(flags.model).terminal_cost = terminal_cost'; % Final cost value from spindle_findGains
fit_data(flags.model).exitflag = exitflag';           % Exit flag from fmincon
fit_data(flags.model).opt_output = output';           % information about iterations from fmincon
fit_data(flags.model).date = date;                    % date code was run
fit_data(flags.model).R_squared = R_squared';         % R_squared for fit

[tempf,tempFfib,tempdF] = kinetics(data,fit_data(flags.model).opt_gains,flags);
subplot(5,1,1)
plot(data.time,data.Force)
plot(data.IFRtimes,tempFfib-(tempFfib(1)))

subplot(5,1,2)
plot(IFR(pert).times,IFR(pert).values,'o')

subplot(5,1,3)
plot(data.IFRtimes,tempf)

subplot(5,1,4)
plot(data.IFRtimes,tempFfib)

subplot(5,1,5)
plot(data.IFRtimes,tempdF)
