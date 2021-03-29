%Created 08/01/2017
%Author: KPB
%Modified by GM 03/18/2021
%Purpose: The purpose of this script is as a bare-bones driver for
%curve-fitting force- and length-related models to recorded IFRs from
%rat muscle spindles. Models are specified in the 'kinetics.m' and
%'kinematics.m' files. This script requires the user to specify which model
%they wish to use for the curve fitting (default models: 1-kinetics; 
%2-kinematics). A 'proc_data' file containing processed experimental data
%should be loaded into the MATLAB workspace prior to running this script. 

clc, tic
disp('Running Code...')
clear
load aff101_sample.mat

%%% STEPS FOR CURVE FITTING SPINDLE DATA %%%
% Load afferent .mat file into MATLAB workspace
% Figure out which models you want to use for curve-fitting (specify below)
% Figure out which trials you want to (specify below)

model = [2]; %Which model do you want to use for fitting?
pert = [1];  %Which perturbation do you want to fit?

flags.model = model; %Set a flag to signal which model to use later
flags.cost = 2; %cost function: 1 SSE/SST centered, 2 SSE/SST uncentered, 3 corrcoeff^2 4 min (max error)
flags.alg = 1; %1 fmincon 2 BADS 

%% step 1 
%%% Initialize optimization parameters %%%
gain_limits(1,:) = [0,0,0,0,0.1,0.1];       %Initial guesses for [k_F,k_dF,o_F,o_dF,A_pas,K_pas]
gain_limits(2,:) = [0,0,-1000,-1000,0.0,0.0]; %Lower bounds for [k_F,k_dF,o_F,o_dF,A_pas,K_pas]
gain_limits(3,:) = [0,0,1000,1000,1,2]; %Upper bounds [k_F,k_dF,o_F,o_dF,A_pas,K_pas]

%%% Resample kinematics data at IFR times %%%
data.IFR = proc_data(pert).firing_rate;                            % Firing Rate
data.time_end = proc_data(pert).time(end)-proc_data(pert).time(1); % Trial Duration
fs = 1/(proc_data(pert).time(2) - proc_data(pert).time(1)); %Sampling frequency

%%% Continuous Force & Length Data (resampled after F_contractile is estimated later) %%%
data.L = proc_data(pert).Length;      % Continuous length (for passive tissue force calculation)
data.Force = proc_data(pert).Force;   % Recorded MT Force

data.time = proc_data(pert).time - proc_data(pert).time(1); %trial time
data.fs = fs;                                               %trial sampling freq
flags.plot = 0;                                             %don't plot data

%%% Spike timing (IFR OCCURS AT SECOND SPIKE OF EACH PAIR) %%%
data.IFRtimes = proc_data(pert).spiketimes(2:end) - proc_data(pert).time(1);
[~, maxY] = max([0; diff(data.Force)]);
[~, maxIFR] = max(data.IFR);
lagtime = data.IFRtimes(maxIFR)-data.time(maxY);
data.IFRtimes = proc_data(pert).spiketimes(2:end) - proc_data(pert).time(1) - lagtime;

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

[tempf,tempFfib,tempdF,kF,kdF] = kinetics2(data,opt_gains(pert,:),flags);
opt_gains(pert,2) = kdF;
opt_gains(pert,1) = kF;

[tempf,tempFfib,tempdF] = kinetics(data,opt_gains(pert,:),flags);

figure
subplot(5,1,1)
plot(data.time,data.Force)
ylabel("Recorded Force")
xlabel("time")

subplot(5,1,2)
plot(IFR(pert).times,IFR(pert).values,'o')
ylabel("IFR")
xlabel("time")

subplot(5,1,2)
plot(data.IFRtimes,tempf)
legend(["recorded","simulated"])

subplot(5,1,4)
plot(data.IFRtimes,tempFfib)
ylabel("Fiber Force")
xlabel("time")

subplot(5,1,5)
plot(data.IFRtimes,tempdF)
ylabel("Fiber Yank")
xlabel("time")

