function [num_perts,ignore,proc_data] = septrials(filename)
% This function is designed to take data from cat 1A afferents
% experiments and separate the trials based on perturbation onsets and ends.
% Because the perturbations are of different temporal lengths, each trial is 
% saved as a "sub-structure." Each sub-structure contains all data for that 
% particular perturbation. After separating the trials, this function
% saves the data into the 'stacked_data' structure and appends it to
% 'filename'. 
%
%
% Created: 3/2013
% Author: Kyle Blum 
% Last Modified: 10/2013
% Modifications: 
%   1) Edited to accommodate triangle perturbation data. These
%      perturbations were performed in two groups of 5 triangles, as 
%      opposed to individual trials. septrials will now determine if the 
%      perturbation type is triangle, and if so, separate the data
%      accordingly. 

load(filename)
proc_data.time_PotMem = PotMem.times();
fs_PotMem = 20000;
proc_data.time = dFdt.times();
fs = 2000;

v2mm = 2; %mm per V
v2vel = 400; % mm/s per V WARNING: LENA'S CONVERSION (200) IS INCORRECT (KPB 10/2016)
v2force = 10; % N per V
v2dFdt = 2000; % N/s per V

if ~exist('RampOn')                     % Ignore bad trials
    num_perts = 0;
    ignore = 1;
    return
elseif strfind(filename,'pyramid')      % Ignore pyramid trials...RampOn.times don't match up
    num_perts = 0;
    ignore = 1;
    return
else
    ignore = 0;
end

if strfind(filename,'triangle') > 0;    % Triangles are special case
    time_start = RampOn.times([1,11]);
    time_end = RampOn.times([10,20]);
    num_perts = length(time_start);
else
    if RampOn.times(1) > 0.0001;             % If the first RampOn value isn't ~0
        time_start = RampOn.times(1:2:end);
        time_end = RampOn.times(2:2:end);
        num_perts = length(time_start)-1;
    elseif RampOn.times(1) < 0.0001;         % If the first RampOn value is ~0
        time_start = RampOn.times(2:2:end);
        time_end = RampOn.times(3:2:end);
        num_perts = length(time_start)-1;
    end
    if time_end(end)-time_start(end) < 0.01; % Get rid of meaningless RampOn values
       time_start = time_start(1:end-1);
       time_end = time_end(1:end-1);
    end
end

% buffer = 0.5;   % Buffer (in seconds) to add before and after perturbation
buffer = 0.05;   % Buffer (in seconds) to add before and after perturbation


ts = round(time_start*fs-buffer*fs); % Perturbation on (2000 Hz)
% if strfind(filename,'triangle')
    tf = round(time_end*fs+buffer*fs);   % Perturabtion off (2000 Hz)
% else
%     tf = round((time_end*fs - time_start*fs)/2 + time_start*fs); % End pert. at midline of plateau
% end


ts_PotMem = round(time_start*fs_PotMem-buffer*fs_PotMem); % Perturbation on (20000 Hz)
% if strfind(filename,'triangle')
    tf_PotMem = round(time_end*fs_PotMem+buffer*fs_PotMem);   % Perturbation off (20000 Hz)
% else
%     tf_PotMem = round((time_end*fs_PotMem - time_start*fs_PotMem)/2 + time_start*fs_PotMem); % End pert. at midline of plateau
% end
% raw_Accel = Accel.values;                                     %Unfiltered Accel from accelerometer
% [Ba,Aa] = butter(3,[10/(2000/2),30/(2000/2)]);                %Band-pass butterworth filter
% filt_Accel = filtfilt(Ba,Aa,raw_Accel);                       %Filtered Accel
    
t = proc_data.time;
V = Velocity.values.*v2vel;
Wnvel = 40/(0.5*fs);
[Bvel,Avel] = butter(5,Wnvel,'low');
filt_velocity = filtfilt(Bvel,Avel,V);


dVdt = diff(filt_velocity)./diff(t);                                      %Estimation of Accel. via difference approximation
dVdt(end+1) = dVdt(end);                                                  %This is to make sure dVdt is same length as other variables
Wnacc = 40/(0.5*fs);
[B_dV,A_dV] = butter(6,Wnacc,'low');
filt_dVdt = filtfilt(B_dV,A_dV,dVdt);

F = Force.values*v2force;
WnF = 50/(0.5*fs);
[BF,AF] = butter(3,WnF,'low');
filt_F = filtfilt(BF,AF,F);

L = Length.values*v2mm;

a = 6; % Exponential stiffness of tendon (hc = high compliance)
LengthT = 1/a*log(a*F); % Estimate of tendon length assuming exponentially elastic tendon
LengthT = LengthT - LengthT(1);
LengthF = L - L(1) - LengthT; % Estimate of change in fascicle length (L - L0) with exponentially elastic tendon

tL = Length.times;
raw_VF = diff(LengthF)./diff(tL);
raw_VF(end+1) = raw_VF(end);
filt_VF = filtfilt(Bvel,Avel,raw_VF);

raw_dVFdt = diff(filt_VF)./diff(tL);
raw_dVFdt(end+1) = raw_dVFdt(end);
filt_dVFdt = filtfilt(B_dV,A_dV,raw_dVFdt);

a_hc = 2; % Exponential stiffness of tendon (hc = high compliance)
LengthT_hc = 1/a_hc*log(a_hc*F); % Estimate of tendon length assuming exponentially elastic tendon
LengthT_hc = LengthT_hc - LengthT_hc(1);
LengthF_hc = L - L(1) - LengthT_hc; % Estimate of change in fascicle length (L - L0) with exponentially elastic tendon

raw_VF_hc = diff(LengthF_hc)./diff(tL);
raw_VF_hc(end+1) = raw_VF_hc(end);
filt_VF_hc = filtfilt(Bvel,Avel,raw_VF_hc);

raw_dVFdt_hc = diff(filt_VF_hc)./diff(tL);
raw_dVFdt_hc(end+1) = raw_dVFdt_hc(end);
filt_dVFdt_hc = filtfilt(B_dV,A_dV,raw_dVFdt_hc);

raw_dFdt = -dFdt.values.*v2dFdt;                                      %Unfiltered dF/dt
raw_dFdt = raw_dFdt - raw_dFdt(1);
WndF = 100/(0.5*fs);
[BdF,AdF] = butter(6,WndF,'low');
filt_dFdt = filtfilt(BdF,AdF,raw_dFdt);
% filt_dFdt = raw_dFdt;




for pert = 1:num_perts 
    proc_data(pert).time = dFdt.times(ts(pert):tf(pert));
    proc_data(pert).time_PotMem = PotMem.times(ts_PotMem(pert):tf_PotMem(pert));
    proc_data(pert).PotMem = PotMem.values(ts_PotMem(pert):tf_PotMem(pert));
%     proc_data(pert).Accel = filt_Accel(ts(pert):tf(pert)) - mean(filt_Accel(ts(pert):ts(pert)+100));
    proc_data(pert).ENG = ENG.values(ts(pert):tf(pert));
    proc_data(pert).Force = filt_F(ts(pert):tf(pert)).*v2force;
    proc_data(pert).Length = Length.values(ts(pert):tf(pert)).*v2mm;
%     proc_data(pert).Velocity = Velocity.values(ts(pert):tf(pert));
    proc_data(pert).Velocity = filt_velocity(ts(pert):tf(pert)) - mean(filt_velocity(ts(pert):ts(pert)+100));
    proc_data(pert).dVdt = filt_dVdt(ts(pert):tf(pert)) - mean(filt_dVdt(ts(pert):ts(pert)+100));
    proc_data(pert).dFdt = filt_dFdt(ts(pert):tf(pert)) - mean(filt_dFdt(ts(pert):ts(pert)+100));
    proc_data(pert).LengthT = LengthT(ts(pert):tf(pert));
    proc_data(pert).LengthF = LengthF(ts(pert):tf(pert));
    proc_data(pert).VelocityF = filt_VF(ts(pert):tf(pert));
    proc_data(pert).dVFdt = filt_dVFdt(ts(pert):tf(pert));
    proc_data(pert).LengthT_hc = LengthT_hc(ts(pert):tf(pert));
    proc_data(pert).LengthF_hc = LengthF_hc(ts(pert):tf(pert));
    proc_data(pert).VelocityF_hc = filt_VF_hc(ts(pert):tf(pert));
    proc_data(pert).dVFdt_hc = filt_dVFdt_hc(ts(pert):tf(pert));
    
    if strfind(filename,'2_accel') > 0;
        try
            proc_data(pert).baseline.bl_length = mean(proc_data(pert).Length(1:round(0.25*fs)));         % Take the mean of the first 0.25 seconds of data to get "baseline" values
            proc_data(pert).baseline.bl_force = mean(proc_data(pert).Force(1:round(0.25*fs)));
        catch
            proc_data(pert).baseline.bl_length = NaN;
            proc_data(pert).baseline.bl_force = NaN;
        end
    end
    

    
end


