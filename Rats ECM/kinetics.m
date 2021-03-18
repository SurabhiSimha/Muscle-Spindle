 function [f,Ffib,dF] = kinetics(data,gains,flags)

F = data.Force;         %Recorded force value (2kHz)
L = data.L - data.L(1); %Relative length change applied to MT (2kHz)


Fpas = gains(5)*exp(gains(6)*L); %Estimate noncontractile tissue force

Fpas = Fpas - Fpas(1);  %Make relative to initial value
Ffib = F - Fpas;        %Estimate contractile tissue force as F_MT - F_NC
Ffib = Ffib + gains(3); %Apply F offset

if gains(3) < 0                        %if the offset is negative...
    Ffib(Ffib <= gains(3)) = gains(3); %half-wave rectify to value of offset
else                                   %otherwise...
    Ffib(Ffib < 0 ) = 0;               %half-wave rectify to 0
end

dF = diff(Ffib)*data.fs; %Calculate derivative of force, dF/dt
[B_df,A_df] = butter(3,50/(data.fs/2),'low'); %Set filter parameters           
dF = filtfilt(B_df,A_df,dF);                  %Filter dF/dt

Ffib = gains(1)*Ffib;    %Apply F gain

% Resample Force variables if we are curve-fitting
if flags.plot == 0          % Are we calling this function for fitting?
    Ffib = Ffib(ceil(data.IFRtimes*data.fs)); %Resample contractile F
    dF = dF(ceil(data.IFRtimes*data.fs));     %Resample contractile dF/dt
elseif flags.plot == 1      % Are we calling this function for plotting?
    dF(end+1) = dF(end);    % Make dF/dt same length as other variables
end

dF = gains(2)*dF + gains(4);  % Apply gains and offsets to dF/dt
dF(dF <= gains(4)) = gains(4);% Half-wave rectify to value of offset


f = Ffib + dF; % Add F and dF/dt components
f(f < 0) = 0;  % This is because we allow dF & F offsets to be < 0, so this makes the min values 0