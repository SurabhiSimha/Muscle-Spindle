function [DI,tdiff] = dynamicIndex(IFR,spiketimes,time,refTime)
% This function calculates dynamic index for an afferent response to a ramp
% and hold perturbation. Classic dynamic index (Matthews, 1963) is the
% difference in firing rate between the peak dynamic response (usually at
% the end of ramp stretch) and a point 0.5s later during the isometric
% "hold" portion of the stretch sequence. This function allows the user to
% adjust this 0.5s parameter to a any point along the plateau. 
%
% 'IFR' is the instantaneous firing rate for a given stretch trial in units
% of spikes/s
%
% 'spiketimes' should contain times at which an action potential was recorded
% for the stretch trial in units of s
%
% 'time' should contain the "continuous" time points for the stretch trial
% in units of s
%
% 'refTime' is the time, relative to the middle of the "hold" phase, to use 
% as the reference for calculating DI. 0 would give the middle of the
% plateau, 0.5 seconds after the end of the ramp. 0.5 would give the end of
% the plateau, 1 second after the end of the ramp. 


if numel(spiketimes) > numel(IFR)
    spiketimes = spiketimes(2:end);
end

spiketimes = spiketimes - time(1);

tRef = (time(end) - time(1))/2 + refTime; % Find reference time relative to trial start


[dynamicResponse,dynamicLoc] = max(IFR(spiketimes > 0.1)); %0.1 chosen to ignore initial burst in IFR
[~,plateauLoc] = min(abs(spiketimes - tRef));
plateauResponse = IFR(plateauLoc);

DI = dynamicResponse - plateauResponse;
tdiff = spiketimes(plateauLoc) - spiketimes(dynamicLoc);

end