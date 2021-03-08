function [proc_data] = findspikes(filename,num_perts)
% This function is a spike-detection algorithm designed 
% to detect action potentials collected from Ia afferents.
% findspikes will load file containing at least Vm and time data. 
%
% The algorithm uses a 2nd order band-pass butterworth filter
% with a cutoff wn = [0.05,0.2] to remove any low frequency 
% components. Minimum amplitude threshold and minimum spike 
% distances are used to perform a search for spikes. This appends ISI,
% firing rate (calculated from ISI), and spiketimes (times at which spikes 
% occur) into 'spikes' data structure and appends it to 'filename'.
%
% Modifications: 
% 
% Created: 3/2013
% Author: Kyle Blum

load(filename)
fs = 20000; %sampling frequency of voltage in Hz

for pert = 1:num_perts;
    time = proc_data(pert).time_PotMem; 
    Vm = proc_data(pert).PotMem;
    
    %%% High-Pass Butterworth Filter %%%
    order = 2;                          % order of filter
    [b,a] = butter(order,[0.05,0.2]); % filter transfer function
    Vm_filt = filtfilt(b,a,Vm);         % filter membrane voltage data 
                                        % (zero-phase)
    spike_thresh = 3*std(Vm_filt);    % Threshold of the noise
    spikewidth = (2e-3)*fs;             % Minimum spike width in # datapoints

    [~,index] = findpeaks(Vm_filt,'MINPEAKHEIGHT',spike_thresh,...
        'MINPEAKDISTANCE',spikewidth); % Search for peaks using specified
                                       % thresholds
   times = time(index);                                 
   
   spiketimes = times;                  %  spiketimes
   ISI = diff(spiketimes); % 'ISI'
   firing_rate = 1./ISI;   % 'firing_rate'
   
%    null = find(firing_rate < 50);
%    
%    spiketimes(null+1) = [];
%    ISI(null) = [];
%    firing_rate(null) = [];
   
   proc_data(pert).PotMem_filt = Vm_filt;
   proc_data(pert).spiketimes = spiketimes;
   proc_data(pert).ISI = ISI;
   proc_data(pert).firing_rate = firing_rate;

end
end