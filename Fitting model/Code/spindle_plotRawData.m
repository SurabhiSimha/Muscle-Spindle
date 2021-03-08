% clc, clear all

% pathname = cd(cd(['..' filesep 'Data' filesep]));
% filename = 'aff32_mn0_02_accel_ramp_proc_no_buffer.mat';          % Ramp and hold
% loadname = [pathname filesep filename];
%
% load(loadname);



% for pert = [93 92 49 48 11];
for pert = [7];
    figure, hold on
    t = proc_data(pert).time - proc_data(pert).time(1);
    ts = proc_data(pert).spiketimes(2:end) - proc_data(pert).time(1);
    ifr = proc_data(pert).firing_rate;
    
    subplot(6,2,1:2),hold on%, title(info.name,'interpreter','none')
    axis([t(1) t(end) 0 300])
    line(ts,ifr,'Marker','.','LineStyle','none','color',[0 0 0])
    plot([ts,ts],[0,0.1*max(ifr)],'k'), ylabel('Firing Rate (Hz)')
    
    subplot(6,2,3:4),hold on
    axis([t(1) t(end) 0 max(proc_data(pert).Force)+0.1])
    line(t,proc_data(pert).Force), ylabel('Force (N)')
    
    subplot(6,2,5:6)
    axis([t(1) t(end) min(proc_data(pert).dFdt)-1 max(proc_data(pert).dFdt)+1])
    line(t,proc_data(pert).dFdt), ylabel('dF/dt (N/s)')
    
    
    subplot(6,2,7:8)
    axis([t(1) t(end) 0 3.5])
    line(t,proc_data(pert).Length), ylabel('Length (mm)')
    
    subplot(6,2,9:10)
    axis([t(1) t(end) -12 12])
    line(t,proc_data(pert).Velocity), ylabel('Velocity (mm/s)')
    
    subplot(6,2,11:12)
    axis([t(1) t(end) -700 700])
    line(t,proc_data(pert).dVdt), ylabel ('Accel. (mm/s^2)'), xlabel('time (s)')
end
