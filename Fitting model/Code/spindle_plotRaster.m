function spindle_plotRaster(binsize,trialID,pathname,affname)

numtrials = length(trialID);
fig;
hold on
xlabel('time')
ylabel('Trial ID')


for trial = 1:numtrials
    filename = [pathname filesep affname '_mn0_' trialID{trial} '_accel_ramp_proc.mat'];
    load(filename)
    for pert = 1:length(proc_data);
        IFR = proc_data(pert).firing_rate;
        spiketimes = proc_data(pert).spiketimes(2:end)-proc_data(pert).time(1);
        
        bin = 1;
        for binstart = 0.05:binsize:2.2
            binIFR((pert+(trial-1)*length(proc_data)),bin) = (length(spiketimes(spiketimes>binstart)) - length(spiketimes(spiketimes>binstart+binsize)))/binsize;
            bin = bin+1;
        end
        plot([spiketimes,spiketimes],[(pert+(trial-1)*length(proc_data))-0.5,(pert+(trial-1)*length(proc_data))+0.5],'k')
    end
    
end

stdevIFR = std(binIFR);

fig;
hold on
time = linspace(0.05,2.2,bin-1);
line(time,mean(binIFR))
line(time,mean(binIFR)+stdevIFR,'linestyle','--')
line(time,mean(binIFR)-stdevIFR,'linestyle','--')



end
