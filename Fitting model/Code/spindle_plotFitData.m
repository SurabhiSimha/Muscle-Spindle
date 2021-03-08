function spindle_plotFitData(filename,ifrfilename,perts)

load(filename)
load(['..' filesep 'Data' filesep 'stats_no_buffer' filesep 'mean_lags.mat']);
ifrdata = load(ifrfilename);

aff = str2double(info.affname(4:5));
for pert = perts;
    
    clearvars -except proc_data info fit_data pert y lagind lagindkin lagindkinc mean_lags aff ifrdata
    
    lag = mean_lags(mean_lags(:,1)==aff,2);
    lagkin = mean_lags(mean_lags(:,1)==aff,3);
    lagprochazka = mean_lags(mean_lags(:,1)==aff,4);
    
    figure, hold on
    
    t = proc_data(pert).time - proc_data(pert).time(1);
    
    index = 1:round(t(end)*2000);
    
    tfit = t(index);
    ts = ifrdata.proc_data(pert).spiketimes(2:end) - proc_data(pert).time(1);
    ifr = ifrdata.proc_data(pert).firing_rate;
    
    subplot(6,2,1:2),hold on%, title(info.name,'interpreter','none')
    axis([-0.5 4 0 300])
    line(ts,ifr,'Marker','.','LineStyle','none','color',[0 0 0])
    plot([ts,ts],[0,0.1*max(ifr)],'k'), ylabel('Firing Rate (Hz)')
    
    F = proc_data(pert).Force(index);
    dF = proc_data(pert).dFdt(index);
    L = proc_data(pert).Length(index);
    V = proc_data(pert).Velocity(index);
    dV = proc_data(pert).dVdt(index);
    
        
        %Need to get rid of short, >0 bumbps in dF/dt
        %without filtering signal
        posdF = [0;dF>0;0]; %index where dF>0
        posneg = find(diff(posdF)==-1); %Find where dF changes from positive to negative value
        negpos = find(diff(posdF)==1); % Find where dF changes from negative to positive value
        ind = posneg-negpos < 80;  % Find sections of dF > 0 that are shorter than 100 samples (50ms)
        if strcmp(info.affname,'aff25')
            ind = posneg-negpos<60; % Had to shorten for aff25 due to short, fast perturbations
        end
        posneg = posneg(ind);       % Remove sections of dF > 0 that are longer than 100 samples
        negpos = negpos(ind);       % "
        if dF(1) > 0
            for k = 1:numel(posneg) % "
                dF(negpos(k):posneg(k)) = 0; % "
            end % "
        elseif dF(1) <= 0
            for k = 1:numel(posneg) % "
                dF(posneg(k):negpos(k)) = 0; % "
            end % "
        end
        
        posdV = [0;dV>0;0]; %index where dF>0
        dVposneg = find(diff(posdV)==-1); %Find where dF changes from positive to negative value
        dVnegpos = find(diff(posdV)==1); % Find where dF changes from negative to positive value
        dVind = dVposneg-dVnegpos < 40;  % Find sections of dF > 0 that are shorter than 80 samples (40ms)
        dVposneg = dVposneg(dVind);       % Remove sections of dF > 0 that are longer than 80 samples
        dVnegpos = dVnegpos(dVind);       % "
        if dV(1) > 0
            for k = 1:numel(dVposneg) % "
                dV(dVnegpos(k):dVposneg(k)) = 0; % "
            end % "
        elseif dV(1) < 0
            for k = 1:numel(dVposneg) % "
                dV(dVposneg(k):dVnegpos(k)) = 0; % "
            end % "
        end
        
        %             F(dF>0) = 0; %Assume dFdt dominates signal where it is nonzeros
        if strcmp(info.affname,'aff25') || strcmp(info.affname,'aff26') || strcmp(info.affname,'aff58') || strcmp(info.affname,'aff60')
            
            dF(dF < 0) = 0;
            F(dF ~=0) = 0;
        end
        
    if numel(dF) > numel(F)
        dF(end - (numel(dF)-numel(F))) = [];
    end
    if numel(dV) > numel(L)
        dV(end - (numel(dV)-numel(L))) = [];
    end
    
    F_component = fit_data(1).opt_gains(pert,1)*F + fit_data(1).opt_gains(pert,3);
    dFdt_component = fit_data(1).opt_gains(pert,2)*dF + fit_data(1).opt_gains(pert,4);
    
    dFdt_component(dFdt_component<=fit_data(1).opt_gains(pert,4)) = 0;
    F_component(F_component<= fit_data(1).opt_gains(pert,3)) = 0;
    
    if strcmp(info.affname,'aff25') || strcmp(info.affname,'aff26') || strcmp(info.affname,'aff58') || strcmp(info.affname,'aff60')
        F_component(dFdt_component > 0) = 0;
    end
    
%     if fit_data(1).opt_gains(pert,3) < 0;
         F_component(F_component <= 0) = 0;
%     elseif fit_data(1).opt_gains(pert,4) < 0;
         dFdt_component(dFdt_component <=0) = 0;
%     end
    
    FR_sim = F_component + dFdt_component; % Estimated firing rate from kinetics information
    
    L_component = fit_data(2).opt_gains(pert,1)*L + fit_data(2).opt_gains(pert,4);
    V_component = fit_data(2).opt_gains(pert,2)*V + fit_data(2).opt_gains(pert,5);
    A_component = fit_data(2).opt_gains(pert,3)*dV + fit_data(2).opt_gains(pert,6);
    
    L_component(L_component <= fit_data(2).opt_gains(pert,4)) = 0;
    V_component(V_component <= fit_data(2).opt_gains(pert,5)) = 0;
    A_component(A_component <= fit_data(2).opt_gains(pert,6)) = 0;
    
    FR_kin = L_component + V_component + A_component; % Estimated firing rate from L,V,A components
    
    FR_kinc = fit_data(3).opt_gains(pert,1).*abs(V).^fit_data(3).opt_gains(pert,2) + fit_data(3).opt_gains(pert,3)*L + fit_data(3).opt_gains(pert,4);
        
    %          line(tfit+lagtimes(lag),F_component,'color',[1 0 0])
    %          line(tfit+lagtimes(lag),dFdt_component,'color',[0 1 0])
    %         line(tfit+0.001*lag,static_component,'color',[0.5 0 0.5])
    line(tfit+lag,FR_sim)
    line(tfit+lagkin,FR_kin,'color',[1 0 0])
    line(tfit+lagprochazka,FR_kinc,'color',[1 0 1])
    text(t(round(length(tfit)/2))-0.5,max(ifr)+0.1*max(ifr),num2str(fit_data(1).R_squared(pert)))
    text(t(round(length(tfit)/2)),max(ifr)+0.1*max(ifr),num2str(fit_data(2).R_squared(pert)))
    text(t(round(length(tfit)/2))+0.5,max(ifr)+0.1*max(ifr),num2str(fit_data(3).R_squared(pert)))
    

    subplot(6,2,3:4),hold on
    axis([-0.5 4 0 450])
    line(t,proc_data(pert).Force/10), ylabel('Force (N)')
    
    subplot(6,2,5:6)
    axis([-0.5 4 -300 300])
    line(t,proc_data(pert).dFdt), ylabel('dF/dt (N/s)')
    
    subplot(6,2,7:8)
    axis([-0.5 4 0 3.5])
    line(t,proc_data(pert).Length), ylabel('Length (mm)')
    
    subplot(6,2,9:10)
    axis([-0.5 4 -20 20])
    line(t,proc_data(pert).Velocity), ylabel('Velocity (mm/s)')
    
    subplot(6,2,11:12)
    axis([-0.5 4 -1200 1200])
    line(t,proc_data(pert).dVdt), ylabel ('Accel. (mm/s^2)'), xlabel('time (s)')
    
    
end
