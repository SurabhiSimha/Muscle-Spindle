clear

load('C:\Users\Giovanni\Desktop\Spindle Group\Fitting model\PLoSCompBioData\aff02_mn0_04_accel_series_proc.mat')
perts = 9;

ISI = diff(data(perts).spiketimes); % 'ISI'
data2.firing_rate = [0; 1./ISI];   % 'firing_rate'

data2.Force = data(perts).Force;
data2.dFdt = [0; diff(data(perts).Force)];
data2.Length = data(perts).Length;
data2.Velocity = data(perts).Velocity;
data2.dVdt = [0; diff(data(perts).Velocity)];

data2.Force = data2.Force(round((data(perts).spiketimes-data(perts).spiketimes(1))*2000)+55);
data2.dFdt = data2.dFdt(round((data(perts).spiketimes-data(perts).spiketimes(1))*2000)+55);
data2.Length = data2.Length(round((data(perts).spiketimes-data(perts).spiketimes(1))*2000)+55);
data2.Velocity = data2.Velocity(round((data(perts).spiketimes-data(perts).spiketimes(1))*2000)+55);
data2.dVdt = data2.dVdt(round((data(perts).spiketimes-data(perts).spiketimes(1))*2000)+55);
    
for model = 1:2          % For all different models
    
    flags.model = model;
    flags.aff = 32;
    
    switch flags.model
        case 1 % Force, dFdt, static model
            
            gain_limits(1,:) = [1,1,0,0];                    %Initial guesses for [k_F,k_dF,Static,F_offset,dF_offset]
            gain_limits(2,:) = [0.0,0.0,-500,-500];                   %Lower bounds for [k_F,k_dF,Static,F_offset,dF_offset]
            gain_limits(3,:) = [1000,1000,0,0];         %Upper bounds [k_F,k_dF,Static,F_offset,dF_offset]
            
            
        case 2 % Kinematic Model (dV/dt,Vel,Len)
            gain_limits = [];
            gain_limits(1,:) = [1,1,1,1,1,1];                       %Initial guesses for [k_p,k_v,k_a,p_offset,v_offset,a_offset]
            gain_limits(2,:) = [0,0,0,-1000,-1000,-1000];           %Lower bounds for [k_p,k_v,k_a,p_offset,v_offset,a_offset]
            gain_limits(3,:) = [10000,10000,10000,1000,1000,1000];  %Upper bounds for [k_p,k_v,k_a,p_offset,v_offset,a_offset]
    end
    
    [opt_gains,terminal_cost,SSE,ybar,SSM,R_squared,exitflag,output] = spindle_findGains(gain_limits,data2,flags);
    
    switch flags.model
        case 1   % Force Model
            opt_fit = kinetics(data2,opt_gains,flags);
        case 2   % Kinematic Model
            opt_fit = kinematics(data2,opt_gains,flags);
    end
    
    figure
    subplot(4,1,1)
    plot(data(perts).spiketimes(2:end),data2.firing_rate(2:end),'o')
    ylim([0 400])
    subplot(4,1,1)
    plot(data(perts).spiketimes(2:end),opt_fit(2:end))
    ylim([0 400])
    switch flags.model
        case 1   % Force Model
            subplot(4,1,2)
            plot(data(perts).spiketimes(2:end),data2.Force(2:end))
            subplot(4,1,3)
            plot(data(perts).spiketimes(2:end),data2.dFdt(2:end))
        case 2   % Kinematic Model
            subplot(4,1,2)
            plot(data(perts).spiketimes(2:end),data2.Length(2:end))
            subplot(4,1,3)
            plot(data(perts).spiketimes(2:end),data2.Velocity(2:end))
            subplot(4,1,4)
            plot(data(perts).spiketimes(2:end),data2.dVdt(2:end))
    end
end