clear

load('C:\Users\Giovanni\Desktop\Spindle Group\Fitting model\PLoSCompBioData\aff02_mn0_04_accel_series_proc.mat')
perts = 9;

ISI = diff(data(perts).spiketimes); % 'ISI'
data3.firing_rate = [0; 1./ISI];   % 'firing_rate'

data2.Force = data(perts).Force;
data2.dFdt = [0; diff(data(perts).Force)];
data2.Length = data(perts).Length;
data2.Velocity = data(perts).Velocity;
data2.dVdt = [0; diff(data(perts).Velocity)];
data2.time = data(perts).time;

data3.Force = data2.Force(round((data(perts).spiketimes-data(perts).spiketimes(1))*2000)+55);
data3.dFdt = data2.dFdt(round((data(perts).spiketimes-data(perts).spiketimes(1))*2000)+55);
data3.Length = data2.Length(round((data(perts).spiketimes-data(perts).spiketimes(1))*2000)+55);
data3.Velocity = data2.Velocity(round((data(perts).spiketimes-data(perts).spiketimes(1))*2000)+55);
data3.dVdt = data2.dVdt(round((data(perts).spiketimes-data(perts).spiketimes(1))*2000)+55);

flags.model = 1;
flags.aff = 32;

gain_limits(1,:) = [1,1,0,0];                    %Initial guesses for [k_F,k_dF,Static,F_offset,dF_offset]
gain_limits(2,:) = [0.0,0.0,-500,-500];                   %Lower bounds for [k_F,k_dF,Static,F_offset,dF_offset]
gain_limits(3,:) = [1000,1000,0,0];         %Upper bounds [k_F,k_dF,Static,F_offset,dF_offset]

[opt_gains,terminal_cost,SSE,ybar,SSM,R_squared,exitflag,output] = spindle_findGains(gain_limits,data3,flags);

a=1;
for i = [0.5 0.75 1 2 3]
    for ii = [0.5 0.75 1 2 3]
        opt_fit = kinetics(data2,opt_gains.*[i ii 1 1],flags);

        subplot(5,5,a)
        [st, IFR] = spiketrain(data2.time,opt_fit/60);
        plot(st(2:end),IFR,'o')

        subplot(5,5,a)
        plot(data2.time,opt_fit/4)        
        ylim([0 400])
        a = a+1;
    end
end

