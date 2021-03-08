function [st, IFR] = spiketrain(t,r)
v_reset = -0.08; %reset voltage
E_l = -0.07; %Nernst leak potential in V
g_l = 2.5e-8; %leak conductance in S
tau = 2e-2; %time constant (s)
v(1) = -0.07; %initial condition
u = (r-0.03)*5e-1; %scale r from spindle model into current
% u = interp1(0:0.001:t(end),u,t); %upsample to higher precision
v_th = -0.0;   %threshold voltage for spike
tsim = t;
dt = 5e-5; %
for i = 1:1:length(t)-1
    v(i+1) = v(i) + dt*((u(i) - g_l*(v(i) - E_l))/tau);
    if v(end) >= v_th
        spike(i+1) = 1;
        v(end) = v_reset;
    else
        spike(i+1) = 0;
    end
end
sIdx = find(spike==1);
st = tsim(sIdx);
ISI = diff(st);
IFR = 1./ISI;
end