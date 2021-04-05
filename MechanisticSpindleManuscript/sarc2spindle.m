function [r,rs,rd] = sarc2spindle(dataB,dataC,kFc,kFb,kYb,occlusion,threshold)
t = dataC.t; %time

% Force-dominant fiber (static fiber)
Fs = dataC.hs_force;
Fs(Fs<0) = 0;



% Yank-dominant fiber (dynamic fiber)
Fd = dataB.hs_force; %force
Fd(Fd<0) = 0; %threshold
Y = diff(Fd)./diff(t); %yank
Y(Y<0) = 0; %threshold
Y(end+1) = Y(end); %make Y same length as F
Y(Fd<8e4) = 0; %Chirp
% Y(Fd<8e4) = 0; %Most of the simulations


rs = Fs*kFc; %static component
rd = Fd*kFb + Y*kYb; %dynamic component

rs(rs<0) = 0;
rs = rs/(10^6);
rd = rd/(10^6);

if occlusion % Hypothesis that branches of Ia ending compete for total firing rate
    rsComp = rs;
    rsComp(rd>=rs) = 0.3*rsComp(rd>=rs);  % 0.3 chosen to match Banks et al. 1997
    
    rdComp = rd;
    rdComp(rs>rd) = 0.3*rdComp(rs>rd);
    
    r = rsComp + rdComp;
%     rs = rsComp;
%     rd = rdComp;
else
    r = rs + rd; %linear sum
end

%Firing Threshold
r = r - threshold;
r(r<0.0) = 0; 

end