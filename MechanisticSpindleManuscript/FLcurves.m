clear

delta_f_activated = [zeros(1,200) 1 zeros(1,2800)];
delta_cdl = zeros(1,3001);
t = 0:0.001:3;

figure
for i = 1:40
    i
    % Make a half-sarcomere
    hsB = halfSarcBag();
    hsC = halfSarcChain();
    
    hsB.cmd_length = 1200+i*20;
    hsC.cmd_length = 1200+i*20;
    hsB.hs_length = 1200+i*20;
    hsC.hs_length = 1200+i*20;
    
    [hsB,dataB,hsC,dataC] = sarcSimDriverMod(t,delta_f_activated,delta_cdl,hsB,hsC);
    
    subplot(2,1,1)
    plot(dataB.hs_length(end),dataB.cb_force(1060),'o')
    subplot(2,1,2)
    plot(t,dataB.cb_force)
end
