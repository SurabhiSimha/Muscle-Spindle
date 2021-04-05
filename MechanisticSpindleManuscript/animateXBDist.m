function animateXBDist(index,hs,data,mode)
time = data.t;
bins = hs.x_bins;
f = hs.f;
g = hs.g;
dist = data.bin_pops;
no_detached = data.no_detached;
cbf = data.cb_force;
hsf = data.hs_force;
pf = data.passive_force;
hsl = data.hs_length;
cmd = data.cmd_length;
f_act = data.f_activated;

hfig = fig;
set(hfig,'Units','Normalized');
set(hfig,'PaperPosition',[0.5 0.25 0.5 0.75],'Position',[0.5 0.25 0.4 0.6])


h1 = subplot(3,8,1:8); hold on;
set(gca,'TickDir','out','FontName','Helvetica','FontSize',10,'xtick',...
    [-50 -25 -20 -19 -18 -17 -16 -15 0 25 50])
title('A'), xlabel('cross bridge length (nm)','FontSize',15), ylabel('Fraction','FontSize',15)
axis([bins(1) bins(end) 0 max(max(dist))])

h2 = subplot(3,2,3:4); hold on;
set(h2,'TickDir','out','FontName','Helvetica','FontSize',10,...
    'xtick',[])
ylabel('Stress (MPa)','FontSize',15)
axis([time(1) time(index(end)) min(hsf)/10^6 max(hsf)/10^6])

h3 = subplot(3,2,5:6); hold on;
set(gca,'TickDir','out','FontName','Helvetica','FontSize',10)
ylabel('hs length (L_0)','FontSize',15), xlabel('time (s)','FontSize',15)
axis([time(1) time(index(end)) min(hsl)/1300 max(hsl)/1300])

align_Ylabels;

% h_1B = subplot(3,8,6); hold on;
% set(gca,'TickDir','out','FontName','Helvetica','FontSize',10,...
%     'xtick',[])
% title('D'),ylabel('Fraction')
% axis([0 2 0 max(max(dist))])
% 
% h_1C = subplot(3,8,7); hold on;
% set(gca,'TickDir','out','FontName','Helvetica','FontSize',10,...
%     'xtick',[])
% title('Binding Sites')
% axis([0 2 0 max(max(f_act))])



h_dist = []; % Attached cross-bridge distribution handle
h_time = []; % Time display handle
h_cbf = [];  % Cross-bridge force handle
h_hsf = [];  % Half-sarcomere force handle
h_pf = [];   % Passive force handle
h_hsl = [];  % Cross-bridge length handle
h_vline2 = []; % Cosmetic moving vertical line (subplot 2)
h_vline3 = []; % Cosmetic moving vertical line (subplot 3)
h_det = []; % Bar plot of unattached myosin heads 
h_act = [];
h_hsl = [];
h_cmd = line(time(1:index(end)),cmd(1:index(end))/1300,'LineWidth',2,...
    'Parent',h3);


for k = index(1):10:index(end)
    h_dist = line(bins,dist(:,k),'Color',[0 0 1],'Parent',h1);
    h_time = text(-45,0.1,['time: ' num2str(time(k))],'FontSize',16,...
        'FontName','Helvetica','Parent',h1);
% Rate Equations
%     if k == index(1)
%         h_f = line(bins,f/1000,'Color',[1 0 0],'Parent',h1);
%         h_g = line(bins,g/1000,'Color',[0 1 0],'Parent',h1);
%     end
    
%     h_det = bar(1,no_detached(k),'Parent',h_1B,'FaceColor',[0 0 0]);
%     h_act = bar(1,f_act(k),'Parent',h_1C,'FaceColor',[0 0 0]);
    
%     h_cbf = line(time(1:k),cbf(1:k)/10^6,'color',[1 0 0],'Parent',h2);
    h_hsf = line(time(1:k),hsf(1:k)/10^6,'color',[0 1 1],...
        'LineWidth',2,'Color','k','Parent',h2);
%     h_pf = line(time(1:k),pf(1:k)/10^6,'color',[1 0 1],'Parent',h2);
    h_vline2 = line([time(k) time(k)],[min(cbf) max(hsf)],...
        'Color','k','Parent',h2);

%     if k == index(1)
%         h_legend = legend(h2,'X-b force','Half-Sarc. force',...
%             'Pas. force','location','southeast');
%     end
   
    h_vline3 = line(time(k),hsl(k)/1300,...
        'Color','r','Marker','o','MarkerFaceColor',[1 0 0],'Parent',h3);
    h_hsl = line(time(1:k),hsl(1:k)/1300,'color',[0 0 0],'LineWidth',2,...
        'Parent',h3);
    
    
    
    drawnow
%     legend boxoff

  if strcmp(mode,'slow')  
     pause(0.04)     
  end
  
    if k < index(end) - 10
        delete(h_dist)
        delete(h_time)
        delete(h_cbf)
        delete(h_hsf)
        delete(h_pf)
        delete(h_hsl)
        delete(h_vline2)
        delete(h_vline3)
        delete(h_det)
        delete(h_act)
    end
end
end