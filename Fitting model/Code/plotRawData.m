function h = plotRawData(proc_data,info,perts,saveImage)
if iscell(perts)
    newPerts = [];
    for i = 1:size(perts,1)
        for j = 1:size(perts,2)
            newPerts = [newPerts; perts{i,j}];
        end
    end
    perts = newPerts';
end
    

for pert = perts
    if saveImage == 1
        h = figure('visible','off'); 
    elseif saveImage == 0
        h = figure;
    end
    h.Color = 'white';
    hold on
    t = proc_data(pert).time - proc_data(pert).time(1);
    ts = proc_data(pert).spiketimes(2:end) - proc_data(pert).time(1);
    ifr = proc_data(pert).firing_rate;
    
    subplot(6,2,1:2),hold on%, title(info.name,'interpreter','none')
    axis([t(1) t(end) 0 350])
    line(ts,ifr,'Marker','.','LineStyle','none','color',[0 0 0])
    line([proc_data(pert).spiketimes(1)-proc_data(pert).time(1),proc_data(pert).spiketimes(1)-proc_data(pert).time(1)],[0,50],'color','k')
    plot([ts,ts]',[zeros(size(ifr)),ifr]','k'), ylabel('ifr (Hz)','FontSize',8)
    set(gca,'xtick',[],'fontsize',10,'xcolor','white',...
        'ytick',[0 350],'fontName','Helvetica')
    
    subplot(6,2,3:4),hold on
    axis tight
    line(t,proc_data(pert).Force), ylabel('F (N)','FontSize',8)
    set(gca,'xtick',[],'fontsize',10,'xcolor','white','fontName','Helvetica')

    
    subplot(6,2,5:6)
    axis tight
    line(t,proc_data(pert).dFdt), ylabel('dF/dt (N/s)','FontSize',8)
    set(gca,'xtick',[],'fontsize',10,'xcolor','white','fontName','Helvetica')

    
    subplot(6,2,7:8)
    axis tight
    line(t,proc_data(pert).Length), ylabel('L (mm)','FontSize',8)
    set(gca,'xtick',[],'fontsize',10,'xcolor','white','fontName','Helvetica')

    
    subplot(6,2,9:10)
    axis tight
    line(t,proc_data(pert).Velocity), ylabel('V (mm/s)','FontSize',8)
    set(gca,'xtick',[],'fontsize',10,'xcolor','white','fontName','Helvetica')
    
    subplot(6,2,11:12)
    axis tight
    line(t,proc_data(pert).dVdt), ylabel ('A (mm/s^2)','FontSize',8), xlabel('time (s)','FontSize',8)
    set(gca,'fontsize',10,'fontName','Helvetica')
    
    set(h,'Units','Inches');
    set(h,'PaperPosition',[0 0 5 5],'Position',[0 0 5 5])
    align_Ylabels(h)

    

end

end