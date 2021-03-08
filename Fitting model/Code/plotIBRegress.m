function h = plotIBRegress(ibData)
h = figure;
h.Color = 'white';
hold on


subplot(5,1,1:2), hold on
axis tight
plot(ibData.dFdt.pk,ibData.IFR.pk,'.')
dF_xdata = min(ibData.dFdt.pk):0.01:max(ibData.dFdt.pk);
dF_ydata = ibData.regressData.dFdt.b(2)*dF_xdata + ibData.regressData.dFdt.b(1);
line(dF_xdata,dF_ydata)
text(min(ibData.dFdt.pk)+10,max(ibData.IFR.pk)-10,...
    ['R^2 = ' num2str(ibData.regressData.dFdt.Rsq)],'FontSize',12,...
    'FontName','Helvetica')
set(gca,'fontsize',10,'fontName','Helvetica')


subplot(5,1,4:5), hold on
axis tight
plot(ibData.acc.pk,ibData.IFR.pk,'.')
acc_xdata = min(ibData.acc.pk):0.01:max(ibData.acc.pk);
acc_ydata = ibData.regressData.acc.b(2)*acc_xdata + ibData.regressData.acc.b(1);
line(acc_xdata,acc_ydata)
text(min(ibData.acc.pk)+100,max(ibData.IFR.pk)-10,...
    ['R^2 = ' num2str(ibData.regressData.acc.Rsq)],'FontSize',12,...
    'FontName','Helvetica')
set(gca,'fontsize',10,'fontName','Helvetica')

set(h,'Units','Inches');
set(h,'PaperPosition',[0 0 5 2],'Position',[0 0 5 2])
    
end
