%% Spindle model - Triangle and RampHold Force
% load(['..' filesep 'data' filesep '2018-08-23' filesep 'triAndRamp 2018-08-23.mat']); % Use this for 
load(['..' filesep 'data' filesep 'TEST.mat']);
time = dataC(1).t;
dataD = dataB(1);
dataS = dataC(1);
L0 = dataS.hs_length(1);

% vertical lines corresponding to xb distributions in figure
% Note: did not use 55e2 or 655e1 in figure (too crowded)
idx_vline = [45e2 50e2 55e2 56e2 57e2 61e2 65e2 655e1 66e2 68e2 1025e1 104e2 11e3 13e3 18e3]/1e3 - 5;


hfig = figure;
set(hfig,'Color','white')
% hfig.RendererMode = 'Painters';

hForce = subplot(2,2,1:2); hold on;
set(hForce,'TickDir','out','FontName','Arial','FontSize',10,...
    'xtick',[],'NextPlot','add')
ylabel('Stress (MPa)')
axis([time(1) time(end) 0 0.2])

hlen = subplot(2,2,3:4); hold on;
set(hlen,'TickDir','out','FontName','Arial','FontSize',10,'NextPlot','add')
ylabel('applied length (L0)'), xlabel('time (s)')
% axis([time(1) time(end) 1300/L0 1400/L0])

lineStyles = linspecer(2);
lineStyles = lineStyles - 0.1;
lineStyles(lineStyles<0) = 0;



h_D = line(time,dataD.hs_force/10^6,'Parent',hForce,'color',lineStyles(1,:),'LineWidth',1.5);
h_S = line(time,dataS.hs_force/10^6,'Parent',hForce,'color',lineStyles(2,:),'LineWidth',1.5);

h_cdl = line(time,dataS.cmd_length/L0,'Parent',hlen,'color','k','LineWidth',1.5);
h_hslS = line(time,dataS.hs_length/L0,'Parent',hlen,'color',lineStyles(2,:),'LineWidth',1.5);
h_hslD = line(time,dataD.hs_length/L0,'Parent',hlen,'color',lineStyles(1,:),'LineWidth',1.5);

for i = 1:numel(idx_vline)
    line([idx_vline(i) idx_vline(i)], [1300/L0 1400/L0],'Parent',hlen,'color','k'), hold on
end
%% Spindle model - General
load(['..' filesep 'data' filesep '2018-08-23' filesep 'DI 2018-08-23.mat']);
idx = 6000:9500;
time = dataS(1).t(idx);
dataY = dataD(1);
dataF = dataS(1);
L0 = dataF.hs_length(1);

hfig = figure;
set(hfig,'Color','White')
% hfig.RendererMode = 'Painters';

haff = subplot(2,2,1:2); hold on;
set(haff,'TickDir','out','FontName','Arial','FontSize',10,...
    'xtick',[],'NextPlot','add')
ylabel('Firing Rate (au)')
axis([time(1) time(end) 0 0.25])

hlen = subplot(2,2,3:4); hold on;
set(hlen,'TickDir','out','FontName','Arial','FontSize',10,'NextPlot','add')
ylabel('Applied Length (L0)'), xlabel('time (s)')
axis([time(1) time(end) 1300/L0 1400/L0])

% Dynamic parameters
kFs_dyn = 0.8;
kFd_dyn = 1.2;
kY_dyn = 0.05;
[r,rs,rd] = sarc2spindle(dataY,dataF,kFs_dyn,kFd_dyn,kY_dyn,1,0);

% Static Parameters
kFs_stat = 1.5;
kFd_stat = 0.8;
kY_stat = 0.01;
[r,rs,rd] = sarc2spindle(dataY,dataF,kFs_stat,kFd_stat,kY_stat,1,0);

lineStyles = linspecer(3);
lineStyles = lineStyles - 0.1;
lineStyles(lineStyles<0) = 0;


h_r = line(time,r(idx),'Parent',haff,'color',lineStyles(1,:),'LineWidth',1.5);
h_sC = line(time,rs(idx),'Parent',haff,'color',lineStyles(2,:),'LineWidth',1.5);
h_dC = line(time,rd(idx),'Parent',haff,'color',lineStyles(3,:),'LineWidth',1.5);

h_cdl = line(time,dataF.cmd_length(idx)/L0,'Parent',hlen,'color','k','LineWidth',1.5);


%% Spindle model - DI
load(['..' filesep 'data' filesep '2018-08-23' filesep 'DI 2018-08-23.mat']);
% load(['..' filesep 'data' filesep 'test.mat']);
time = dataS(1).t;
L0 = dataS(1).hs_length(1);
idx = 6500:8000;


hfig = figure;
set(hfig,'Color','White')


kFs = 1;
kFd = 1;
kY = 0.03;

haff = subplot(2,2,1:2); hold on;
set(haff,'TickDir','out','FontName','Arial','FontSize',10,...
    'xtick',[],'NextPlot','add')
ylabel('Firing Rate (au)')
axis([time(idx(1)) time(idx(end)) 0 0.35])

hlen = subplot(2,2,3:4); hold on;
set(hlen,'TickDir','out','FontName','Arial','FontSize',10,'NextPlot','add')
ylabel('Applied Length (L0)'), xlabel('time (s)')
axis([time(idx(1)) time(idx(end)) 1300/L0 1400/L0])


lineStyles = linspecer(numel(dataS));
lineStyles = lineStyles - 0.1;
lineStyles(lineStyles<0) = 0;


for a = numel(dataS):-1:1
    dataY = dataD(a);
    dataF = dataS(a);
    [r,rs,rd] = sarc2spindle(dataY,dataF,kFs,kFd,kY,1,0.09);
    h_r = line(time(idx),r(idx),'Parent',haff,'Color',lineStyles(a,:),'linewidth',1.5); hold on;
    h_cdl = line(time(idx),dataF.cmd_length(idx)/L0,'Parent',hlen,'Color',lineStyles(a,:),'linewidth',1.5); hold on;
end

load(['..' filesep 'data' filesep '2018-08-23' filesep 'DIandIBdata_20180914.mat']);
hfig2 = figure;
set(hfig2,'Color','White')

hDI = subplot(2,2,1:2); hold on;
set(hDI,'TickDir','out','FontName','Arial','FontSize',10,...
    'NextPlot','add')
plot(vel,DI,'.','MarkerSize',10)

hIB = subplot(2,2,3:4); hold on;
set(hIB,'TickDir','out','FontName','Arial','FontSize',10,'NextPlot','add')
plot(acc,IB,'.','MarkerSize',10)

%% Spindle model - THD
load(['..' filesep 'data' filesep '2018-08-23' filesep 'THD 2018-08-23.mat']);
idx = 4000:17000;
time = dataS(1).t(idx);
L0 = dataS(1).hs_length(1);

hfig = figure;
set(hfig,'Color','White')


kFs = 1;
kFd = 1;
kY = 0.03;

haff = subplot(2,2,1:2); hold on;
set(haff,'TickDir','out','FontName','Arial','FontSize',10,...
    'xtick',[],'NextPlot','add')
ylabel('Firing Rate (au)')
axis([time(1) time(end) 0 0.1])

hlen = subplot(2,2,3:4); hold on;
set(hlen,'TickDir','out','FontName','Arial','FontSize',10,'NextPlot','add')
ylabel('Applied Length (L0)'), xlabel('time (s)')
axis([time(1) time(end) 1300/L0 1400/L0])


lineStyles = linspecer(numel(dataS));
lineStyles = lineStyles - 0.1;
lineStyles(lineStyles<0) = 0;


for a = 1:numel(dataS)
    dataY = dataD(a);
    dataF = dataS(a);
    [r,rs,rd] = sarc2spindle(dataY,dataF,kFs,kFd,kY,1,0.09);
    h_r = line(time,r(idx),'Parent',haff,'Color',lineStyles(a,:),'linewidth',1.5);
    h_cdl = line(time,dataF.cmd_length(idx)/L0,'Parent',hlen,'Color',lineStyles(a,:),'linewidth',1.5);
%     h_sC = line(time,sC,'Parent',h1);
%     h_dC = line(time,dC,'Parent',h1);

end
%% Spindle model - ACT
clear, clc

load(['..' filesep 'data' filesep '2018-08-23' filesep 'ActRH 2018-10-23.mat']);
% load(['..' filesep 'data' filesep 'test4.mat']);

time = dataS(1).t;
L0 = dataS(1).hs_length(1);

hfig = figure;
set(hfig,'Color','White')


kFs = 1.5;
kFd = 0.8;
kY = 0.03;


haff = subplot(2,2,1:2); hold on;
set(haff,'TickDir','out','FontName','Arial','FontSize',10,...
    'xtick',[],'NextPlot','add')
ylabel('Firing Rate (au)')
axis([-1 2 0.0 0.3])

hlen = subplot(2,2,3:4); hold on;
set(hlen,'TickDir','out','FontName','Arial','FontSize',10,'NextPlot','add')
ylabel('Applied Length (L0)'), xlabel('time (s)')
axis([-1 2 1300/L0 1400/L0])

for a = [2 4 6 8]
    for b = [3]
        dataY = dataD(a);
        dataF = dataS(b);
        [r,rs,rd] = sarc2spindle(dataY,dataF,kFs,kFd,kY,1,0.09);
        h_r = line(time,r,'Parent',haff,'Color',[a/8 0 b/8],'LineWidth',1.5);
        h_cdl = line(time,dataF.cmd_length/L0,'Parent',hlen,'Color','k','LineWidth',1.5);
%         h_sC = line(time,rs,'Parent',h1,'Color','c');
%         h_dC = line(time,rd,'Parent',h1,'Color','r');
    end
end

%% Spindle model - AHD
load(['..' filesep 'data' filesep '2018-08-23' filesep 'AHD 2018-08-23.mat']);
time = dataS(1).t;
L0 = dataF.hs_length(1);

hfig = figure;
set(hfig,'Color','White')


kFs = 1.0;
kFd = 1.0;
kY = 0.03;


haff = subplot(2,2,1:2); hold on;
set(haff,'TickDir','out','FontName','Arial','FontSize',10,...
    'xticklabel',[],'NextPlot','add')
ylabel('Firing Rate (au)')
axis([-2 3 0.0 0.1])

hlen = subplot(2,2,3:4); hold on;
set(hlen,'TickDir','out','FontName','Arial','FontSize',10,'NextPlot','add')
ylabel('hs length (nm)'), xlabel('time (s)')
axis([-2 3 1300/L0 1400/L0])

lineStyles = linspecer(numel(dataS));
lineStyles = lineStyles - 0.1;
lineStyles(lineStyles<0) = 0;


% for a = [1 2 3 4 5 8]
for a = 20:-1:1
    dataY = dataD(a);
    dataF = dataS(a);
    [r,rs,rd] = sarc2spindle(dataY,dataF,kFs,kFd,kY,1,0.09);
    h_r = line(time,r,'Parent',haff,'Color',lineStyles(a,:),'linewidth',1.5);
    h_cdl = line(time,dataF.cmd_length/L0,'Parent',hlen,'Color',lineStyles(a,:),'linewidth',1.5);
%     h_sC = line(time,sC,'Parent',h1);
%     h_dC = line(time,dC,'Parent',h1);

end

%% Supplemental - Day sinusoids

% load(['data' filesep 'DaySinusoid.mat']);
load(['..' filesep 'data' filesep 'test4.mat']);
time = dataS(1).t;

hfig = figure;
set(hfig,'Color','White')


kFs = 1.8;
kFd = 2;
kY = 0.15;

% hmus = subplot(3,2,1:2); hold on;
% set(hmus,'TickDir','out','FontName','Arial','FontSize',10,...
%     'xticklabel',[],'NextPlot','add')
% ylabel('Stress (MPa)')
% axis([-1 11.5 0 2.5])

haff = subplot(1,2,1:2); hold on;
set(haff,'TickDir','out','FontName','Arial','FontSize',10,'NextPlot','add')
ylabel('Normalized firing rate')
axis([-1 11.5 -0.5 0.5])

% hlen = subplot(2,2,3:4); hold on;
% set(hlen,'TickDir','out','FontName','Arial','FontSize',10,'NextPlot','add')
% ylabel('hs length (nm)'), xlabel('time (s)')
% axis([-1 11.5 1300 1500])

% for a = [1 2 3 4 5 8]
    dataY = dataD;
    dataF = dataS;
    [r,rs,rd] = sarc2spindle(dataY,dataF,kFs,kFd,kY,0,0.09);
%     h_f = line(time,dataF.hs_force/10^6,'Parent',hmus);
    h_r = line(time,r,'Parent',haff);
    h_cdl = line(time,dataF.cmd_length/1300-1,'Parent',haff);
%     h_sC = line(time,sC,'Parent',h1);
%     h_dC = line(time,dC,'Parent',h1);

%% Supplemental - Haftel RR

load(['data' filesep 'Haftel2004.mat']);
time = data(1).t;

hfig = figure;
set(hfig,'Color','White')


kFs = 1.5;
kFd = 1;
kY = 0.15;

% hmus = subplot(3,2,1:2); hold on;
% set(hmus,'TickDir','out','FontName','Arial','FontSize',10,...
%     'xticklabel',[],'NextPlot','add')
% ylabel('Stress (MPa)')
% axis([-1 11.5 0 2.5])

haff = subplot(1,2,1:2); hold on;
set(haff,'TickDir','out','FontName','Arial','FontSize',10,'NextPlot','add')
ylabel('Normalized firing rate')
axis([-1 11.5 -0.5 0.5])

% hlen = subplot(2,2,3:4); hold on;
% set(hlen,'TickDir','out','FontName','Arial','FontSize',10,'NextPlot','add')
% ylabel('hs length (nm)'), xlabel('time (s)')
% axis([-1 11.5 1300 1500])

% for a = [1 2 3 4 5 8]
    dataY = data;
    dataF = data;
    [r,rs,rd] = sarc2spindle(dataY,dataF,kFs,kFd,kY,0);
%     h_f = line(time,dataF.hs_force/10^6,'Parent',hmus);
    h_r = line(time,r,'Parent',haff);
    h_cdl = line(time,dataF.cmd_length/1300-1,'Parent',haff);
%     h_sC = line(time,sC,'Parent',h1);
%     h_dC = line(time,dC,'Parent',h1);


%% Spindle model - XB Distribution
load(['..' filesep 'data' filesep '2018-08-23' filesep 'THD 2018-08-23.mat']);
hfig = figure;
set(hfig,'Color','White')

%Dynamic
data = dataD(1);
hs = hsD(1);

%Static
data = dataS(1);
hs = hsS(1);

xbl = hs.x_bins;
pops = data.bin_pops(:,1:100:end);

lineStyles = linspecer(2);
lineStyles = lineStyles - 0.1;
lineStyles(lineStyles<0) = 0;

idx = [45e2 50e2 55e2; 56e2 57e2 61e2; 65e2 655e1 66e2; 68e2 1025e1 104e2; 11e3 13e3 18e3];
counter = 0;
for i = 1:5
    for j = 1:3
        counter = counter + 1;
        pops = data.bin_pops(:,idx(i,j));
        subplot(5,3,counter), hold on
        set(gca,'FontName','Arial','FontSize',8)
        axis([-20 20 0 0.2])
        title(num2str(idx(i,j)/1000))
        plot(xbl,pops,'Color',lineStyles(2,:))
    end
end


%% Spindle Model - Alpha-gamma coactivation
% load(['..' filesep 'data' filesep '2018-08-23' filesep 'Alpha-gamma ramp 2018-09-21.mat'])
load(['..' filesep 'data' filesep '2018-08-23' filesep 'alphaGamma2018-10-23.mat'])
idx = 4000:9000;
time = dataS(1).t(idx);
L0 = dataS(1).hs_length(1);

hfig = figure;
set(hfig,'Color','White')

kFs = 1;
kFd = 1;
kY = 0.03;

haff = subplot(3,2,1:2); hold on;
set(haff,'TickDir','out','FontName','Arial','FontSize',10,...
    'xtick',[],'NextPlot','add')
ylabel('Firing Rate (au)')
axis([time(1) time(end) 0 0.25])

hlen = subplot(3,2,3:4); hold on;
set(hlen,'TickDir','out','FontName','Arial','FontSize',10,'NextPlot','add')
ylabel('Applied Length (L0)'), xlabel('time (s)')
axis([time(1) time(end) 1250/L0 1300/L0])

hact = subplot(3,2,5:6); hold on;
set(hact,'TickDir','out','FontName','Arial','FontSize',10,'NextPlot','add')
ylabel('Gamma Activation (normalized)'), xlabel('time (s)')
axis([time(1) time(end) 0 1])

lineStyles = linspecer(numel(dataS));
lineStyles = lineStyles - 0.1;
lineStyles(lineStyles<0) = 0;


for a = 1:numel(dataS)
    dataY = dataD(a);
    dataF = dataS(a);
    [r,rs,rd] = sarc2spindle(dataY,dataF,kFs,kFd,kY,1,0.05);
    h_r = line(time,r(idx),'Parent',haff,'Color',lineStyles(a,:),'linewidth',1.5);
    h_cdl = line(time,dataF.cmd_length(idx)/L0,'Parent',hlen,'Color',lineStyles(a,:),'linewidth',1.5);
    h_act = line(time,dataF.f_activated(idx),'Parent',hact,'Color',lineStyles(a,:),'linewidth',1.5);
%     h_sC = line(time,sC,'Parent',h1);
%     h_dC = line(time,dC,'Parent',h1);

end



