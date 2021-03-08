% Last Modified: 08/16/2017
% Author: KPB
% This script takes processed and analyzed data from the cat muscle spindle
% experiments and generates the MATLAB figures used in the 2017 PLoS
% Computational Biology publication. 

%% Figure 1 Raster Plot
pathname = cd(cd(['..' filesep 'Data' filesep 'proc_data_no_fits' filesep 'aff60_proc']));
affname = 'aff60';
trialID = {'02','10','14','15'};
spindle_plotRaster(0.02,trialID,pathname,affname);

load([pathname filesep affname '_mn0_02_accel_ramp_proc.mat'])
fig;
hold on
subplot(3,1,1),hold on
pert = 1;
plot(proc_data(pert).spiketimes(2:end) - proc_data(pert).time(1),proc_data(pert).firing_rate,'.')
plot([proc_data(pert).spiketimes - proc_data(pert).time(1), proc_data(pert).spiketimes - proc_data(pert).time(1)],[0,20],'k')
subplot(3,1,2)
plot(proc_data(pert).time - proc_data(pert).time(1),proc_data(pert).Force)
subplot(3,1,3)
plot(proc_data(pert).time - proc_data(pert).time(1),proc_data(pert).dFdt)


%% Figure 3: Model Fits

% 3B: Effects of Varying Velocity
datapath = ['..' filesep 'Data' filesep 'proc_data_no_buffer'];
filenames = {[datapath filesep 'aff32_proc' filesep 'aff32_mn0_58_trad_series_proc.mat']};%...
ifrdatapath = ['..' filesep 'Data' filesep 'proc_data_no_fits'];
ifrfilenames = {[ifrdatapath filesep 'aff32_proc' filesep 'aff32_mn0_58_trad_series_proc.mat']};%...
perts = [9,13,15,16];
for i = 1:size(perts,1)
    spindle_plotFitData(filenames{i},ifrfilenames{i},perts(i,:))
end
    


% 3C: Effects of Varying Acceleration
datapath = ['..' filesep 'Data' filesep 'proc_data_no_buffer'];
filenames = {[datapath filesep 'aff25_proc' filesep 'aff25_mn0_04_accel_series_proc.mat']};
perts = [11,16];
ifrdatapath = ['..' filesep 'Data' filesep 'proc_data_no_fits'];
ifrfilenames = {[ifrdatapath filesep 'aff25_proc' filesep 'aff25_mn0_04_accel_series_proc.mat']};

for i = 1:size(perts,1)
    spindle_plotFitData(filenames{i},ifrfilenames{i},perts(i,:))
end

% 3D: Ramp & Release Perturbations
datapath = ['..' filesep 'Data' filesep 'proc_data_no_buffer'];
filenames = {[datapath filesep 'aff32_proc' filesep 'aff32_mn0_25_triangle_series_proc.mat']...
    [datapath filesep 'aff32_proc' filesep 'aff32_mn0_28_triangle_series_proc.mat']};
perts = [2;2];
ifrdatapath = ['..' filesep 'Data' filesep 'proc_data_no_fits'];
ifrfilenames = {[ifrdatapath filesep 'aff32_proc' filesep 'aff32_mn0_25_triangle_series_proc.mat']...
    [ifrdatapath filesep 'aff32_proc' filesep 'aff32_mn0_28_triangle_series_proc.mat']};

for i = 1:size(perts,1)
    spindle_plotFitData(filenames{i},ifrfilenames{i},perts(i,:))
end

%% Figure 4: Tendon Compliance
datapath = ['..' filesep 'Data' filesep 'proc_data_no_buffer'];
filenames = {[datapath filesep 'aff32_proc' filesep 'aff32_mn0_02_accel_ramp_proc.mat']...
    [datapath filesep 'aff32_proc' filesep 'aff32_mn0_10_triangle_series_proc.mat']};
perts = [1;1];
RHdata = load(filenames{1});
RRdata = load(filenames{2});

RHtime = RHdata.proc_data(perts(1)).time - RHdata.proc_data(perts(1)).time(1);
RRtime = RRdata.proc_data(perts(2)).time - RRdata.proc_data(perts(1)).time(1);
RHForce = RHdata.proc_data(perts(1)).Force/10;
RRForce = RRdata.proc_data(perts(2)).Force/10;
RHLength = RHdata.proc_data(perts(1)).Length;
RRLength = RRdata.proc_data(perts(2)).Length;
RHLengthF = RHdata.proc_data(perts(1)).LengthF;
RRLengthF = RRdata.proc_data(perts(2)).LengthF;
RHLengthF_hc = RHdata.proc_data(perts(1)).LengthF_hc;
RRLengthF_hc = RRdata.proc_data(perts(2)).LengthF_hc;
RHVelocity = RHdata.proc_data(perts(1)).Velocity;
RRVelocity = RRdata.proc_data(perts(2)).Velocity;
RHVelocityF = RHdata.proc_data(perts(1)).VelocityF;
RRVelocityF = RRdata.proc_data(perts(2)).VelocityF;
RHVelocityF_hc = RHdata.proc_data(perts(1)).VelocityF_hc;
RRVelocityF_hc = RRdata.proc_data(perts(2)).VelocityF_hc;
RHdVdt = RHdata.proc_data(perts(1)).dVdt;
RRdVdt = RRdata.proc_data(perts(2)).dVdt;
RHdVFdt = RHdata.proc_data(perts(1)).dVFdt;
RRdVFdt = RRdata.proc_data(perts(2)).dVFdt;
RHdVFdt_hc = RHdata.proc_data(perts(1)).dVFdt_hc;
RRdVFdt_hc = RRdata.proc_data(perts(2)).dVFdt_hc;
RHspiketimes = RHdata.proc_data(perts(1)).spiketimes(2:end) - RHdata.proc_data(perts(1)).time(1);
RHIFR = RHdata.proc_data(perts(1)).firing_rate;
RRspiketimes = RRdata.proc_data(perts(2)).spiketimes(2:end) - RRdata.proc_data(perts(2)).time(1);
RRIFR = RRdata.proc_data(perts(2)).firing_rate;

RHxlim = [0 RHtime(end)];
RRxlim = [0 RRtime(end)];

fig, hold on

subplot(10,4,[1:2,5:6]), hold on 
axis([RHxlim 0 max(RHdata.proc_data(perts(1)).Force/10)])
set(gca,'xtick',[],'clipping','off','FontSize',12,'FontName','Helvetica')
line(RHtime,RHForce)

subplot(10,4,[9:10,13:14]), hold on
axis([RHxlim 0 max(RHdata.proc_data(perts(1)).Length)])
set(gca,'xtick',[],'clipping','off','FontSize',12,'FontName','Helvetica')
line(RHtime,RHLength,'color',[0 0 0])
line(RHtime,RHLengthF,'color',[0.25 0.25 0.25])
line(RHtime,RHLengthF_hc,'color',[0.7 0.7 0.7])

subplot(10,4,[17:18,21:22]), hold on
axis([RHxlim -max(RHdata.proc_data(perts(1)).VelocityF_hc)...
    max(RHdata.proc_data(perts(1)).VelocityF_hc)])
set(gca,'xtick',[],'clipping','off','FontSize',12,'FontName','Helvetica')
line(RHtime,RHVelocity,'color',[0 0 0])
line(RHtime,RHVelocityF,'color',[0.25 0.25 0.25])
line(RHtime,RHVelocityF_hc,'color',[0.7 0.7 0.7])

subplot(10,4,[25:26,29:30]), hold on
axis([RHxlim -max(RHdata.proc_data(perts(1)).dVFdt_hc)...
    max(RHdata.proc_data(perts(1)).dVFdt_hc)])
set(gca,'clipping','off','FontSize',12,'FontName','Helvetica')
line(RHtime,RHdVdt,'color',[0 0 0])
line(RHtime,RHdVFdt,'color',[0.25 0.25 0.25])
line(RHtime,RHdVFdt_hc,'color',[0.7 0.7 0.7])

subplot(10,4,[33:34,37:38]),hold on
axis([RHxlim 0 ...
    max(RHdata.proc_data(perts(1)).firing_rate)])
set(gca,'clipping','off','FontSize',12,'FontName','Helvetica')
plot(RHspiketimes,RHIFR,'.')

subplot(10,4,[3:4,7:8]), hold on 
axis([RHxlim 0 max(RHdata.proc_data(perts(1)).Force/10)])
set(gca,'xtick',[],'ytick',[],'clipping','off','FontSize',12,'FontName','Helvetica')
line(RRtime,RRForce)

subplot(10,4,[11:12,15:16]), hold on
axis([RHxlim 0 max(RHdata.proc_data(perts(1)).Length)])
set(gca,'xtick',[],'ytick',[],'clipping','off','FontSize',12,'FontName','Helvetica')
line(RRtime,RRLength,'color',[0 0 0])
line(RRtime,RRLengthF,'color',[0.25 0.25 0.25])
line(RRtime,RRLengthF_hc,'color',[0.7 0.7 0.7])

subplot(10,4,[19:20,23:24]), hold on
axis([RHxlim -max(RHdata.proc_data(perts(1)).VelocityF_hc)...
    max(RHdata.proc_data(perts(1)).VelocityF_hc)])
set(gca,'xtick',[],'ytick',[],'clipping','off','FontSize',12,'FontName','Helvetica')
line(RRtime,RRVelocity,'color',[0 0 0])
line(RRtime,RRVelocityF,'color',[0.25 0.25 0.25])
line(RRtime,RRVelocityF_hc,'color',[0.7 0.7 0.7])

subplot(10,4,[27:28,31:32]), hold on
axis([RHxlim -max(RHdata.proc_data(perts(1)).dVFdt_hc)...
    max(RHdata.proc_data(perts(1)).dVFdt_hc)])
set(gca,'ytick',[],'clipping','off','FontSize',12,'FontName','Helvetica')
line(RRtime,RRdVdt,'color',[0 0 0])
line(RRtime,RRdVFdt,'color',[0.25 0.25 0.25])
line(RRtime,RRdVFdt_hc,'color',[0.7 0.7 0.7])

subplot(10,4,[35:36,39:40]),hold on
axis([RHxlim 0 ...
    max(RHdata.proc_data(perts(1)).firing_rate)])
set(gca,'ytick',[],'clipping','off','FontSize',12,'FontName','Helvetica')
plot(RRspiketimes,RRIFR,'.')

%% Figure 5: DI & MLEs for Force model - individual and training datasets
% data = readtable(['..' filesep 'Data' filesep 'stats' ...
%     filesep 'all_affs_stat_array_no_buffer_lags2.csv']);
data = readtable(['..' filesep 'Data' filesep 'stats' ...
    filesep 'all_affs_stat_array_no_buffer_All.csv']);
% load(['..' filesep 'Data' filesep 'stats' ...
%     filesep 'MLEtable_force']);
load(['..' filesep 'Data' filesep 'stats' ...
    filesep 'xvalDataAll20-JUL-2016.mat']);

% MLEtable_force = generateMLETable(xvalDataAll);

model = data.model;

maxV = data.maxV(model==1);
maxL = data.maxL(model==1);



% Dynamic Index
afferents = [24,25,26,27,32,56,57,58,60,64];
cond_Vels = [97.50, 81.25, 90.28, 91.98, 84.78, 82.67, 87.89, 92.78, 75.00, 80.10]; % These are from Lena's notebook
model = data.model;
affnum = data.affnum(model==1);
DI = data.DI_1s(model==1);    %DI = dynamic response - response 1s later
DIc = data.DI_500ms(model==1);%Classic DI = dynamic response - response 0.5s later
bgF = data.bgF(model==1);
kF = data.kF(model==1)/10; %This fixes a conversion error
kdF = data.kdF(model==1);
pert_type = data.pert_type(model==1);


% plotVel = 20; %Velocity in mm/s for which to plot dynamic index
plotVelInd = [1:5];
velbinsize = 5;
numBins = 5;

n = zeros(numBins,length(afferents));


plotDI = struct('data',cell(1,20),'grp',cell(1,20),'newgrp',cell(1,20),'newdata',cell(1,20));
plotDIc = struct('data',cell(1,20),'grp',cell(1,20),'newgrp',cell(1,20),'newdata',cell(1,20));
plotkF = struct('data',cell(1,20),'newdata',cell(1,20));
plotkdF = struct('data',cell(1,20),'newdata',cell(1,20));

Vindex = [];
for aff = 1:10;%:length(afferents)
    index = affnum==afferents(aff) & ~strcmp(pert_type,'triangle_series') & maxL >3; %Trial characteristics to include
    
    for j = 1:numBins; %velocity bins
        if j < 10
            Vindex = maxV>= 0.5*velbinsize*((2*j-1)-1) & maxV<0.5*velbinsize*((2*j-1)-1)+velbinsize & index;
        else
            Vindex = maxV>= 0.5*velbinsize*((2*j-1)-1) & index;
        end
    

        n(j,aff) = numel(DI(Vindex)); % Number of datapoints for this afferent/velocity bin
        affplotDI(j).data = DI(Vindex); %DI values for this afferent/velocity bin
        meanplotDI(aff,j) = mean(affplotDI(j).data); %mean DI value for this afferent/velocity bin
        plotDI(j).data = [plotDI(j).data; affplotDI(j).data]; %Concatenate this round's DI values with all values
        affplotDI(j).grp = aff*ones(size(affplotDI(j).data)); %afferent indexing for this afferent
        plotDI(j).grp = [plotDI(j).grp; affplotDI(j).grp]; %Concatenate this rounds's afferent indexing
        
        affplotDIc(j).data = DIc(Vindex); %DIc values for this afferent/velocity bin
        meanplotDIc(aff,j) = mean(affplotDIc(j).data); %mean DIc value for this afferent/velocity bin
        plotDIc(j).data = [plotDIc(j).data; affplotDIc(j).data]; %Concatenate this round's DIc values with all values
        affplotDIc(j).grp = aff*ones(size(affplotDIc(j).data)); %afferent indexing for this afferent
        plotDIc(j).grp = [plotDIc(j).grp; affplotDIc(j).grp]; %Concatenate this rounds's afferent indexing
        
        affplotkF(j).data = kF(Vindex);
        affplotkdF(j).data = kdF(Vindex);
        plotkF(j).data = [plotkF(j).data; affplotkF(j).data];
        plotkdF(j).data = [plotkdF(j).data; affplotkdF(j).data];
    end
end

clearvars affplotDI affplotDIc affplotkF affplotkdF

% Restructure in descending order of DIc
[sorted,plotOrder] = sort(meanplotDIc(:,plotVelInd(3)),1,'descend');

for i = 1:length(plotOrder)
    for j = 1:numBins
    plotDI(j).newgrp = [plotDI(j).newgrp; i*ones(size(plotDI(j).grp(plotDI(j).grp==plotOrder(i))))];
    plotDI(j).newdata = [plotDI(j).newdata; plotDI(j).data(plotDI(j).grp==plotOrder(i))];
    plotDIc(j).newdata = [plotDIc(j).newdata; plotDIc(j).data(plotDI(j).grp==plotOrder(i))];
    plotkF(j).newdata = [plotkF(j).newdata; plotkF(j).data(plotDI(j).grp==plotOrder(i))];
    plotkdF(j).newdata = [plotkdF(j).newdata; plotkdF(j).data(plotDI(j).grp==plotOrder(i))];
    end
end
fig; hold on;
ylabels = {'D.I. - 1s (imp/s)', 'D.I. - 0.5s (imp/s)', 'k_F (N/imp/s)',...
    'k_{dF} (N/s/imp/s)','D.I. 20mm/s','k_{F} (N/imp/s)','k_{dF} (N/s/imp/s)'};

counter = 0;

for sp = [1 3 5 7]
    counter = counter + 1;
    subplot(7,2,sp:sp+1), hold on
    set(gca,'FontSize',12,'FontName','Helvetica')
    ylabel(ylabels{counter})
    
    for i = 1:length(plotVelInd)
        if sp == 1
            notBoxPlot(plotDI(plotVelInd(i)).newdata,plotDI(plotVelInd(i)).newgrp+(0.15*i-0.15),0.05,'patch');
        elseif sp == 3
            notBoxPlot(plotDIc(plotVelInd(i)).newdata,plotDI(plotVelInd(i)).newgrp+(0.15*i-0.15),0.05,'patch');
        elseif sp == 5
            %Scale by 10 because force was scaled wrong
            notBoxPlot(plotkF(plotVelInd(i)).newdata*10,plotDI(plotVelInd(i)).newgrp+(0.15*i-0.15),0.05,'patch');
        elseif sp == 7
            notBoxPlot(plotkdF(plotVelInd(i)).newdata,plotDI(plotVelInd(i)).newgrp+(0.15*i-0.15),0.05,'patch');
        end
    end
    
    xtickpos = [];
    for i = 1:10
        for j = 1:5
            xtickpos(end+1) = i + (0.15*j-0.15);
        end
    end
    set(gca,'xtick',xtickpos,'xticklabel',[])
end



subplot(7,2,9:10), hold on
set(gca,'FontSize',12,'FontName','Helvetica')
notBoxPlot(plotDI(plotVelInd(3)).newdata,plotDI(plotVelInd(3)).newgrp);
set(gca,'xtick',[])

paramNames = {'k_F','k_{dF}'};

MLEdata = zeros(100,10);

for param = 1:2 %Force gain and dF/dt gain
subplot(7,2,param*2+9:param*2+10),hold on
set(gca,'FontSize',12,'FontName','Helvetica')
ylabel(ylabels{param+5})


for aff = 1:10
    MLEdata(:,aff) = xvalDataAll.MLEout(1,plotOrder(aff)).MLEs.all(:,param);
end

if param == 1
    %Scale force gain by 10 because of earlier conversion error
    MLEdata = MLEdata*10;
end
notBoxPlot(MLEdata)

if param == 1
set(gca,'xtick',[],'xticklabel',[]), hold on
end
end
xlabel('Afferent number')
set(gca,'xticklabel',plotOrder,'FontSize',12,'FontName','Helvetica')

%% Figure 6: Cross-validation/ prediction results
load(['..' filesep 'Data' filesep 'stats' ...
    filesep 'xvalDataAll20-Jul-2016.mat'])
data = readtable(['..' filesep 'Data' filesep 'stats' ...
    filesep 'all_affs_stat_array_no_buffer_no_classic_model.csv']);

AICc_L = xvalDataAll.test.AICcdata.L';
AICc_w = xvalDataAll.test.AICcdata.w';
BIC_L = xvalDataAll.test.BICdata.L';
BIC_w = xvalDataAll.test.BICdata.w';
R2_mean = xvalDataAll.test.R2data.mean';
R2_stdev = xvalDataAll.test.R2data.stdev';
AICc_mean = xvalDataAll.test.AICcdata.mean';
AICc_stdev = xvalDataAll.test.AICcdata.stdev';
BIC_mean = xvalDataAll.test.BICdata.mean';
BIC_stdev = xvalDataAll.test.BICdata.stdev';
numGroups = size(R2_mean,1);%Number of groups of bars (same as number of afferents)
numBars = size(R2_mean,2); %Number of bars in each group
groupWidth = min(0.8, numBars/(numBars + 1.5)); %This is how Matlab determines the width of the group


fig, hold on
subplot(14,1,1:2), hold on
bar(R2_mean), ylabel('R^2')
set(gca,'box','off')
for i = 1:numBars
      % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
      x = (1:numGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*numBars);  % Aligning error bar with individual bar
      h = errorbar(x, R2_mean(:,i), R2_stdev(:,i));
      h.LineStyle = 'none';
      h.Color = [0 0 0];
end
set(gca,'xtick',[])

subplot(14,1,3:4), hold on
bar(AICc_mean), ylabel('mean AICc')
set(gca,'box','off','xtick',[])
for i = 1:numBars
      % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
      x = (1:numGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*numBars);  % Aligning error bar with individual bar
      h = errorbar(x, AICc_mean(:,i), AICc_stdev(:,i));
      h.LineStyle = 'none';
      h.Color = [0 0 0];
end


subplot(14,1,5:6)
bar(AICc_L), ylabel('AICc L')
set(gca,'box','off','xtick',[])


subplot(14,1,7:8)
bar(AICc_w), ylabel('AICc w')
set(gca,'box','off','xtick',[])


subplot(14,1,9:10), hold on
bar(BIC_mean), ylabel('mean BIC')
set(gca,'box','off','xtick',[])
for i = 1:numBars
      % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
      x = (1:numGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*numBars);  % Aligning error bar with individual bar
      h = errorbar(x, BIC_mean(:,i), BIC_stdev(:,i));
      h.LineStyle = 'none';
      h.Color = [0 0 0];
end


subplot(14,1,11:12)
bar(BIC_L), ylabel('BIC L')
set(gca,'box','off','xtick',[])


subplot(14,1,13:14)
bar(BIC_w), ylabel('BIC w')
set(gca,'box','off')


popAICc_mean = mean(AICc_mean);
popAICc_stdev = sqrt(sum(AICc_stdev.^2)/10);
popR2_mean = mean(R2_mean);
popR2_stdev = sqrt(sum(R2_stdev.^2)/10);
popAICc_delta = popAICc_mean - min(popAICc_mean);
popAICc_L = exp(-1/2*popAICc_delta);
popAICc_w = popAICc_L/sum(popAICc_L);


fig, hold on
subplot(8,1,1:2), hold on
bar(popR2_mean)
set(gca,'box','off','xtick',[])
h = errorbar(popR2_mean,popR2_stdev);
h.LineStyle = 'none';

subplot(8,1,3:4), hold on
bar(popAICc_mean)
set(gca,'box','off','xtick',[])
h = errorbar(popAICc_mean,popAICc_stdev);
h.LineStyle = 'none';

subplot(8,1,5:6), hold on
bar(popAICc_L)
set(gca,'box','off','xtick',[])
h.LineStyle = 'none';

subplot(8,1,7:8), hold on
bar(popAICc_w)
set(gca,'box','off','xtick',[])
h.LineStyle = 'none';




