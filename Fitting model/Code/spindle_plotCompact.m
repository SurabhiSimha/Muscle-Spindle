% This function takes all relevant data from spindle_operate and plots it in
% a single figure for each perturbation. Data plotted includes the trial's spike train, firing
% rate, Force component of optimal fit (k_F*Force), dF/dt component of
% optimal fit (k_dF*dFdt), static component of optimal fit, as well as the
% optimal fit itself. Axes labels, legend, and additional text (e.g.
% feedback gains, r-squared, etc.) are included underneath the figure title.
%
% Modifications:
% 1) Modified to take inputs 'affname' (e.g.
% 'aff24') and load all data for the specified afferent.
%
% Created 04/2013
% Last Modified: 06/2015
% Author: Kyle Blum
function spindle_plotCompact
tic
% affs = [38,51,61];
%affs = {'23','24','25','26','27','32','38','51','56','57','58','60','64'};
affs = {'60'};

for affnum = 1:length(affs);
    affname = ['aff' num2str(affs{affnum})]
    pathname = cd(cd(['..' filesep 'Data' filesep 'proc_data' filesep affname '_proc' filesep]));
    
    directory = dir(pathname);
    for m = 1:length(directory) % For all files in directory
        filename = directory(m).name % Current file
        if strfind(filename,'proc.mat')    % If current file is .mat file
            loadname = [pathname filesep filename];
            load(loadname)
            for model = 1:size(fit_data,2); % For each model
                if length(proc_data) <= 16;
                    xtrans = [0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3]*4-1.8;    % Horiz. Translation for 16 plots
                    ytrans = [2 2 2 2 1 1 1 1 0 0 0 0 -1 -1 -1 -1]*500; % Vert. Translation for 16 plots
                else
                    disp([filename ' does not work'])
                    break
                end
                figure()
                title([filename ' Model: ' num2str(model)],...
                    'Interpreter','none','fontsize',10);
                set(gcf,'PaperUnits','normalized','PaperPosition',[0 0 1 1])
                set(gcf,'visible','off')
                hold on
                for pert = 1:length(proc_data);
                    %%%%% FIRST,SECOND,THIRD MODEL START %%%%%
                    if model == 1; % If we are plotting 1st model
                        
                        %%% Fit Data %%%
                        opt_gains = fit_data(model).opt_gains;                  % Optimal gains from fitting algorithm
                        %                             if model == 1;
                        Force = proc_data(pert).Force;
                        dFdt = proc_data(pert).dFdt;                                        % dFdt
                        time = proc_data(pert).time-proc_data(pert).time(1);          % Time relative to start of pert.
                        FR_times = fit_data(model).FR(pert).times;        % Times at which firing rate is defined
                        FR_recorded = fit_data(model).FR(pert).values;    % Values for FR from ISI
                        
                        fit_comps.Force_comp = opt_gains(pert,1)*Force;                               % Force component
                        fit_comps.dFdt_comp = opt_gains(pert,2)*dFdt;                                 % dFdt component
                        fit_comps.Static_comp = opt_gains(pert,3)*ones(size(fit_comps.Force_comp));             % Reconstructed static component
                        opt_fit = fit_comps.Force_comp + fit_comps.dFdt_comp + fit_comps.Static_comp;                     %
                        
                        
                    elseif model == 2;
                        
                        %%% Kinematics Data %%%
                        Length = proc_data(pert).Length;     % 2000 Hz Length
                        Velocity = proc_data(pert).Velocity; % 2000 Hz Velocity
                        Accel = proc_data(pert).dVdt;
                        
                        %%% Fit Data %%%
                        opt_gains = fit_data(model).opt_gains;                  % Optimal gains from fitting algorithm
                        fit_comps.Length_comp = opt_gains(pert,1)*Length;
                        fit_comps.Velocity_comp = opt_gains(pert,2)*Velocity;
                        fit_comps.Accel_comp = opt_gains(pert,3)*Accel;
                        fit_comps.Static_comp = opt_gains(pert,4)*ones(size(Length));
                        opt_fit = fit_comps.Length_comp + fit_comps.Velocity_comp + fit_comps.Accel_comp + fit_comps.Static_comp;
                        
                        %%% Timing %%%
                        time = proc_data(pert).time-proc_data(pert).time(1);          % Time relative to start of pert.
                        FR_times = fit_data(model).FR(pert).times;        % Times at which firing rate is defined
                        FR_recorded = fit_data(model).FR(pert).values;    % Values for FR from ISI
                        
                    elseif model == 3;
                        
                        %%% Kinematics Data %%%
                        Length = proc_data(pert).Length;     % 2000 Hz Length
                        Velocity = proc_data(pert).Velocity; % 2000 Hz Velocity
                        Accel = proc_data(pert).dVdt;
                        Force = proc_data(pert).Force;
                        dFdt = proc_data(pert).dFdt;                                        % dFdt
                        
                        %%% Fit Data %%%
                        opt_gains = fit_data(model).opt_gains;                  % Optimal gains from fitting algorithm
                        fit_comps.Length_comp = opt_gains(pert,3)*Length;
                        fit_comps.Velocity_comp = opt_gains(pert,4)*Velocity;
                        fit_comps.Accel_comp = opt_gains(pert,5)*Accel;
                        fit_comps.Static_comp = opt_gains(pert,6)*ones(size(Length));
                        fit_comps.Force_comp = opt_gains(pert,1)*Force;                               % Force component
                        fit_comps.dFdt_comp = opt_gains(pert,2)*dFdt;                                 % dFdt component
                        opt_fit = fit_comps.Force_comp + fit_comps.dFdt_comp + fit_comps.Length_comp + fit_comps.Velocity_comp + fit_comps.Accel_comp + fit_comps.Static_comp;
                        
                        %%% Timing %%%
                        time = proc_data(pert).time-proc_data(pert).time(1);          % Time relative to start of pert.
                        FR_times = fit_data(model).FR(pert).times;        % Times at which firing rate is defined
                        FR_recorded = fit_data(model).FR(pert).values;    % Values for FR from ISI
                        
                    end
                    
                    %%% Creat graphics objects %%%
                    ymin = -100;
                    ymax = 350;
                    xlim([0 3].*4)
                    ylim([ymin ymax].*4)
                    xmin = time(1);
                    xmax = time(end);
                    R_squared = fit_data(model).R_squared(pert);         % R-squared value for trial
                    
                    h = plotcode(xmin,xmax,ymin,ymax,time,fit_comps,opt_fit,FR_times,FR_recorded,R_squared,pert,model);
                    transobj(h,xtrans(pert),ytrans(pert))
                    
                    axis off
                end
            end
            %%% Save Figure %%%
            savepath = [pathname filesep...                   % Path in which to save file
                'Plots' filesep...
                'Compact'];
            if ~exist(savepath,'dir')                         % If savepath doesn't exist...
                mkdir(savepath)                               % make it.
            end
            savename = [savepath filesep affname];                     %         *
            print(gcf,savename,'-dpsc','-append')      % Append figure to PostScript color file
            close
        end
    end
end


toc
end

function h = plotcode(xmin,xmax,ymin,ymax,time,fit_comps,opt_fit,FR_times,FR_recorded,R_squared,pert,model)
h = [];
h(end+1) = line([xmin xmax],[ymin ymin],'Color',[0 0 0],'Clipping','off'); % x-axis
h(end+1) = line([xmin xmin],[ymin ymax],'Color',[0 0 0],'Clipping','off'); % y-axis
h(end+1) = line(time,opt_fit,'LineStyle','-','Color',[0.85 0 0.4],'Clipping','off','LineWidth',1); % Optimal Fit
h(end+1) = line(FR_times,FR_recorded,'LineStyle','None','Marker','.','Color',[0 0 0],'Clipping','off','MarkerSize',3); % Recorded FR
% h(end+1) = line([0.5 0.5],[ymin ymax],'Color',[0 0 0]); %Onset of perturbation
h(end+1) = text(xmax/2+0.08*xmax,ymax-110,...                          % R-squared text
    ['R^2: ' num2str(R_squared,3)],...
    'fontsize',8,...
    'HorizontalAlignment','center',...
    'Interpreter','none');
if model == 1;
    h(end+1) = line(time,fit_comps.Force_comp,'LineStyle','-','Color',[1 0 0],'Clipping','off','LineWidth',0.1);
    h(end+1) = line(time,fit_comps.dFdt_comp,'LineStyle','-','Color',[0 1 0],'Clipping','off','LineWidth',0.1);
    h(end+1) = line(time,fit_comps.Static_comp,'LineStyle','-','Color',[1 1 1],'Clipping','off','LineWidth',0.1);
elseif model == 2;
    h(end+1) = line(time,fit_comps.Length_comp,'LineStyle','-','Color',[1 1 0],'Clipping','off','LineWidth',0.1);
    h(end+1) = line(time,fit_comps.Velocity_comp,'LineStyle','-','Color',[0 1 1],'Clipping','off','LineWidth',0.1);
    h(end+1) = line(time,fit_comps.Accel_comp,'LineStyle','-','Color',[1 0 0],'Clipping','off','LineWidth',0.1);
    h(end+1) = line(time,fit_comps.Static_comp,'LineStyle','-','Color',[1 1 1],'Clipping','off','LineWidth',0.1);
elseif model == 3;
    h(end+1) = line(time,fit_comps.Length_comp,'LineStyle','-','Color',[1 1 0],'Clipping','off','LineWidth',0.1);
    h(end+1) = line(time,fit_comps.Velocity_comp,'LineStyle','-','Color',[0 1 1],'Clipping','off','LineWidth',0.1);
    h(end+1) = line(time,fit_comps.Accel_comp,'LineStyle','-','Color',[1 0 0],'Clipping','off','LineWidth',0.1);
    h(end+1) = line(time,fit_comps.Force_comp,'LineStyle','-','Color',[0 0 1],'Clipping','off','LineWidth',0.1);
    h(end+1) = line(time,fit_comps.dFdt_comp,'LineStyle','-','Color',[0 1 0],'Clipping','off','LineWidth',0.1);
    h(end+1) = line(time,fit_comps.Static_comp,'LineStyle','-','Color',[1 1 1],'Clipping','off','LineWidth',0.1);
end
end


function []=transobj(h,xtrans,ytrans)
% function []=transobj(h,xtrans,ytrans)
% This functioon translates the graphics in the vector of handles h by the
% translations (in units of the current figure) in xtrans and ytrans.
h=h(:);
for hand=1:length(h)
    try
        temp=get(h(hand),'XData');
        set(h(hand),'XData',temp+xtrans);
        temp=get(h(hand),'YData');
        set(h(hand),'YData',temp+ytrans);
    catch
        temp=get(h(hand),'Position');
        set(h(hand),'Position',temp+[xtrans ytrans 0]);
    end
end
end