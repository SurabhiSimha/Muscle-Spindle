% This script takes all relevant data from spindle_operate and plots it in
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
% Created 4/2013
% Last Modified: 04/2015
% Author: Kyle Blum

tic
% affs = [38,51,61];
 %affs = {'23','24','25','26','27','32','38','51','56','57','58','60','64'};         
affs = {'60'};

for affnum = 1:length(affs);
    affname = ['aff' num2str(affs{affnum})]
      pathname = cd(cd(['..' filesep 'Data' filesep 'proc_data' filesep affname '_proc' filesep]));
%       pathname = cd(cd(['..' filesep 'Rats' filesep 'Data' filesep 'proc_data' filesep affname '_proc' filesep]));

    directory = dir(pathname);
    for m = 1:length(directory) % For all files in directory
        filename = directory(m).name % Current file
        if strfind(filename,'proc.mat')    % If current file is .mat file
            loadname = [pathname filesep filename];
            load(loadname)
            for model = 1;%3;%1:size(fit_data,2) % For each model
                for pert = 1:length(proc_data); % For each perturbation
                    
                    %%%%% FORCE MODELS START %%%%%
                    if model == 1; % If we are plotting 1st model
                        
                        %%% Fit Data %%%
                        opt_gains = fit_data(model).opt_gains;                  % Optimal gains from fitting algorithm
                        Force = proc_data(pert).Force;
                        dFdt = proc_data(pert).dFdt;                                        % dFdt
                        Length = proc_data(pert).Length;     % 2000 Hz Length
                        Velocity = proc_data(pert).Velocity; % 2000 Hz Velocity
                        Accel = proc_data(pert).dVdt;        % 2000 Hz Acceleration
                        time = proc_data(pert).time-proc_data(pert).time(1);          % Time relative to start of pert.
                        FR_times = fit_data(model).FR(pert).times;        % Times at which firing rate is defined
                        FR_recorded = fit_data(model).FR(pert).values;    % Values for FR from ISI
                        
                        
                        Force_comp = opt_gains(pert,1)*(Force) + opt_gains(pert,4);         % Force component
                        dFdt_comp = opt_gains(pert,2)*dFdt;                                 % dFdt component
                        dFdt_comp(dFdt_comp<0) = 0;
                        static_comp = opt_gains(pert,3)*ones(size(Force_comp));             % Reconstructed static component
                        opt_fit = Force_comp + dFdt_comp + static_comp;                     %
                        
                        
                        %%% Fix potential rounding error %%%
                        if length(dFdt_comp)-length(Force_comp)==1
                            dFdt_comp = dFdt_comp(2:end);
                        end
                        
                        %%% MAIN SUBPLOT %%%
                        currentfig = figure(); % Call new figure
                        set(currentfig,'PaperOrientation','landscape')              % Sets figure to landscape
                        set(currentfig,'PaperUnits','normalized')                   % Normalize paper units
                        set(currentfig,'PaperPosition',[0 0 1 1])                   % Maximize figure size
                        set(currentfig,'visible','off')
                        axes('fontsize',14)                                         % Set fontsize
                        hold on                                                     % Turn hold on
                        ymin = -100;                                                % Y axis: lower limit
                        ymax = 350;                                                 % Y axis: upper limit
                        xmin = time(1);                                             % X axis: lower limit
                        xmax = time(end);                                           % X axis: upper limit
                        subplot(18,1,1:8)                                           % Set up subplot; designate 1:4 for main plot
                        axis([xmin xmax ymin ymax])                                 % Set up axes
                        hold on                                                     % Seems redundant, but it is necessary
                        plot(time,Force_comp,'b--','linewidth',2)                   % Plot Force component
                        plot(time,dFdt_comp,'r--','linewidth',2)                    % Plot dF/dt component
                        plot(time,static_comp,'g--','linewidth',2)                  % Plot Static component
                        plot(time,opt_fit,...    % Plot optimal fit
                            'Color',[1 0 0.6],...                                   %       *
                            'linewidth',3)                                          %       *
                        plot(FR_times,FR_recorded,'k.')                             % Plot recorded FR
                        plot([FR_times,FR_times],[ymin,ymin+18],'k')                % Plot spiketrain
                        hold on                                                     % Again, redundant, but necessary
                        title([filename ': perturbation ' num2str(pert)],...        % Figure title
                            'fontsize',14,...                                       %       *
                            'Interpreter','none')                                   %       *
                        ylabel('Firing Rate (spikes/s)','fontsize',14)              % Y axis: label
                        set(gca,'xtick',[],'color','none')                          % Remove xticks
                        leg = legend('Force Component',...                          % Figure legend
                            'dF/dt Component',...                                   %       *
                            'Static Component',...                                  %       *
                            'Optimal Fit',...                                       %       *
                            'Recorded FR',...                                       %       *
                            'Recorded Spikes');                                     %       *
                        set(leg,'FontSize',12)
                        text(xmax/2,ymax-20,...                                     % Force gain text
                            ['Force Gain: ' num2str(opt_gains(pert,1),4)...         %       *
                            ' (spikes/N*s)'],...                                    %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        text(xmax/2,ymax-50,...                                     % dF/dt gain text
                            ['dF/dt Gain: ' num2str(opt_gains(pert,2),4)...         %       *
                            ' (spikes/N)'],...                                      %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        text(xmax/2,ymax-80,...                                     % Static gain text
                            ['Static Gain: ' num2str(opt_gains(pert,3),4)...        %       *
                            ' (spikes/s)'],...                                      %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        R_squared = fit_data(model).R_squared(pert);    % R-squared value for trial
                        text(xmax/2,ymax-110,...                                    % R-squared text
                            ['R^2: ' num2str(R_squared,3)],...                      %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center',...                      %       *
                            'Interpreter','tex')                                    %       *
                        
                        %%% LENGTH SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(Length) - min(Length);                 % Y axis: range of values
                        ymin = min(Length) - 0.1*yamp;                    % Y axis: lower limit
                        ymax = max(Length) + 0.1*yamp;                    % Y axis: upper limit
                        subplot(18,1,9:10)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,Length)
                        ylabel('l (mm)','Fontsize',14)
                        set(gca,'xtick',[],'color','none')
                        
                        %%% VELOCITY SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(Velocity) - min(Velocity);             % Y axis: range of values
                        ymin = min(Velocity) - 0.1*yamp;                  % Y axis: lower limit
                        ymax = max(Velocity) + 0.1*yamp;                  % Y axis: upper limit
                        subplot(18,1,11:12)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,Velocity)
                        ylabel('v (mm/s)','Fontsize',14)
                        set(gca,'xtick',[],'color','none')
                        
                        %%% ACCEL SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(Accel) - min(Accel);                   % Y axis: range of values
                        ymin = min(Accel) - 0.1*yamp;                     % Y axis: lower limit
                        ymax = max(Accel) + 0.1*yamp;                     % Y axis: upper limit
                        subplot(18,1,13:14)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,Accel)
                        ylabel('a (mm/s^2)','Fontsize',14)
                        set(gca,'xtick',[],'color','none')
                        
                        %%% FORCE SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(Force) - min(Force);                   % Y axis: range of values
                        ymin = min(Force) - 0.1*yamp;                     % Y axis: lower limit
                        ymax = max(Force) + 0.1*yamp;                     % Y axis: upper limit
                        subplot(18,1,15:16)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,Force)
                        ylabel('F (N)','Fontsize',14)
                        set(gca,'xtick',[],'color','none')
                        
                        %%% dF/dt SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(dFdt) - min(dFdt);                     % Y axis: range of values
                        ymin = min(dFdt) - 0.1*yamp;                      % Y axis: lower limit
                        ymax = max(dFdt) + 0.1*yamp;                      % Y axis: upper limit
                        subplot(18,1,17:18)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,dFdt)
                        ylabel('dF/dt (N/s)','Fontsize',14)
                        xlabel('time(s)','fontsize',14)
                        set(gca,'color','none')
                        
                        
                        %%% Save Figure %%%
                        savepath = [pathname filesep...                         % Path in which to save file
                            'Plots' filesep...
                            'Expanded' filesep...
                            'model_' num2str(model) filesep];                 %         *
                        if ~exist(savepath,'dir')                               % If savepath doesn't exist...
                            mkdir(savepath)                                     % make it.
                        end                                                     %         *
                        savename = [savepath filesep filename(1:end-4)...       % Name with which to save the file
                            '_m' num2str(model)];                      %         *
                        print(currentfig,savename,'-dpsc','-append')            % Append figure to PostScript color file
                        close                                                   % Close the figure
%                         if pert == length(proc_data)                            % If this is the last pert. in trial...
%     %                             ps2pdf('psfile',[savename '.ps'],...                 % Convert .ps to .pdf (3rd party function)
% %                                 'pdffile',[savename '.pdf'],...               %         *
% %                                 'deletepsfile',1)                             %         *
%                         end
                        
                        %%%%% KINEMATICS MODEL START %%%%%
                    elseif model == 2;
                        
                        %%% Kinematics Data %%%
                        Length = proc_data(pert).Length;     % 2000 Hz Length
                        Velocity = proc_data(pert).Velocity; % 2000 Hz Velocity
                        Accel = proc_data(pert).dVdt;
                        
                        Force = proc_data(pert).Force;
                        dFdt = proc_data(pert).dFdt;                                        % dFdt
                        %%% Fit Data %%%
                        opt_gains = fit_data(model).opt_gains;                  % Optimal gains from fitting algorithm
                        Length_comp = opt_gains(pert,1)*Length;
                        Velocity_comp = opt_gains(pert,2)*Velocity;
                        Accel_comp = opt_gains(pert,3)*Accel;
                        Static_comp = opt_gains(pert,4)*ones(size(Length));
                        opt_fit = Length_comp + Velocity_comp + Accel_comp + Static_comp;
                        
                        %%% Timing %%%
                        time = proc_data(pert).time-proc_data(pert).time(1);          % Time relative to start of pert.
                        FR_times = fit_data(model).FR(pert).times;        % Times at which firing rate is defined
                        FR_recorded = fit_data(model).FR(pert).values;    % Values for FR from ISI
                        
                        %%% MAIN SUBPLOT %%%
                        currentfig = figure(); % Call new figure
                        set(currentfig,'PaperOrientation','landscape')    % Sets figure to landscape
                        set(currentfig,'PaperUnits','normalized')         % Normalize paper units
                        set(currentfig,'PaperPosition',[0 0 1 1])         % Maximize figure size
                        set(currentfig,'visible','off')
                        axes('fontsize',14)                               % Set fontsize
                        hold on                                           % Turn hold on
                        ymin = -100;                                      % Y axis: lower limit
                        ymax = 350;                                       % Y axis: upper limit
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        subplot(18,1,1:8)                                 % Set up subplot; designate 1:4 for main plot
                        axis([xmin xmax ymin ymax])                       % Set up axes
                        hold on                                           % Seems redundant, but it is necessary
                        plot(time,Length_comp,'b--','linewidth',2)        % Plot Length component
                        plot(time,Velocity_comp,'r--','linewidth',2)      % Plot Velocity component
                        plot(time,Accel_comp,'g--','linewidth',2)         % Plot Acceleration component
                        plot(time,Static_comp,'c--','linewidth',2)        % Plot Static component
                        plot(time,opt_fit,...                             % Plot optimal fit
                            'Color',[1 0 0.6],...                         %       *
                            'linewidth',3)                                %       *
                        plot(FR_times,FR_recorded,'k.')                   % Plot recorded FR
                        plot([FR_times,FR_times],[ymin,ymin+18],'k')      % Plot spiketrain
                        hold on                                           % Again, redundant, but necessary
                        title([filename ': perturbation ' num2str(pert)],...        % Figure title
                            'fontsize',14,...                                       %       *
                            'Interpreter','none')                                   %       *
                        ylabel('Firing Rate (spikes/s)','fontsize',14)              % Y axis: label
                        set(gca,'xtick',[],'color','none')                                         % Remove xticks
                        leg = legend('Length Component',...                         % Figure legend
                            'Velocity Component',...                                %       *
                            'Acceleration Component',...                            %       *
                            'Static Component',...
                            'Optimal Fit',...                                       %       *
                            'Recorded FR',...                                       %       *
                            'Recorded Spikes');                                     %       *
                        set(leg,'FontSize',12)
                        text(xmax/2,ymax-20,...                                     % Force gain text
                            ['Length Gain: ' num2str(opt_gains(pert,1),4)...        %       *
                            ' (spikes/mm*s)'],...                                    %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        text(xmax/2,ymax-50,...                                     % dF/dt gain text
                            ['Velocity Gain: ' num2str(opt_gains(pert,2),4)...      %       *
                            ' (spikes/mm)'],...                                      %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        text(xmax/2,ymax-80,...                                     % Static gain text
                            ['Acceleration Gain: ' num2str(opt_gains(pert,3),4)...  %       *
                            ' (spikes*s/mm)'],...                                      %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        text(xmax/2,ymax-110,...                                    % Static gain text
                            ['Static Gain: ' num2str(opt_gains(pert,4),4)...        %       *
                            ' (spikes/s)'],...                                      %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        R_squared = fit_data(model).R_squared(pert);    % R-squared value for trial
                        text(xmax/2,ymax-140,...                                    % R-squared text
                            ['R^2: ' num2str(R_squared,3)],...                      %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center',...                      %       *
                            'Interpreter','tex')                                    %       *
                        
                        %%% LENGTH SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(Length) - min(Length);                 % Y axis: range of values
                        ymin = min(Length) - 0.1*yamp;                    % Y axis: lower limit
                        ymax = max(Length) + 0.1*yamp;                    % Y axis: upper limit
                        subplot(18,1,9:10)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,Length)
                        ylabel('l (mm)','Fontsize',14)
                        set(gca,'xtick',[])
                        
                        %%% VELOCITY SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(Velocity) - min(Velocity);             % Y axis: range of values
                        ymin = min(Velocity) - 0.1*yamp;                  % Y axis: lower limit
                        ymax = max(Velocity) + 0.1*yamp;                  % Y axis: upper limit
                        subplot(18,1,11:12)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,Velocity)
                        ylabel('v (mm/s)','Fontsize',14)
                        set(gca,'xtick',[],'color','none')
                        
                        %%% ACCEL SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(Accel) - min(Accel);                   % Y axis: range of values
                        ymin = min(Accel) - 0.1*yamp;                     % Y axis: lower limit
                        ymax = max(Accel) + 0.1*yamp;                     % Y axis: upper limit
                        subplot(18,1,13:14)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,Accel)
                        ylabel('a (mm/s^2)','Fontsize',14)
                        set(gca,'xtick',[],'color','none')
                        
                        %%% FORCE SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(Force) - min(Force);                   % Y axis: range of values
                        ymin = min(Force) - 0.1*yamp;                     % Y axis: lower limit
                        ymax = max(Force) + 0.1*yamp;                     % Y axis: upper limit
                        subplot(18,1,15:16)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,Force)
                        ylabel('F (N)','Fontsize',14)
                        set(gca,'xtick',[],'color','none')
                        
                        %%% dF/dt SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(dFdt) - min(dFdt);                     % Y axis: range of values
                        ymin = min(dFdt) - 0.1*yamp;                      % Y axis: lower limit
                        ymax = max(dFdt) + 0.1*yamp;                      % Y axis: upper limit
                        subplot(18,1,17:18)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,dFdt)
                        ylabel('dF/dt (N/s)','Fontsize',14)
                        xlabel('time(s)','fontsize',14)
                        set(gca,'color','none')
                        
                        %%% Save Figure %%%
                        savepath = [pathname filesep...                         % Path in which to save file
                            'Plots' filesep...
                            'Expanded' filesep...
                            'model_' num2str(model) filesep];                 %         *
                        if ~exist(savepath,'dir')                               % If savepath doesn't exist...
                            mkdir(savepath)                                     % make it.
                        end
                        savename = [savepath filesep filename(1:end-4)...       % Name with which to save the file
                            '_m' num2str(model)];                      %         *
                        print(currentfig,savename,'-dpsc','-append')            % Append figure to PostScript color file
                        close                                                   % Close the figure
%                         if pert == length(proc_data) 
%                             % If this is the last pert. in trial...
% %                             ps2pdf('psfile',[savename '.ps'],...                 % Convert .ps to .pdf (3rd party function)
% %                                 'pdffile',[savename '.pdf'],...               %         *
% %                                 'deletepsfile',1)                             %         *
%                         end
                        %%%%% KINEMATICS MODEL START %%%%%
                    elseif model == 2;
                        
                        %%% Kinematics Data %%%
                        Length = proc_data(pert).Length;     % 2000 Hz Length
                        Velocity = proc_data(pert).Velocity; % 2000 Hz Velocity
                        Accel = proc_data(pert).dVdt;
                        
                        Force = proc_data(pert).Force;
                        dFdt = proc_data(pert).dFdt;                                        % dFdt
                        %%% Fit Data %%%
                        opt_gains = fit_data(model).opt_gains;                  % Optimal gains from fitting algorithm
                        Length_comp = opt_gains(pert,1)*Length;
                        Velocity_comp = opt_gains(pert,2)*Velocity;
                        Accel_comp = opt_gains(pert,3)*Accel;
                        Static_comp = opt_gains(pert,4)*ones(size(Length));
                        opt_fit = Length_comp + Velocity_comp + Accel_comp + Static_comp;
                        
                        %%% Timing %%%
                        time = proc_data(pert).time-proc_data(pert).time(1);          % Time relative to start of pert.
                        FR_times = fit_data(model).FR(pert).times;        % Times at which firing rate is defined
                        FR_recorded = fit_data(model).FR(pert).values;    % Values for FR from ISI
                        
                        %%% MAIN SUBPLOT %%%
                        currentfig = figure(); % Call new figure
                        set(currentfig,'PaperOrientation','landscape')    % Sets figure to landscape
                        set(currentfig,'PaperUnits','normalized')         % Normalize paper units
                        set(currentfig,'PaperPosition',[0 0 1 1])         % Maximize figure size
                        set(currentfig,'visible','off')
                        axes('fontsize',14)                               % Set fontsize
                        hold on                                           % Turn hold on
                        ymin = -100;                                      % Y axis: lower limit
                        ymax = 350;                                       % Y axis: upper limit
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        subplot(18,1,1:8)                                 % Set up subplot; designate 1:4 for main plot
                        axis([xmin xmax ymin ymax])                       % Set up axes
                        hold on                                           % Seems redundant, but it is necessary
                        plot(time,Length_comp,'b--','linewidth',2)        % Plot Length component
                        plot(time,Velocity_comp,'r--','linewidth',2)      % Plot Velocity component
                        plot(time,Accel_comp,'g--','linewidth',2)         % Plot Acceleration component
                        plot(time,Static_comp,'c--','linewidth',2)        % Plot Static component
                        plot(time,opt_fit,...                             % Plot optimal fit
                            'Color',[1 0 0.6],...                         %       *
                            'linewidth',3)                                %       *
                        plot(FR_times,FR_recorded,'k.')                   % Plot recorded FR
                        plot([FR_times,FR_times],[ymin,ymin+18],'k')      % Plot spiketrain
                        hold on                                           % Again, redundant, but necessary
                        title([filename ': perturbation ' num2str(pert)],...        % Figure title
                            'fontsize',14,...                                       %       *
                            'Interpreter','none')                                   %       *
                        ylabel('Firing Rate (spikes/s)','fontsize',14)              % Y axis: label
                        set(gca,'xtick',[],'color','none')                                         % Remove xticks
                        leg = legend('Length Component',...                         % Figure legend
                            'Velocity Component',...                                %       *
                            'Acceleration Component',...                            %       *
                            'Static Component',...
                            'Optimal Fit',...                                       %       *
                            'Recorded FR',...                                       %       *
                            'Recorded Spikes');                                     %       *
                        set(leg,'FontSize',12)
                        text(xmax/2,ymax-20,...                                     % Force gain text
                            ['Length Gain: ' num2str(opt_gains(pert,1),4)...        %       *
                            ' (spikes/mm*s)'],...                                    %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        text(xmax/2,ymax-50,...                                     % dF/dt gain text
                            ['Velocity Gain: ' num2str(opt_gains(pert,2),4)...      %       *
                            ' (spikes/mm)'],...                                      %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        text(xmax/2,ymax-80,...                                     % Static gain text
                            ['Acceleration Gain: ' num2str(opt_gains(pert,3),4)...  %       *
                            ' (spikes*s/mm)'],...                                      %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        text(xmax/2,ymax-110,...                                    % Static gain text
                            ['Static Gain: ' num2str(opt_gains(pert,4),4)...        %       *
                            ' (spikes/s)'],...                                      %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        R_squared = fit_data(model).R_squared(pert);    % R-squared value for trial
                        text(xmax/2,ymax-140,...                                    % R-squared text
                            ['R^2: ' num2str(R_squared,3)],...                      %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center',...                      %       *
                            'Interpreter','tex')                                    %       *
                        
                        %%% LENGTH SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(Length) - min(Length);                 % Y axis: range of values
                        ymin = min(Length) - 0.1*yamp;                    % Y axis: lower limit
                        ymax = max(Length) + 0.1*yamp;                    % Y axis: upper limit
                        subplot(18,1,9:10)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,Length)
                        ylabel('l (mm)','Fontsize',14)
                        set(gca,'xtick',[])
                        
                        %%% VELOCITY SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(Velocity) - min(Velocity);             % Y axis: range of values
                        ymin = min(Velocity) - 0.1*yamp;                  % Y axis: lower limit
                        ymax = max(Velocity) + 0.1*yamp;                  % Y axis: upper limit
                        subplot(18,1,11:12)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,Velocity)
                        ylabel('v (mm/s)','Fontsize',14)
                        set(gca,'xtick',[],'color','none')
                        
                        %%% ACCEL SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(Accel) - min(Accel);                   % Y axis: range of values
                        ymin = min(Accel) - 0.1*yamp;                     % Y axis: lower limit
                        ymax = max(Accel) + 0.1*yamp;                     % Y axis: upper limit
                        subplot(18,1,13:14)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,Accel)
                        ylabel('a (mm/s^2)','Fontsize',14)
                        set(gca,'xtick',[],'color','none')
                        
                        %%% FORCE SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(Force) - min(Force);                   % Y axis: range of values
                        ymin = min(Force) - 0.1*yamp;                     % Y axis: lower limit
                        ymax = max(Force) + 0.1*yamp;                     % Y axis: upper limit
                        subplot(18,1,15:16)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,Force)
                        ylabel('F (N)','Fontsize',14)
                        set(gca,'xtick',[],'color','none')
                        
                        %%% dF/dt SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(dFdt) - min(dFdt);                     % Y axis: range of values
                        ymin = min(dFdt) - 0.1*yamp;                      % Y axis: lower limit
                        ymax = max(dFdt) + 0.1*yamp;                      % Y axis: upper limit
                        subplot(18,1,17:18)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,dFdt)
                        ylabel('dF/dt (N/s)','Fontsize',14)
                        xlabel('time(s)','fontsize',14)
                        set(gca,'color','none')
                        
                        %%% Save Figure %%%
                        savepath = [pathname filesep...                         % Path in which to save file
                            'Plots' filesep...
                            'Expanded' filesep...
                            'model_' num2str(model) filesep];                 %         *
                        if ~exist(savepath,'dir')                               % If savepath doesn't exist...
                            mkdir(savepath)                                     % make it.
                        end
                        savename = [savepath filesep filename(1:end-4)...       % Name with which to save the file
                            '_m' num2str(model)];                      %         *
                        print(currentfig,savename,'-dpsc','-append')            % Append figure to PostScript color file
                        close                                                   % Close the figure
%                         if pert == length(proc_data) 
%                             % If this is the last pert. in trial...
% %                             ps2pdf('psfile',[savename '.ps'],...                 % Convert .ps to .pdf (3rd party function)
% %                                 'pdffile',[savename '.pdf'],...               %         *
% %                                 'deletepsfile',1)                             %         *
%                         end

                        %%%%% MIXED MODEL START %%%%%
                    elseif model == 3;
                        
                        %%% Kinematics Data %%%
                        Length = proc_data(pert).Length;     % 2000 Hz Length
                        Velocity = proc_data(pert).Velocity; % 2000 Hz Velocity
                        Accel = proc_data(pert).dVdt;
                        
                        Force = proc_data(pert).Force;
                        dFdt = proc_data(pert).dFdt;                                        % dFdt
                        %%% Fit Data %%%
                        opt_gains = fit_data(model).opt_gains;                  % Optimal gains from fitting algorithm
                        Force_comp = opt_gains(pert,1)*Force;
                        dFdt_comp = opt_gains(pert,2)*dFdt;
                        Length_comp = opt_gains(pert,3)*Length;
                        Velocity_comp = opt_gains(pert,4)*Velocity;
                        Accel_comp = opt_gains(pert,5)*Accel;
                        Static_comp = opt_gains(pert,6)*ones(size(Length));
                        opt_fit = Force_comp + dFdt_comp + Length_comp + Velocity_comp + Accel_comp + Static_comp;
                        
                        %%% Timing %%%
                        time = proc_data(pert).time-proc_data(pert).time(1);          % Time relative to start of pert.
                        FR_times = fit_data(model).FR(pert).times;        % Times at which firing rate is defined
                        FR_recorded = fit_data(model).FR(pert).values;    % Values for FR from ISI
                        
                        %%% MAIN SUBPLOT %%%
                        currentfig = figure(); % Call new figure
                        set(currentfig,'PaperOrientation','landscape')    % Sets figure to landscape
                        set(currentfig,'PaperUnits','normalized')         % Normalize paper units
                        set(currentfig,'PaperPosition',[0 0 1 1])         % Maximize figure size
                        set(currentfig,'visible','off')
                        axes('fontsize',14)                               % Set fontsize
                        hold on                                           % Turn hold on
                        ymin = -100;                                      % Y axis: lower limit
                        ymax = 350;                                       % Y axis: upper limit
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        subplot(18,1,1:8)                                 % Set up subplot; designate 1:4 for main plot
                        axis([xmin xmax ymin ymax])                       % Set up axes
                        hold on                                           % Seems redundant, but it is necessary
                        plot(time,Length_comp,'b--','linewidth',2)        % Plot Length component
                        plot(time,Velocity_comp,'r--','linewidth',2)      % Plot Velocity component
                        plot(time,Accel_comp,'g--','linewidth',2)         % Plot Acceleration component
                        plot(time,Static_comp,'c--','linewidth',2)        % Plot Static component
                        plot(time,Force_comp,'--','linewidth',2)
                        plot(time,dFdt_comp,'--','linewidth',2)
                        plot(time,opt_fit,...                             % Plot optimal fit
                            'Color',[1 0 0.6],...                         %       *
                            'linewidth',3)                                %       *
                        plot(FR_times,FR_recorded,'k.')                   % Plot recorded FR
                        plot([FR_times,FR_times],[ymin,ymin+18],'k')      % Plot spiketrain
                        hold on                                           % Again, redundant, but necessary
                        title([filename ': perturbation ' num2str(pert)],...        % Figure title
                            'fontsize',14,...                                       %       *
                            'Interpreter','none')                                   %       *
                        ylabel('Firing Rate (spikes/s)','fontsize',14)              % Y axis: label
                        set(gca,'xtick',[],'color','none')                                         % Remove xticks
                        leg = legend('Length Component',...                         % Figure legend
                            'Velocity Component',...                                %       *
                            'Acceleration Component',...                            %       *
                            'Static Component',...
                            'Force Component',...
                            'dFdt Component',...
                            'Optimal Fit',...                                       %       *
                            'Recorded FR',...                                       %       *
                            'Recorded Spikes');                                     %       *
                        set(leg,'FontSize',12)
                        text(xmax/2,ymax-20,...                                     % Force gain text
                            ['Force Gain: ' num2str(opt_gains(pert,1),4)...        %       *
                            ' (spikes/g*s)'],...                                    %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        text(xmax/2,ymax-50,...                                     % Force gain text
                            ['dF/dt Gain: ' num2str(opt_gains(pert,3),4)...        %       *
                            ' (spikes/g)'],...                                    %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        text(xmax/2,ymax-80,...                                     % Force gain text
                            ['Length Gain: ' num2str(opt_gains(pert,3),4)...        %       *
                            ' (spikes/mm*s)'],...                                    %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        text(xmax/2,ymax-110,...                                     % dF/dt gain text
                            ['Velocity Gain: ' num2str(opt_gains(pert,4),4)...      %       *
                            ' (spikes/mm)'],...                                      %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        text(xmax/2,ymax-140,...                                     % Static gain text
                            ['Acceleration Gain: ' num2str(opt_gains(pert,5),4)...  %       *
                            ' (spikes*s/mm)'],...                                      %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        text(xmax/2,ymax-170,...                                    % Static gain text
                            ['Static Gain: ' num2str(opt_gains(pert,6),4)...        %       *
                            ' (spikes/s)'],...                                      %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center')                         %       *
                        R_squared = fit_data(model).R_squared(pert);    % R-squared value for trial
                        text(xmax/2,ymax-200,...                                    % R-squared text
                            ['R^2: ' num2str(R_squared,3)],...                      %       *
                            'fontsize',14,...                                       %       *
                            'HorizontalAlignment','center',...                      %       *
                            'Interpreter','tex')                                    %       *
                        
                        %%% LENGTH SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(Length) - min(Length);                 % Y axis: range of values
                        ymin = min(Length) - 0.1*yamp;                    % Y axis: lower limit
                        ymax = max(Length) + 0.1*yamp;                    % Y axis: upper limit
                        subplot(18,1,9:10)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,Length)
                        ylabel('l (mm)','Fontsize',14)
                        set(gca,'xtick',[])
                        
                        %%% VELOCITY SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(Velocity) - min(Velocity);             % Y axis: range of values
                        ymin = min(Velocity) - 0.1*yamp;                  % Y axis: lower limit
                        ymax = max(Velocity) + 0.1*yamp;                  % Y axis: upper limit
                        subplot(18,1,11:12)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,Velocity)
                        ylabel('v (mm/s)','Fontsize',14)
                        set(gca,'xtick',[],'color','none')
                        
                        %%% ACCEL SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(Accel) - min(Accel);                   % Y axis: range of values
                        ymin = min(Accel) - 0.1*yamp;                     % Y axis: lower limit
                        ymax = max(Accel) + 0.1*yamp;                     % Y axis: upper limit
                        subplot(18,1,13:14)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,Accel)
                        ylabel('a (mm/s^2)','Fontsize',14)
                        set(gca,'xtick',[],'color','none')
                        
                        %%% FORCE SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(Force) - min(Force);                   % Y axis: range of values
                        ymin = min(Force) - 0.1*yamp;                     % Y axis: lower limit
                        ymax = max(Force) + 0.1*yamp;                     % Y axis: upper limit
                        subplot(18,1,15:16)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,Force)
                        ylabel('F (N)','Fontsize',14)
                        set(gca,'xtick',[],'color','none')
                        
                        %%% dF/dt SUBPLOT %%%
                        xmin = time(1);                                   % X axis: lower limit
                        xmax = time(end);                                 % X axis: upper limit
                        yamp = max(dFdt) - min(dFdt);                     % Y axis: range of values
                        ymin = min(dFdt) - 0.1*yamp;                      % Y axis: lower limit
                        ymax = max(dFdt) + 0.1*yamp;                      % Y axis: upper limit
                        subplot(18,1,17:18)
                        axis([xmin xmax ymin ymax])
                        hold on
                        plot(time,dFdt)
                        ylabel('dF/dt (N/s)','Fontsize',14)
                        xlabel('time(s)','fontsize',14)
                        set(gca,'color','none')
                        
                        %%% Save Figure %%%
                        savepath = [pathname filesep...                         % Path in which to save file
                            'Plots' filesep...
                            'Expanded' filesep...
                            'model_' num2str(model) filesep];                 %         *
                        if ~exist(savepath,'dir')                               % If savepath doesn't exist...
                            mkdir(savepath)                                     % make it.
                        end
                        savename = [savepath filesep filename(1:end-4)...       % Name with which to save the file
                            '_m' num2str(model)];                      %         *
                        print(currentfig,savename,'-dpsc','-append')            % Append figure to PostScript color file
                        close                                                   % Close the figure

                    end % First Model end
                end % Perturbation end
                
            end % Models end
        end % ".mat file?" end
    end % Directory end
end
toc