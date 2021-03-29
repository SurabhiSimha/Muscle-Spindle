%% Single set of Triangles & RH
    clear,clc
    tic
    time_step = 0.001; %Temporal precision
    t = -5:time_step:25; % Time vector
    pertStart = 5000;
    RHStart = 18000;
    numSims = 1;       % Number of simulations to run in parallel
    delta_cdl = zeros(numSims,numel(t)); % change in command length for all sims
%     delta_Ca = zeros(numSims,numel(t)); % change in [Ca] for all sims
    delta_f_activated = zeros(numSims,numel(t));
    strDur = 600; % duration of stretch period
    lsf = 0.8; % length scaling factor to account for pinnation & elastic attachment of fibers
    
    for a = 1:numSims
        for i = 1:numel(t)
            
            if i == 1
                delta_f_activated(a,i) = 0.3;
            elseif i > pertStart && i < pertStart + strDur
                delta_cdl(a,i) = 0.1182*lsf;
            elseif i > pertStart + strDur && i < pertStart + 2*(strDur)
                delta_cdl(a,i) = -0.1182*lsf;
            elseif i > pertStart + 2*(strDur) && i < pertStart + 3*(strDur)
                delta_cdl(a,i) = 0.1182*lsf;
            elseif i > pertStart + 3*(strDur) && i < pertStart + 4*(strDur)
                delta_cdl(a,i) = -0.1182*lsf;
            elseif i > pertStart + 4*(strDur) + 1e3*(a-1) && i < pertStart + 5*(strDur) + 1e3*(a-1)
                delta_cdl(a,i) = 0.1182*lsf;
            elseif i > pertStart + 5*(strDur) + 1e3*(a-1) && i < pertStart + 6*(strDur) + 1e3*(a-1)
                delta_cdl(a,i) = -0.1182*lsf;
            elseif i > RHStart && i < RHStart + strDur/a
                delta_cdl(a,i) = 0.1182*a*lsf;
            elseif i > RHStart + strDur/a + 1000 && i < RHStart + 2*(strDur/a) + 1000 
                delta_cdl(a,i) = -0.1182*a*lsf;
            end
        end
    end
    
    parfor a = 1:numSims
        [hsB(a),dataB(a),hsC(a),dataC(a)] = sarcSimDriver(t,delta_f_activated(a,:),delta_cdl(a,:));
        disp(['Done with simulation number ' num2str(a)])
    end
    
    beep; toc;


    
    %% DI
    clear,clc
    tic
    time_step = 0.001; %Temporal precision
    t = -7:time_step:6; % Time vector
    pertStart = 7000;
    numSims = 10;       % Number of simulations to run in parallel
    delta_cdl = zeros(numSims,numel(t)); % change in command length for all sims
%     delta_Ca = zeros(numSims,numel(t)); % change in [Ca] for all sims
    delta_f_activated = zeros(numSims,numel(t));
    strDur = 600; % duration of stretch period
    lsf = 0.8; % length scaling factor to account for pinnation & elastic attachment of fibers
    
    for a = 1:numSims
        for i = 1:numel(t)
            if i == 1
                delta_f_activated(a,i) = 0.3;
            elseif i > pertStart && i < pertStart + strDur/a
                delta_cdl(a,i) = 0.1182*a*lsf;
            elseif i > pertStart + strDur/a + 1000 && i < pertStart + 2*(strDur/a) + 1000 
                delta_cdl(a,i) = -0.1182*a*lsf;
            end
        end
    end
    
    parfor a = 1:numSims
        [hsB(a),dataB(a),hsC(a),dataC(a)] = sarcSimDriver(t,delta_f_activated(a,:),delta_cdl(a,:));
        disp(['Done with simulation number ' num2str(a)])
    end
    
    beep; toc;
    

%% THD
    clear,clc
    tic
    time_step = 0.001; %Temporal precision
    t = -5:time_step:20; % Time vector
    pertStart = 5e3;
    numSims = 16;       % Number of simulations to run in parallel
    delta_cdl = zeros(numSims,numel(t)); % change in command length for all sims
%     delta_Ca = zeros(numSims,numel(t)); % change in [Ca] for all sims
    delta_f_activated = zeros(numSims,numel(t));
    strDur = 600; % duration of stretch period
    lsf = 0.8; % length scaling factor to account for pinnation & elastic attachment of fibers
    
    for a = 1:numSims
        for i = 1:numel(t)
            if i == 1
                delta_f_activated(a,i) = 0.3;
            elseif i > pertStart && i < pertStart + strDur
                delta_cdl(a,i) = 0.1182*lsf;
            elseif i > pertStart + strDur && i < pertStart + 2*(strDur)
                delta_cdl(a,i) = -0.1182*lsf;
            elseif i > pertStart + 2*(strDur) && i < pertStart + 3*(strDur)
                delta_cdl(a,i) = 0.1182*lsf;
            elseif i > pertStart + 3*(strDur) && i < pertStart + 4*(strDur)
                delta_cdl(a,i) = -0.1182*lsf;
            elseif i > pertStart + 4*(strDur) + 5e2*(a-1) && i < pertStart + 5*(strDur) + 5e2*(a-1)
                delta_cdl(a,i) = 0.1182*lsf;
            elseif i > pertStart + 5*(strDur) + 5e2*(a-1) && i < pertStart + 6*(strDur) + 5e2*(a-1)
                delta_cdl(a,i) = -0.1182*lsf;
            end
        end
    end
    
    parfor a = 1:numSims
        [hsB(a),dataB(a),hsC(a),dataC(a)] = sarcSimDriver(t,delta_f_activated(a,:),delta_cdl(a,:));
        
        disp(['Done with simulation number ' num2str(a)])
    end
    
    beep; toc;

%% AHD
    clear,clc
    tic
    time_step = 0.001; %Temporal precision
    t = -5:time_step:5; % Time vector
    pertStart = 5e3;
    numSims = 20;       % Number of simulations to run in parallel
    delta_cdl = zeros(numSims,numel(t)); % change in command length for all sims
%     delta_Ca = zeros(numSims,numel(t)); % change in [Ca] for all sims
    delta_f_activated = zeros(numSims,numel(t));
    strDur = 600; % duration of stretch period
    lsf = 0.8; % length scaling factor to account for pinnation & elastic attachment of fibers
    
    amps = 0.0:0.05:numSims*0.05 ;% 1 = 5.46% stretch (0.1182/ms for 600ms)
    
    for a = 1:numSims
        for i = 1:numel(t)
            if i == 1
                delta_f_activated(a,i) = 0.3;
            elseif i > pertStart - 2*amps(a)*strDur && i < pertStart - amps(a)*strDur
                delta_cdl(a,i) = 0.1182*lsf;
            elseif i > pertStart - amps(a)*strDur && i < pertStart
                delta_cdl(a,i) = -0.1182*lsf;
            elseif i > pertStart && i < pertStart + strDur
                delta_cdl(a,i) = 0.1182*lsf;
            elseif i > pertStart + strDur && i < pertStart + 2*(strDur)
                delta_cdl(a,i) = -0.1182*lsf;
            elseif i > pertStart + 2*(strDur) && i < pertStart + 3*(strDur)
                delta_cdl(a,i) = 0.1182*lsf;
            elseif i > pertStart + 3*(strDur) && i < pertStart + 4*(strDur)
                delta_cdl(a,i) = -0.1182*lsf;
            end
        end
    end
    
    parfor a = 1:numSims
        [hsB(a),dataB(a),hsC(a),dataC(a)] = sarcSimDriver(t,delta_f_activated(a,:),delta_cdl(a,:));
        disp(['Done with simulation number ' num2str(a)])
    end
    
    beep; toc;


%% Act RH
    clear,clc
    tic
    time_step = 0.001; %Temporal precision
    t = -5:time_step:5; % Time vector
    pertStart = 5e3;
    numSims = 8;       % Number of simulations to run in parallel
    delta_cdl = zeros(numSims,numel(t)); % change in command length for all sims
%     delta_Ca = zeros(numSims,numel(t)); % change in [Ca] for all sims
    delta_f_activated = zeros(numSims,numel(t));
    strDur = 600; % duration of stretch period
    lsf = 0.8; % length scaling factor to account for pinnation & elastic attachment of fibers
    
%     Ca_levels = [10^(-7) 10^(-6.5) 10^(-6.0) 10^(-4.5) 10^(-4.5)] - 10^(-6.5);
    
    for a = 1:numSims
        for i = 1:numel(t)
            if i == 1
%                 delta_Ca(a,i) = Ca_levels(a);
                delta_f_activated(a,i) = 0.1*(a)+0.1;
            elseif i > pertStart && i < pertStart + strDur
                delta_cdl(a,i) = 0.1182*lsf;
%                 delta_Ca(a,i) = 0;
            end
        end
    end
    
    parfor a = 1:numSims
        [hsB(a),dataB(a),hsC(a),dataC(a)] = sarcSimDriver(t,delta_f_activated(a,:),delta_cdl(a,:));
        disp(['Done with simulation number ' num2str(a)])
    end
    
    beep; toc;
 %% sinusoid (Day 2017)
    clear,clc
    tic
    time_step = 0.001; %Temporal precision
    t = -7:time_step:25; % Time vector
    pertStart = 7000;
    numSims = 1;       % Number of simulations to run in parallel
    delta_cdl = zeros(numSims,numel(t)); % change in command length for all sims
    delta_f_activated = zeros(numSims,numel(t));
    strDur = 2e4; % duration of stretch period
    lsf = cos(pi/6); % length scaling factor to account for pinnation & elastic attachment of fibers
    freq = 1.57; %rad/s taken from Day et al. 2017
    
    for a = 1:numSims
        for i = 1:numel(t)
            if i == 1
                delta_f_activated(a,i) = 0.3;
            elseif i > pertStart && i < pertStart + strDur
                delta_cdl(a,i) = 0.05*sin(freq*(t(i)-t(pertStart))+pi/4);
            
            end
        end
    end
    
    parfor a = 1:numSims
        [hsB(a),dataB(a),hsC(a),dataC(a)] = sarcSimDriver(t,delta_f_activated(a,:),delta_cdl(a,:));
        disp(['Done with simulation number ' num2str(a)])
    end
    
    beep; toc;
    %% History dependence (Haftel 2004)
    
    clear,clc
    tic
    time_step = 0.001; %Temporal precision
    t = -5:time_step:10; % Time vector
    pertStart = 5e3;
    numSims = 1;       % Number of simulations to run in parallel
    delta_cdl = zeros(numSims,numel(t)); % change in command length for all sims
%     delta_Ca = zeros(numSims,numel(t)); % change in [Ca] for all sims
    delta_f_activated = zeros(numSims,numel(t));
    strDur = 750; % duration of stretch period
    lsf = 1; % length scaling factor to account for pinnation & elastic attachment of fibers
    
    for a = 1:numSims
        for i = 1:numel(t)
            if i == 1
                delta_f_activated(a,i) = 0.3;
            elseif i > pertStart && i < pertStart + strDur
                delta_cdl(a,i) = 0.1182*lsf; %0.1182 nm/ms is 0.0909L0/s (3mm/s for a 44mm muscle)
            elseif i > pertStart + strDur && i < pertStart + 2*(strDur)
                delta_cdl(a,i) = -0.1182*lsf;
            elseif i > pertStart + 2*(strDur) && i < pertStart + 3*(strDur)
                delta_cdl(a,i) = 0.1182*lsf;
            elseif i > pertStart + 3*(strDur) && i < pertStart + 4*(strDur)
                delta_cdl(a,i) = -0.1182*lsf;
            elseif i > pertStart + 4*(strDur) + 1e3*(a-1) && i < pertStart + 5*(strDur) + 1e3*(a-1)
                delta_cdl(a,i) = 0.1182*lsf;
            elseif i > pertStart + 5*(strDur) + 1e3*(a-1) && i < pertStart + 6*(strDur) + 1e3*(a-1)
                delta_cdl(a,i) = -0.1182*lsf;
            end
        end
    end
    
    parfor a = 1:numSims
        [hsB(a),dataB(a),hsC(a),dataC(a)] = sarcSimDriver(t,delta_f_activated(a,:),delta_cdl(a,:));
        disp(['Done with simulation number ' num2str(a)])
    end
    
    beep; toc;


    %% Alpha-gamma isometric ramp
    
    clear,clc
    tic
    time_step = 0.001; %Temporal precision
    t = -5:time_step:6; % Time vector
    pertStart = 5e3;
    numSims = 8;       % Number of simulations to run in parallel
    delta_cdl = zeros(numSims,numel(t)); % change in command length for all sims
%     delta_Ca = zeros(numSims,numel(t)); % change in [Ca] for all sims
    delta_f_activated = zeros(numSims,numel(t));
    strDur = 2000; % duration of activation period
    lsf = 0.8; % length scaling factor to account for pinnation & elastic attachment of fibers
    
    for a = 1:numSims
        for i = 1:numel(t)
            if i == 1
                delta_f_activated(a,i) = 0.1;
            elseif i > pertStart && i < pertStart + strDur
                delta_cdl(a,i) = -0.0015*a*lsf + 0.0015*lsf;
                delta_f_activated(a,i) = 0.0004;
%             elseif i > pertStart+strDur && i < pertStart + strDur + 10
%                 if a > 1
%                     delta_cdl(a,i) = 0.1*lsf; %elastic rebound of tendon (none if tendon is rigid)
%                 end
            elseif i == pertStart+strDur %Difference between 10 m/s and 100 m/s axon delays
                delta_f_activated(a,i) = -0.0004*strDur;                
            end
        end
    end
    
    parfor a = 1:numSims
        [hsB(a),dataB(a),hsC(a),dataC(a)] = sarcSimDriver(t,delta_f_activated(a,:),delta_cdl(a,:));
        disp(['Done with simulation number ' num2str(a)])
    end
    
    beep; toc;

    