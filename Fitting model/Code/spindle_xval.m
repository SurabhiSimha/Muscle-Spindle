%% Cross Validation of Kinetics and Kinematics Spindle Firing Models (KPB 04/10/2015)
% This code tests the ability of each model to reproduce (training) and
% predict (testing) randomly selected datapoints of the complete dataset
% for each afferent.
%%% External functions used
% <spindle_MLEcost.html spindle_MLEcost.m spindle_getMLE.m>

%% Start Clock, Clear Variables, and Specify Afferents for Analysis
tic
clear
affnums = [24,25,26,27,32,56,57,58,60,64];

%% Preallocate Variables for Each Iteration (Afferent)
for aff = 1:length(affnums)
    affname = ['aff' num2str(affnums(aff))];
    pathname = cd(cd(['..' filesep 'Data' filesep 'proc_data_no_buffer' filesep affname '_proc']));
    directory = dir(pathname);
    xvaldata = struct('Length',[],'Velocity',[],'Accel',[],'Force',[],'dFdt',[],...
        'LengthF','VelocityF','dVFdt','LengthF_hc','VelocityF_hc','dVFdt_hc','IFR',[]);
    flags.aff = affnums(aff);
    concatLength = [];
    concatVel = [];
    concatAccel = [];
    concatForce = [];
    concatdFdt = [];
    concatLengthF = [];
    concatVelocityF = [];
    concatdVFdt = [];
    concatLengthF_hc = [];
    concatVelocityF_hc = [];
    concatdVFdt_hc = [];
    concatIFR = [];
    concatPertEnd = [];
    
    %% Go to Each File to Check if it is .mat
    for trial = 1:length(directory) % For all files in directory
        filename = directory(trial).name; % Current file
        if strfind(filename,'.mat')    % If current file is .mat file 
            loadname = [pathname filesep filename];
            load(loadname); % Load file
    %% Concatenate All Trials for an Afferent Into a Single Array
        % This is data used for fitting in proc_data_no_buffer (lagged
        % values, using mean lags from individual fits). 
            for pert = 1:size(fit_data(1).opt_gains,1)
                Length = fit_data(2).fitting_data(pert).Length;
                Velocity = fit_data(2).fitting_data(pert).Velocity;
                Accel = fit_data(2).fitting_data(pert).dVdt;
                Force = fit_data(1).fitting_data(pert).Force;
                dFdt = fit_data(1).fitting_data(pert).dFdt;
                LenF = fit_data(4).fitting_data(pert).LengthF;
                VelF = fit_data(4).fitting_data(pert).VelocityF;
                dVFdt = fit_data(4).fitting_data(pert).dVFdt;
                LenFhc = fit_data(5).fitting_data(pert).LengthF_hc;
                VelFhc = fit_data(5).fitting_data(pert).VelocityF_hc;
                dVFdthc = fit_data(5).fitting_data(pert).dVFdt_hc;
                IFR = fit_data(1).FR(pert).values;
                n = numel(IFR);
                PertEnd = [zeros(n-1,1); 1];

                concatPertEnd = vertcat(concatPertEnd,PertEnd);
                concatLength = vertcat(concatLength,Length);
                concatVel = vertcat(concatVel,Velocity);
                concatAccel = vertcat(concatAccel,Accel);
                concatForce = vertcat(concatForce,Force);
                concatdFdt = vertcat(concatdFdt,dFdt);
                concatLengthF = vertcat(concatLengthF,LenF);
                concatVelocityF = vertcat(concatVelocityF,VelF);
                concatdVFdt = vertcat(concatdVFdt,dVFdt);
                concatLengthF_hc = vertcat(concatLengthF_hc,LenFhc);
                concatVelocityF_hc = vertcat(concatVelocityF_hc,VelFhc);
                concatdVFdt_hc = vertcat(concatdVFdt_hc,dVFdthc);
                concatIFR = vertcat(concatIFR,IFR);
                
            end
        end
    end
%% Set up Initial Parameters for Optimization
% These variables will be fed into _fmincon.m_ and used as boundaries and
% conditions for curve fitting and information criteria fitting


                                
    numParam = [6,8,6,9,9,28]; % Number of free parameters in each model (gains + offsets + 2 for variance & lag, +1 for tendon)
%     numObs = length(concatIFR); % Number of total observations (action potentials)
    numObs = sum(concatPertEnd); % Number of total observations (each perturbation; action potentials aren't independent)
    PertStart = [1; find(concatPertEnd(1:end-1) == 1) + 1];
    PertEnd = find(concatPertEnd == 1);
%% Populate Data Structure for Optimization   
    xvaldata.Length = concatLength;
    xvaldata.Velocity = concatVel;
    xvaldata.Accel = concatAccel;
    xvaldata.Force = concatForce;
    xvaldata.dFdt = concatdFdt;
    xvaldata.LengthF = concatLengthF;
    xvaldata.VelocityF = concatVelocityF;
    xvaldata.dVFdt = concatdVFdt;
    xvaldata.LengthF_hc = concatLengthF_hc;
    xvaldata.VelocityF_hc = concatVelocityF_hc;
    xvaldata.dVFdt_hc = concatdVFdt_hc;
    xvaldata.IFR = concatIFR;
%%  Run Training and Testing Curve Fitting
% This section of the code will run 100 times to ensure confidence in the
% resulting fitting statistics. 
% * First, a random permutation is generated at the same size as the 
%   concatenated data for the afferent. 
% * Then, the training and testing indeces are specified for each set, and 
%   the data structures for each are populated accordingly.
% * The gains that minimize the cost function specified in 


    for loop = 1:100
        clc
        disp([num2str(affnums(aff)) 'loop: ' num2str(loop)])
        random_index = randperm(numObs); % Generate random permutation of 1 to number of observations
        
        training_index.ind = [];
        test_index.ind = [];
        
        for i = 1:floor(0.25*numel(random_index))
            ind = (PertStart(random_index(i)):PertEnd(random_index(i)))';
            training_index.ind = vertcat(training_index.ind,ind);
        end
        for i = floor(0.25*numel(random_index))+1:numel(random_index)
            ind = (PertStart(random_index(i)):PertEnd(random_index(i)))';
            test_index.ind = vertcat(test_index.ind,ind);
        end
        
        
  
        traindata.Length = concatLength(training_index().ind);
        traindata.Velocity = concatVel(training_index().ind);
        traindata.dVdt = concatAccel(training_index().ind);
        traindata.Force = concatForce(training_index().ind);
        traindata.dFdt = concatdFdt(training_index().ind);
        traindata.LengthF = concatLengthF(training_index().ind);
        traindata.VelocityF = concatVelocityF(training_index().ind);
        traindata.dVFdt = concatdVFdt(training_index().ind);
        traindata.LengthF_hc = concatLengthF_hc(training_index().ind);
        traindata.VelocityF_hc = concatVelocityF_hc(training_index().ind);
        traindata.dVFdt_hc = concatdVFdt_hc(training_index().ind);
        traindata.firing_rate = concatIFR(training_index().ind);
        
        for model = 1:6
            flags.model = model;
            switch flags.model
                case 1 % Force, dFdt, static model
                    gain_limits(1,:) = [1,1,-65,41];                    %Initial guesses for [k_F,k_dF,Static,F_offset,dF_offset]
                    gain_limits(2,:) = [0.0,0.0,-500,-500];                   %Lower bounds for [k_F,k_dF,Static,F_offset,dF_offset]
                    gain_limits(3,:) = [2,2,1000,1000];         %Upper bounds [k_F,k_dF,Static,F_offset,dF_offset]
                    
                    %Exceptions
%                     if affnums(aff) == 58
%                         gain_limits(2,1) = 1.2; %Force the lower bound for k_F = 1. (Discrepancy between individual trial fits and xval). 
%                         gain_limits(3,1) = 1.5;
%                     end
                    
                case 2 % Kinematic Model (dV/dt,Vel,Len)
                    gain_limits(1,:) = [0.001,0.001,0.001,1,1,1];                       %Initial guesses for [k_p,k_v,k_a,p_offset,v_offset,a_offset]
                    gain_limits(2,:) = [0,0,0,-1000,-1000,-1000];           %Lower bounds for [k_p,k_v,k_a,p_offset,v_offset,a_offset]
                    gain_limits(3,:) = [1000,100,10,1000,1000,1000];  %Upper bounds for [k_p,k_v,k_a,p_offset,v_offset,a_offset]
                    
%                 case 3 % Classic Kinematic Model (Vel,Len)
%                     gain_limits(1,:) = [1,1,1,1];                           %Initial guesses for [k_p,k_v,p_offset,v_offset]
%                     gain_limits(2,:) = [0,0,-1000,-1000];                   %Lower bounds for [k_p,k_v,p_offset,v_offset]
%                     gain_limits(3,:) = [10000,10000,1000,1000];             %Upper bounds for [k_p,k_v,p_offset,v_offset]
                    
                case 4 % Fascicle Estimates Model (dVF/dt, VelF, LenF)
                    gain_limits(1,:) = [1,1,1,1,1,1];                       %Initial guesses for [k_p,k_v,k_a,p_offset,v_offset,a_offset]
                    gain_limits(2,:) = [0,0,0,-1000,-1000,-1000];           %Lower bounds for [k_p,k_v,k_a,p_offset,v_offset,a_offset]
                    gain_limits(3,:) = [10000,10000,10000,1000,1000,1000];  %Upper bounds for [k_p,k_v,k_a,p_offset,v_offset,a_offset]
                    
                case 5 % Fascicle Estimates Model w/ highly compliant tendon (dVF_hc/dt, VelF_hc, LenF_hc)
                    gain_limits(1,:) = [1,1,1,1,1,1];                       %Initial guesses for [k_p,k_v,k_a,p_offset,v_offset,a_offset]
                    gain_limits(2,:) = [0,0,0,-1000,-1000,-1000];           %Lower bounds for [k_p,k_v,k_a,p_offset,v_offset,a_offset]
                    gain_limits(3,:) = [10000,10000,10000,1000,1000,1000];  %Upper bounds for [k_p,k_v,k_a,p_offset,v_offset,a_offset]
                    
                case 6 % Free regression with all variables
                    gain_limits(1,:) = [ones(1,25)];             %Initial guesses: [kF,kdF,kp,kv,ka,kpf,kvf,kaf,kpfhc,kpvhc,kpahc,static,F_offset,p_offset,pF_offset,pFhc_offset]
                    gain_limits(2,:) = [zeros(1,12),-1000*ones(1,12),0];       %Lower Bounds
                    gain_limits(3,:) = [1e4*ones(1,12),1e3*ones(1,12),1];    %Upper bounds
                    
                case 3 % Prochazka model
                    gain_limits(1,:) = [0.5*ones(1,4)];
                    gain_limits(2,:) = [zeros(1,4)];
                    gain_limits(3,:) = [1000,1,1000,1000];
            end
            flags.weights = [1; 10];
            
            
            [train_MLEs,train_fval,train_SSE,train_ybar,train_SSM,train_R_squared,train_exitflag,train_output] = spindle_findGains(gain_limits,traindata,flags);
            xvaldata.MLEout(model).train_MLEs(loop,:) = train_MLEs;
            xvaldata.MLEout(model).train_fval(loop) = train_fval/n;
            xvaldata.MLEout(model).train_R_squared(loop) = train_R_squared;
            xvaldata.MLEout(model).train_exitflag(loop) = train_exitflag;
            xvaldata.MLEout(model).train_output(loop) = train_output;
            xvaldata.MLEout(model).train_n(loop) = floor(0.25*numel(random_index));
            xvaldata.MLEout(model).train_nspikes(loop) = numel(traindata.firing_rate);
            xvaldata.MLEout(model).train_numParam = numParam(model);
            xvaldata.MLEout(model).train_SSE(loop) = train_SSE;
            xvaldata.MLEout(model).train_ybar(loop) = train_ybar;
            xvaldata.MLEout(model).train_SSM(loop) = train_SSM;
            
            xvaldata.MLEout(model).train_logL(loop) = log(1 - train_SSE/train_SSM);
            [xvaldata.MLEout(model).train_AIC(loop),xvaldata.MLEout(model).train_BIC(loop)] = aicbic(xvaldata.MLEout(model).train_logL(loop),numParam(model),xvaldata.MLEout(model).train_n(loop));
            xvaldata.MLEout(model).train_AICc(loop) = xvaldata.MLEout(model).train_AIC(loop) + 2*numParam(model)*(numParam(model) + 1)/(xvaldata.MLEout(model).train_n(loop) - numParam(model) - 1);
            
            % TEST
            test_gains = train_MLEs;
            testdata.Length = concatLength(test_index().ind);
            testdata.Velocity = concatVel(test_index().ind);
            testdata.dVdt = concatAccel(test_index().ind);
            testdata.Force = concatForce(test_index().ind);
            testdata.dFdt = concatdFdt(test_index().ind);
            testdata.LengthF = concatLengthF(test_index().ind);
            testdata.VelocityF = concatVelocityF(test_index().ind);
            testdata.dVFdt = concatdVFdt(test_index().ind);
            testdata.LengthF_hc = concatLengthF_hc(test_index().ind);
            testdata.VelocityF_hc = concatVelocityF_hc(test_index().ind);
            testdata.dVFdt_hc = concatdVFdt_hc(test_index().ind);
            testdata.firing_rate = concatIFR(test_index().ind);
            
            % 
           
            
            switch model

                   
                case 1 % Force Model
                    opt_fit = kinetics(testdata,test_gains,flags);
                case 2   % Kinematic Model with Length and Velocity and Accel.
                    opt_fit = kinematics(testdata,test_gains,flags);
%                 case 3   % Classic Kinematics
%                     opt_fit = kinematics(testdata,test_gains,flags);
                case 4
                    opt_fit = kinematics(testdata,test_gains,flags);
                case 5
                    opt_fit = kinematics(testdata,test_gains,flags);
                case 6
                    opt_fit = mixedKin(testdata,test_gains);                    
                case 3 
                    opt_fit = prochazka(testdata,test_gains);
            end
            
            FR_recorded = testdata.firing_rate;
            R = corrcoef(opt_fit,FR_recorded);
            if numel(R) == 1
                xvaldata.MLEout(model).test_R_squared(loop) = R^2;
            else
                xvaldata.MLEout(model).test_R_squared(loop) = R(1,2)^2;
            end
            
            xvaldata.MLEout(model).test_SSE(loop) = sum((FR_recorded-opt_fit).^2);
            xvaldata.MLEout(model).test_ybar(loop) = mean(FR_recorded);
            xvaldata.MLEout(model).test_SSM(loop) = sum((FR_recorded-xvaldata.MLEout(model).test_ybar(loop)).^2);

            test_fval = spindle_cost(test_gains,testdata,flags);
            xvaldata.MLEout(model).test_fval(loop) = test_fval/n;
            xvaldata.MLEout(model).test_n(loop) = numel(random_index) - (floor(0.25*numel(random_index))+1);
            xvaldata.MLEout(model).test_nspikes(loop) = numel(testdata.firing_rate);
            xvaldata.MLEout(model).test_logL(loop) = log(1/xvaldata.MLEout(model).test_SSE(loop)/xvaldata.MLEout(model).test_SSM(loop));
            [xvaldata.MLEout(model).test_AIC(loop),xvaldata.MLEout(model).test_BIC(loop)] = aicbic(xvaldata.MLEout(model).test_logL(loop),numParam(model),xvaldata.MLEout(model).test_n(loop));
            xvaldata.MLEout(model).test_AICc(loop) = xvaldata.MLEout(model).test_AIC(loop) + 2*numParam(model)*(numParam(model) + 1)/(xvaldata.MLEout(model).test_n(loop) - numParam(model) - 1);
            
            clear gain_limits 
        end
        clear traindata testdata train_MLEs test_gains
    end
 %% Saves all data used in  cross-validation to a single .mat file   
    savename = [pathname filesep affname 'xvaldata.mat'];
    if exist(savename,'file')
        delete(savename)
    end
    save(savename,'xvaldata');
    
 %% Set up summary array for all afferents
 for model = 1:6
 % Training Set
    xvalDataAll.train.AICdata.mean(model,aff) = mean(xvaldata.MLEout(model).train_AIC);
    xvalDataAll.train.AICcdata.mean(model,aff) = mean(xvaldata.MLEout(model).train_AICc);
    xvalDataAll.train.BICdata.mean(model,aff) = mean(xvaldata.MLEout(model).train_BIC);
    xvalDataAll.train.R2data.mean(model,aff) = mean(xvaldata.MLEout(model).train_R_squared);
    xvalDataAll.train.SSEdata.mean(model,aff) = mean(xvaldata.MLEout(model).train_SSE);
    xvalDataAll.train.ybardata.mean(model,aff) = mean(xvaldata.MLEout(model).train_ybar);
    xvalDataAll.train.AICdata.stdev(model,aff) = std(xvaldata.MLEout(model).train_AIC);
    xvalDataAll.train.AICcdata.stdev(model,aff) = std(xvaldata.MLEout(model).train_AICc);
    xvalDataAll.train.BICdata.stdev(model,aff) = std(xvaldata.MLEout(model).train_BIC);
    xvalDataAll.train.R2data.stdev(model,aff) = std(xvaldata.MLEout(model).train_R_squared);
    xvalDataAll.train.SSEdata.stdev(model,aff) = std(xvaldata.MLEout(model).train_SSE);
    xvalDataAll.train.ybardata.stdev(model,aff) = std(xvaldata.MLEout(model).train_ybar);
    xvalDataAll.train.ndata.numperts(model,aff) = xvaldata.MLEout(model).train_n(1); % This is constant
    xvalDataAll.train.nspikedata.mean(model,aff) = mean(xvaldata.MLEout(model).train_nspikes);
    xvalDataAll.train.nspikedata.stdev(model,aff) = std(xvaldata.MLEout(model).train_nspikes);
    
    
 % Testing Set   
    xvalDataAll.test.AICdata.mean(model,aff) = mean(xvaldata.MLEout(model).test_AIC);
    xvalDataAll.test.AICcdata.mean(model,aff) = mean(xvaldata.MLEout(model).test_AICc);
    xvalDataAll.test.BICdata.mean(model,aff) = mean(xvaldata.MLEout(model).test_BIC);
    xvalDataAll.test.R2data.mean(model,aff) = mean(xvaldata.MLEout(model).test_R_squared);
    xvalDataAll.test.SSEdata.mean(model,aff) = mean(xvaldata.MLEout(model).test_SSE);
    xvalDataAll.test.ybardata.mean(model,aff) = mean(xvaldata.MLEout(model).test_ybar);
    xvalDataAll.test.AICdata.stdev(model,aff) = std(xvaldata.MLEout(model).test_AIC);
    xvalDataAll.test.AICcdata.stdev(model,aff) = std(xvaldata.MLEout(model).test_AICc);
    xvalDataAll.test.BICdata.stdev(model,aff) = std(xvaldata.MLEout(model).test_BIC);
    xvalDataAll.test.R2data.stdev(model,aff) = std(xvaldata.MLEout(model).test_R_squared);
    xvalDataAll.test.SSEdata.stdev(model,aff) = std(xvaldata.MLEout(model).test_SSE);
    xvalDataAll.test.ybardata.stdev(model,aff) = std(xvaldata.MLEout(model).test_ybar);
    xvalDataAll.test.ndata.numperts(model,aff) = xvaldata.MLEout(model).test_n(1); % This is constant
    xvalDataAll.test.nspikedata.mean(model,aff) = mean(xvaldata.MLEout(model).test_nspikes);
    xvalDataAll.test.nspikedata.stdev(model,aff) = std(xvaldata.MLEout(model).test_nspikes);
    
 % Gain Data

    
    
    if model == 6 % For the last loop
        % Train AIC
        xvalDataAll.train.AICdata.delta(:,aff) = ...
            xvalDataAll.test.AICdata.mean(:,aff) - min(xvalDataAll.test.AICdata.mean(:,aff));
        xvalDataAll.train.AICdata.L(:,aff) = exp(-1/2*xvalDataAll.train.AICdata.delta(:,aff));
        xvalDataAll.train.AICdata.w(:,aff) = xvalDataAll.train.AICdata.L(:,aff)./sum(xvalDataAll.train.AICdata.L(:,aff));
        % Train AICc
        xvalDataAll.train.AICcdata.delta(:,aff) = ...
            xvalDataAll.test.AICcdata.mean(:,aff) - min(xvalDataAll.test.AICcdata.mean(:,aff));
        xvalDataAll.train.AICcdata.L(:,aff) = exp(-1/2*xvalDataAll.train.AICcdata.delta(:,aff));
        xvalDataAll.train.AICcdata.w(:,aff) = xvalDataAll.train.AICcdata.L(:,aff)./sum(xvalDataAll.train.AICcdata.L(:,aff));
        % Train BIC
        xvalDataAll.train.BICdata.delta(:,aff) = ...
            xvalDataAll.test.BICdata.mean(:,aff) - min(xvalDataAll.test.BICdata.mean(:,aff));
        xvalDataAll.train.BICdata.L(:,aff) = exp(-1/2*xvalDataAll.train.BICdata.delta(:,aff));
        xvalDataAll.train.BICdata.w(:,aff) = xvalDataAll.train.BICdata.L(:,aff)./sum(xvalDataAll.train.BICdata.L(:,aff));
       
        % Test AIC
        xvalDataAll.test.AICdata.delta(:,aff) = ...
            xvalDataAll.test.AICdata.mean(:,aff) - min(xvalDataAll.test.AICdata.mean(:,aff));
        xvalDataAll.test.AICdata.L(:,aff) = exp(-1/2*xvalDataAll.test.AICdata.delta(:,aff));
        xvalDataAll.test.AICdata.w(:,aff) = xvalDataAll.test.AICdata.L(:,aff)./sum(xvalDataAll.test.AICdata.L(:,aff));
        % Test AICc
        xvalDataAll.test.AICcdata.delta(:,aff) = ...
            xvalDataAll.test.AICcdata.mean(:,aff) - min(xvalDataAll.test.AICcdata.mean(:,aff));
        xvalDataAll.test.AICcdata.L(:,aff) = exp(-1/2*xvalDataAll.test.AICcdata.delta(:,aff));
        xvalDataAll.test.AICcdata.w(:,aff) = xvalDataAll.test.AICcdata.L(:,aff)./sum(xvalDataAll.test.AICcdata.L(:,aff));
        % Test BIC
        xvalDataAll.test.BICdata.delta(:,aff) = ...
            xvalDataAll.test.BICdata.mean(:,aff) - min(xvalDataAll.test.BICdata.mean(:,aff));
        xvalDataAll.test.BICdata.L(:,aff) = exp(-1/2*xvalDataAll.test.BICdata.delta(:,aff));
        xvalDataAll.test.BICdata.w(:,aff) = xvalDataAll.test.BICdata.L(:,aff)./sum(xvalDataAll.test.BICdata.L(:,aff));

        
        
        xvalDataAll = compileMLEmeans(xvalDataAll);

    end
 end
    
    clearvars -except aff affnums xvalDataAll

end
savepath = cd(cd(['..' filesep 'Data' filesep 'stats']));
savename = [savepath filesep 'xvalDataAll' date '.mat'];
if exist(savename,'file')
    delete(savename)
end
save(savename,'xvalDataAll');
toc
clearvars