tic

affs = [24,25,26,27,32,56,57,58,60,64];
savefile = 'all_affs_stat_array_no_buffer_lags2';
savepath = cd(cd(['..' filesep 'Data' filesep 'stats' filesep]));
savename = [savepath filesep savefile];
fid = fopen([savename '.csv'],'wt');
stat_array = {'affnum','model','pert_type','trial','R_squared',...
    'kF','kdF','o_F','o_dF','ka','kv','kl','o_a','o_v','o_l','bgF',...
    'bgIFR','maxF','maxIFR',...
    'maxL','maxV','maxA','maxdF','AdjR_squared','Lag','AIC','BIC','AICc','J','DI_1s','DI_500ms'};

for affnum = 1:length(affs) 
    affname = ['aff' num2str(affs(affnum))];
    %     pathname = cd(cd(['..' filesep 'Data' filesep 'proc_data_dynamic' filesep affname '_proc']));
    pathname = cd(cd(['..' filesep 'Data' filesep 'proc_data_no_buffer_lags' filesep affname '_proc']));
    directory = dir(pathname);
    for trial = 1:length(directory) % For all files in directory
        filename = directory(trial).name; % Current file
        if strfind(filename,'proc.mat')    % If current file is .mat file
            loadname = [pathname filesep filename];
            load(loadname)
            disp(loadname)
            pert_type = info.trial_type;
            for model = 1:length(fit_data) % For each model
                for pert = 1:size(fit_data(model).opt_gains,1)
                    
                    lagtimes = 0:0.001:0.015;
                    [R_squared,lag] = max(fit_data(model).R_squared(:,pert));
                    lagtime = lagtimes(lag);
                    bgF = mean(proc_data(pert).Force(1));
                    bgIFR = mean(proc_data(pert).firing_rate(1:3));
                    maxF = max(proc_data(pert).Force);
                    maxL = max(proc_data(pert).Length);
                    maxV = max(proc_data(pert).Velocity);
                    maxA = max(proc_data(pert).dVdt);
                    maxdF = max(proc_data(pert).dFdt);
                    maxIFR = max(proc_data(pert).firing_rate);
                    [DI1,~] = dynamicIndex(proc_data(pert).firing_rate,proc_data(pert).spiketimes,proc_data(pert).time,0.5);
                    [DI2,~] = dynamicIndex(proc_data(pert).firing_rate,proc_data(pert).spiketimes,proc_data(pert).time,0.0);
                    
                    
                    
                    if model == 1
                        k_F = fit_data(model).opt_gains(pert,lag,1);
                        k_dF = fit_data(model).opt_gains(pert,lag,2);
                        o_F = fit_data(model).opt_gains(pert,lag,3);
                        o_dF = fit_data(model).opt_gains(pert,lag,4);
                        k_a = 0;
                        k_v = 0;
                        k_l = 0;
                        o_a = 0;
                        o_v = 0;
                        o_l = 0;
                    elseif model == 2
                        k_F = 0;
                        k_dF = 0;
                        o_F = 0;
                        o_dF = 0;
                        k_a = fit_data(model).opt_gains(pert,lag,3);
                        k_v = fit_data(model).opt_gains(pert,lag,2);
                        k_l = fit_data(model).opt_gains(pert,lag,1);
                        o_a = fit_data(model).opt_gains(pert,lag,6);
                        o_v = fit_data(model).opt_gains(pert,lag,5);
                        o_l = fit_data(model).opt_gains(pert,lag,4);
                    elseif model == 3
                        k_F = 0;
                        k_dF = 0;
                        o_F = 0;
                        o_dF = 0;
                        k_a = 0;
                        k_v = fit_data(model).opt_gains(pert,lag,2);
                        k_l = fit_data(model).opt_gains(pert,lag,1);
                        o_a = 0;
                        o_v = fit_data(model).opt_gains(pert,lag,4);
                        o_l = fit_data(model).opt_gains(pert,lag,3);
                    elseif model == 4
                        k_F = 0;
                        k_dF = 0;
                        o_F = 0;
                        o_dF = 0;
                        k_a = fit_data(model).opt_gains(pert,lag,3);
                        k_v = fit_data(model).opt_gains(pert,lag,2);
                        k_l = fit_data(model).opt_gains(pert,lag,1);
                        o_a = fit_data(model).opt_gains(pert,lag,6);
                        o_v = fit_data(model).opt_gains(pert,lag,5);
                        o_l = fit_data(model).opt_gains(pert,lag,4);
                    elseif model == 5;
                        k_F = 0;
                        k_dF = 0;
                        o_F = 0;
                        o_dF = 0;
                        k_a = fit_data(model).opt_gains(pert,lag,3);
                        k_v = fit_data(model).opt_gains(pert,lag,2);
                        k_l = fit_data(model).opt_gains(pert,lag,1);
                        o_a = fit_data(model).opt_gains(pert,lag,6);
                        o_v = fit_data(model).opt_gains(pert,lag,5);
                        o_l = fit_data(model).opt_gains(pert,lag,4);
                    elseif model == 6
                        k_F = fit_data(model).opt_gains(pert,lag,1);
                        k_dF = fit_data(model).opt_gains(pert,lag,2);
                        o_F = fit_data(model).opt_gains(pert,lag,6);
                        o_dF = fit_data(model).opt_gains(pert,lag,7);
                        k_a = fit_data(model).opt_gains(pert,lag,5);
                        k_v = fit_data(model).opt_gains(pert,lag,4);
                        k_l = fit_data(model).opt_gains(pert,lag,3);
                        o_a = fit_data(model).opt_gains(pert,lag,10);
                        o_v = fit_data(model).opt_gains(pert,lag,9);
                        o_l = fit_data(model).opt_gains(pert,lag,8);
                    elseif model == 7
                        k_F = 0;
                        k_dF = 0;
                        o_F = 0;
                        o_dF = 0;
                        k_a = 0;
                        k_v = fit_data(model).opt_gains(pert,lag,1);
                        k_l = fit_data(model).opt_gains(pert,lag,3);
                        o_a = 0;
                        o_v = fit_data(model).opt_gains(pert,lag,2);
                        o_l = fit_data(model).opt_gains(pert,lag,4);
                    end
                   
                    % Adjusted R squared & Information Criterion
                    J = fit_data(model).terminal_cost(lag,pert);
                    logL = log(1/J);
                    numParam = numel(fit_data(model).opt_gains(pert,lag,:)) + 2;  % Estimated Variance + lag are necessarily included as a parameter
                    numObs = numel(proc_data(pert).firing_rate);
                    [AIC,BIC] = aicbic(logL,numParam,numObs);
                    AICc = AIC + (2 * numParam * (numParam + 1))/(numObs - numParam - 1);
                    AdjR_squared = R_squared - (1 - R_squared)*(numParam/(numObs - numParam - 1));
                    

                    
                    %%% Generate stat_array %%%
                    stat_array{end+1,1} = affs(affnum);         % Also generates empty cells for rest of column
                    stat_array{end,2} = model;
                    stat_array{end,3} = pert_type;
                    stat_array{end,4} = str2double(filename(11:12));
                    stat_array{end,5} = R_squared;
                    stat_array{end,6} = k_F;
                    stat_array{end,7} = k_dF;
                    stat_array{end,8} = o_F;
                    stat_array{end,9} = o_dF;
                    stat_array{end,10} = k_a;
                    stat_array{end,11} = k_v;
                    stat_array{end,12} = k_l;
                    stat_array{end,13} = o_a;
                    stat_array{end,14} = o_v;
                    stat_array{end,15} = o_l;
                    stat_array{end,16} = bgF;
                    stat_array{end,17} = bgIFR;
                    stat_array{end,18} = maxF;
                    stat_array{end,19} = maxIFR;
                    stat_array{end,20} = maxL;
                    stat_array{end,21} = maxV;
                    stat_array{end,22} = maxA;
                    stat_array{end,23} = maxdF;
                    stat_array{end,24} = AdjR_squared;
                    stat_array{end,25} = lagtime;
                    stat_array{end,26} = AIC;
                    stat_array{end,27} = BIC;
                    stat_array{end,28} = AICc;
                    stat_array{end,29} = J;
                    stat_array{end,30} = DI1;
                    stat_array{end,31} = DI2;
                    
                    
                end
            end
        end
    end
end

for i=1:size(stat_array,1)
    if i == 1
        fprintf(fid,[repmat('%s,',[1,size(stat_array,2)-1]),'%s\n'],stat_array{i,:});
    else
        fprintf(fid, ['%f,%f,%s,',repmat('%f,',[1,size(stat_array,2)-4]),'%f\n'], stat_array{i,:});
    end
end
fclose(fid);

spindle_getMeanLags([cd(cd(['..' filesep 'Data' filesep 'stats' filesep])) filesep 'all_affs_stat_array_no_buffer_lags2.csv'])

toc