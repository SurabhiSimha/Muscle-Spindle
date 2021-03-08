% Author: Kyle P. Blum
% Date of last edit: June, 12 2015
%
% Description: This script will parse through all datacl files contained in
% '../Data' and decide whether to process it or not. User can specify which
% afferents to process data for in the 'good_affs' cell array. 
%
% User-defined dependencies: 'septrials.m', 'findspikes.m'

tic
clear

pathname = cd(cd(['..' filesep 'Data' filesep]));
proc_dataName = 'proc_data_all_afferents';
% good_affs = {'24','25','26','27','32','56','57','58','60','64'};
% good_affs =
% {'25','26','58','60','12','13','15','17','23','31','38','51','54'}; % %Potential IB afferents
good_affs = {'12','13','15','17','23','24','25','26','27','32','56','57','58','60','64'};

main_dir = dir(pathname);

trial_types = {'accel_ramp';...
    'accel_series';...
    'trad_ramp';...
    'trad_series';...
    'triangle_series';...
    'pyramid';...
    'zap'};                   % List of trial types

for i = 1:length(main_dir)                   % All items in pathname
    folder_name = main_dir(i).name;          % Current item
    if strfind(folder_name,'20')>0;        % If current item's name contains '20', e.g. '2009/'
        sub_dir = dir([pathname filesep folder_name]);    % Use this folder
        for j = 1:length(sub_dir)                   % For all files in current item
            filename = sub_dir(j).name;               % Current file
            if strfind(filename,'.mat')               % If current file contains '.mat'
                for k = 1:length(good_affs)              % For all good afferents
                    affname = ['aff' num2str(good_affs{k}) '_'];  % Current afferent name
                    %                     if length(affname) == 4;
                    %                         affname = ['aff' '0' num2str(good_affs(k))];
                    %                     end
                    if strfind(filename,affname)              % If current file is from current afferent
                        if strfind(filename,'descente')          % If current file also contains 'descente'
                            break                                 % Break out (want to ignore descente files)
                        else                                     % If not...
%                             if strfind(filename,'triangle')
                            
                            loadname = [pathname filesep folder_name filesep filename];    % Load this file
                            
                            [num_perts,ignore,proc_data] = septrials(loadname);   % Separate perturbations
%                             [num_perts,ignore,proc_data] = rat_septrials(loadname);   % Separate rat perturbations
                            if length(affname) == 4
                                savefolder = [pathname filesep proc_dataName filesep affname(1:3) '0' affname(4) 'proc' filesep];% Save them to this folder
                            else
%                                 savefolder = [pathname filesep 'proc_data' filesep affname 'proc' filesep];% Save them to this folder
                                savefolder = [pathname filesep proc_dataName filesep affname 'proc' filesep];% Save them to this folder
                            end
                            if ignore == 0;                                       % If we don't ignore this file
                                if ~exist(savefolder,'dir')                             % If folder doesn't exist...
                                    mkdir(savefolder)                                   % Make it
                                end
                                if filename(12) == '_'
                                    filename = [filename(1:10) '0' filename(11:end)];
                                end
                                if filename(5) == '_'
                                    filename = [filename(1:3) '0' filename(4:end)];
                                end
                                
                                savename = [savefolder filename(1:end-4) '_proc.mat']; % Save new file as this
                                disp(savename)
                                save(savename,'proc_data')
                                proc_data = findspikes(savename,num_perts);           % Use spike-detection algorithm
                                save(savename,'proc_data','-append');                 % Append new variables
                                for n = 1:length(trial_types)                         % For all trial types
                                    if strfind(savename,trial_types{n})               % If this file is current trial type
                                        info.trial_type = trial_types{n};                       % Set 'trial_type'
                                    end
                                end
                                info.name = [filename(1:end-4) '_proc.mat'];          % Name of new file
                                info.date = sub_dir(j).date;                          % Use folder to find date
                                info.affname = savefolder(end-10:end-6);
                                
                                save(savename,'info','-append')                       % Append these variables
                            end
%                             end
                        end
                    end
                end
            end
        end
    end
end


toc




