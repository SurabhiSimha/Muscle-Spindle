% Created 08/16/2017
% Author: KPB
%
% Description: This script takes processed data from
% "proc_data_all_afferents" folder, stores data into a new struct called
% "data" and saves files for PLoS Computational Biology submission. 

% Afferents in order of PLoS Comp Bio manuscript:
affs = [58 60 26 32 25 57 56 24 27 64]; %Model fitting afferents
affs = [15 38]; %IB-only afferents (filenames will have to be changed manually)

%% Create data files for PLoS Computational Biology
% Loops cycle through folders, data files, and trials. 

for affnum = 1:length(affs)
    affname = ['aff' num2str(affs(affnum))];
    pathname = cd(cd(['..' filesep 'Data' filesep 'proc_data_all_afferents' filesep affname '_proc']));
%     pathname = cd(cd(['..' filesep 'Data' filesep 'proc_data_all_afferents' filesep affname '_proc'])); %  Have to use different pathname for IB-only afferents. They were not included in proc_data_no_fits     
    directory = dir(pathname);
    for trial = 1:length(directory) % For all files in directory
        filename = directory(trial).name; % Current file
        if strfind(filename,'proc.mat')    % If current file is .mat file
            loadname = [pathname filesep filename];
            load(loadname)
            disp(loadname)
            
            for pert = 1:numel(proc_data)
                
                data(pert).Force = proc_data(pert).Force;
                data(pert).Length = proc_data(pert).Length;
                data(pert).Velocity = proc_data(pert).Velocity; % Velocity was off by factor of 2 in proc_data_no_buffer THIS IS NOT THE CASE FOR AFF15 & 38 
                data(pert).time = proc_data(pert).time - proc_data(pert).time(1);
                data(pert).spiketimes = proc_data(pert).spiketimes;
                
                savefolder = [pathname filesep 'PLoSCompBioData' filesep];
                if ~exist(savefolder,'dir') % If folder doesn't exist...
                    mkdir(savefolder)       % Make it
                end
                
                if affnum <10 % Need to add '0' in front of affnum
                    savename = [savefolder filesep 'aff0' num2str(affnum) filename(6:end-4) '.mat']; % Save new file as this
                else
                    savename = [savefolder filesep 'aff' num2str(affnum) filename(6:end-4) '.mat']; % Save new file as this
                end
                

            end
            save(savename,'data')
            disp(savename)
            clear data
        end
    end
end

%% Turn all data files into a single file for Dryad



