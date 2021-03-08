% Created 06/16/2016
% Author: KPB
% This function loads a .csv file from the 'stats' data folder created by
% the spindle_statFormat_lags.m script, which selected the lagtime with the
% highest associated R2 value. This function does not return any values,
% but instead saves the 'mean_lags.mat' file, which is loaded and used
% by the 'spindle_operate' script to fit the models. 

function spindle_getMeanLags(filename)
A = readtable(filename);
affs = [24;25;26;27;32;56;57;58;60;64];
mean_lags = zeros(10,8);
mean_lags(:,1) = affs;
for i = 1:10
    for j = 2:8
        mean_lags(i,j) = mean(A.Lag(A.affnum==affs(i) & A.model==(j-1)));
    end
end

savepathname = cd(cd(['..' filesep 'Data' filesep 'stats']));
save([savepathname filesep 'mean_lags2.mat'],'mean_lags')
end