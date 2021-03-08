% Created 07/13/2016
% Author: KPB
% This function is used at the end of the 'spindle_xval.m' script to append
% a new field to the 'xvalDataAll' struct. For each afferent and model,
% this function will append all MLEs (model parameters) for each of the 100
% loops, as well as the mean and standard deviation of these values. These
% data are used for plotting and statistics later. 


function xvalDataAll = compileMLEmeans(xvalDataAll)

affnums = [24,25,26,27,32,56,57,58,60,64];
for aff = 1:length(affnums)
    affname = ['aff' num2str(affnums(aff))];
    pathname = cd(cd(['..' filesep 'Data' filesep 'proc_data_no_buffer' filesep affname '_proc']));
    load([pathname filesep affname 'xvaldata.mat'])
    
    for model = 1:6
        xvalDataAll.MLEout(model,aff).MLEs.all = xvaldata.MLEout(model).train_MLEs;
        xvalDataAll.MLEout(model,aff).MLEs.mean = mean(xvaldata.MLEout(model).train_MLEs);
        xvalDataAll.MLEout(model,aff).MLEs.stdev = std(xvaldata.MLEout(model).train_MLEs);
    end
    
end