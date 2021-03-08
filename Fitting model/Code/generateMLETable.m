function MLEtable = generateMLETable(xvalDataAll)

for aff = 1:size(xvalDataAll.MLEout,2);
    MLEmeans(aff,:) = xvalDataAll.MLEout(1,aff).MLEs.mean;
    MLEstdevs(aff,:) = xvalDataAll.MLEout(1,aff).MLEs.stdev;
end

kFmean = MLEmeans(:,1);
bFmean = MLEmeans(:,3);
kdFmean = MLEmeans(:,2);
bdFmean = MLEmeans(:,4);
kFstd = MLEstdevs(:,1);
bFstd = MLEstdevs(:,3);
kdFstd = MLEstdevs(:,2);
bdFstd = MLEstdevs(:,4);

MLEtable = table(kFmean,bFmean,kdFmean,bdFmean,kFstd,bFstd,kdFstd,bdFstd);


RowNames = {'1','2','3','4','5','6','7','8','9','10'};
MLEtable.Properties.RowNames = RowNames;
end