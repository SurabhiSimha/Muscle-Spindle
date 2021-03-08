function ibData = get_ibData(proc_data,perts)
%Author: KPB
%Date created: 04/26/2017
%Description: This function takes the 'proc_data' structure for a given
%afferent and finds the initial burst amplitude and location, along with
%the peaks in accel and dF/dt and their respective locations. The data is
%returned to the user as a structure, with fields, "IFR", "dFdt", and "acc"
%which all have their own "pk" and "lc" fields corresponding to each
%variable's peak value and location. 


for i = 1:length(perts)
    pert = perts(i);
    %Sometimes the initial burst is the first IFR datapoint
    [IFRPeaks,IFRLocs] = findpeaks(proc_data(pert).firing_rate);
    if proc_data(pert).firing_rate(1) > IFRPeaks(1)
        IFRPeaks = proc_data(pert).firing_rate(1);
        IFRLocs = 1;
    end
    [dFPeaks,dFLocs] = findpeaks(proc_data(pert).dFdt,'minpeakheight',50);
    [aPeaks,aLocs] = findpeaks(proc_data(pert).dVdt,'minpeakheight',200);
    
    
    ibData.IFR.pk(i) = IFRPeaks(1);
    ibData.IFR.lc(i) = IFRLocs(1);
    ibData.IFR.t(i) = proc_data(pert).spiketimes(ibData.IFR.lc(i)+1);
    ibData.dFdt.pk(i) = dFPeaks(1);
    ibData.dFdt.lc(i) = dFLocs(1);
    ibData.dFdt.t(i) = proc_data(pert).time(ibData.dFdt.lc(i));
    ibData.acc.pk(i) = aPeaks(1);
    ibData.acc.lc(i) = aLocs(1);
    ibData.acc.t(i) = proc_data(pert).time(ibData.acc.lc(i));

end
