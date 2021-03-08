function ibData = ibRegress(ibData)
%Author: KPB
%Date created: 04/26/2017
%Description: This function takes the 'ibData' struct from get_ibData.m and
%uses regression to determine the best fit linear function of dF/dt or
%accel. peaks to the peak initial burst amplitude. Specifically, this
%function uses back-division of the dF/dt or accel. data (X) by IFR data 
%(y) to determine a gain and offset for the linear function relating X to
%y.

x1 = ibData.dFdt.pk';
X1 = [ones(size(x1)) x1];
x2 = ibData.acc.pk';
X2 = [ones(size(x2)) x2];
y = ibData.IFR.pk';

b1 = X1\y;
b2 = X2\y;

yCalc1 = X1*b1;
yCalc2 = X2*b2;

ybar = mean(y); 

SSE1 = sum((y - yCalc1).^2);
SSE2 = sum((y - yCalc2).^2);

SSM = sum((y - ybar).^2);

Rsq1 = 1 - SSE1/SSM;
Rsq2 = 1 - SSE2/SSM;

[~,p1] = corrcoef(x1,y);
[~,p2] = corrcoef(x2,y);

ibData.regressData.dFdt.b = b1;
ibData.regressData.dFdt.Rsq = Rsq1;
ibData.regressData.dFdt.p = p1(2);
ibData.regressData.dFdt.SSE = SSE1;
ibData.regressData.dFdt.SSM = SSM;

ibData.regressData.acc.b = b2;
ibData.regressData.acc.Rsq = Rsq2;
ibData.regressData.acc.p = p2(2);
ibData.regressData.acc.SSE = SSE2;
ibData.regressData.acc.SSM = SSM;



end