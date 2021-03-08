function f = kinetics(data,gains,flags)

F = data.Force;
dF = data.dFdt;

% Cancel negative dF/dt values, assume remaining dF/dt values dominate F
% if flags.aff == 25 || flags.aff == 26 || flags.aff == 58 || flags.aff == 60
%     dF(dF<0) = 0;
%     F(dF~=0) = 0;
% end

% Assign gains and offsets
F = gains(1)*F + gains(3);
dF = gains(2)*dF + gains(4);

% Intermediate step for half-wave rectification for resultant F, dF/dt components
F(F <= gains(3)) = 0;
dF(dF <= gains(4)) = 0;

% dF/dt dominance over F, half-wave rectify
% if flags.aff == 25 || flags.aff == 26 || flags.aff == 58 || flags.aff == 60
%     F(dF > 0) = 0;
%     F(F <= 0) = 0;
%     dF(dF <= 0) = 0;
% end

f = F + dF;