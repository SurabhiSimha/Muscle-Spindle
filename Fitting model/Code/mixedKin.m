function f = mixedKin(data,gains)

F = gains(1)*data.Force + gains(13);
dF = gains(2)*data.dFdt + gains(14);
L = gains(3)*data.Length + gains(15);
V = gains(4)*data.Velocity + gains(16);
A = gains(5)*data.dVdt + gains(17);
Lf = gains(6)*data.LengthF + gains(18);
Vf = gains(7)*data.VelocityF + gains(19);
Af = gains(8)*data.dVFdt + gains(20);
Lf_hc = gains(9)*data.LengthF_hc + gains(21);
Vf_hc = gains(10)*data.VelocityF_hc + gains(22);
Af_hc = gains(11)*data.dVFdt_hc + gains(23);
V_power = gains(12).*abs(data.Velocity).^gains(25) + gains(24);

F(F <= gains(12)) = 0;
dF(dF <= gains(13)) = 0;
L(L <= gains(14)) = 0;
V(V <= gains(15)) = 0;
A(A <= gains(16)) = 0;
Lf(Lf <= gains(17)) = 0;
Vf(Vf <= gains(18)) = 0;
Af(Af <= gains(19)) = 0;
Lf_hc(Lf_hc <= gains(20)) = 0;
Vf_hc(Vf_hc <= gains(21)) = 0;
Af_hc(Af_hc <= gains(22)) = 0;
V_power(V_power <= gains(24)) = 0;

f = F + dF + L + V + A + Lf + Vf + Af + Lf_hc + Vf_hc + Af_hc + V_power;

