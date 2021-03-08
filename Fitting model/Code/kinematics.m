function f = kinematics(data,gains,flags)


if flags.model == 2
    L = gains(1)*data.Length + gains(4);
    V = gains(2)*data.Velocity + gains(5);
    A = gains(3)*data.dVdt + gains(6);
    
    L(L <= gains(4)) = 0;
    V(V <= gains(5)) = 0;
    A(A <= gains(6)) = 0;
    
elseif flags.model == 3
    L = gains(1)*data.Length + gains(3);
    V = gains(2)*data.Velocity + gains(4);
    
    L(L <= gains(3)) = 0;
    V(V <= gains(4)) = 0;

elseif flags.model == 4 
    L = gains(1)*data.LengthF + gains(4);
    V = gains(2)*data.VelocityF + gains(5);
    A = gains(3)*data.dVFdt + gains(6);
    
    L(L <= gains(4)) = 0;
    V(V <= gains(5)) = 0;
    A(A <= gains(6)) = 0;
    
    
elseif flags.model == 5
    L = gains(1)*data.LengthF_hc + gains(4);
    V = gains(2)*data.VelocityF_hc + gains(5);
    A = gains(3)*data.dVFdt_hc + gains(6);
    
    L(L <= gains(4)) = 0;
    V(V <= gains(5)) = 0;
    A(A <= gains(6)) = 0;
end

if flags.model == 3
    f = L + V;
else
    f = L + V + A;
end