function f = kinematics(data,gains,flags)

if flags.model == 2 %add 'elseif' statements for other length-based models
    L = gains(1)*data.Length + gains(4);     %Length component
    V = gains(2)*data.Velocity + gains(5);   %Velocity component
    A = gains(3)*data.dVdt + gains(6);       %Acceleration component 
    
    %Set components below offset to value of offset (removes sharp changes)
    L(L <= gains(4)) = gains(4);  
    V(V <= gains(5)) = gains(5);
    A(A <= gains(6)) = gains(6);
end
    
f = L + V + A;  % Compute firing rate from L,V,A components