function f = prochazka(data,gains)

f = gains(1).*abs(data.Velocity).^(gains(2)) + gains(3)*data.Length + gains(4);