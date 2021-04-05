function x = find_hsl_from_force(obj,F_isotonic)
x_adj = fzero(@(x) tempForce(x,obj,F_isotonic), 0, optimset('display','off'));



x = obj.hs_length + x_adj;
end


function zeroF = tempForce(x,obj,F_isotonic)
% Adjust for filament compliance
delta_x = x * obj.compliance_factor;
% Shift populations by interpolation
interp_positions = obj.x_bins - delta_x;
% THIS IS WHERE YOU SOLVE FOR BIN POPULATIONS %
temp_bin_pops = interp1(obj.x_bins,obj.bin_pops,interp_positions, ...
    'linear',0)';
cbF = obj.cb_number_density * obj.k_cb * 1e-9 * ...
    sum((obj.x_bins + obj.power_stroke).* temp_bin_pops');

pF = obj.k_passive * (obj.hs_length - obj.hsl_slack);
hsF = cbF + pF;


zeroF = (hsF - F_isotonic);

end