function new_xe = update_xe(M, x_r, x_e, tau_re, tau_ei)
    new_xe = x_e + (M*(x_r/tau_re)) - (x_e/tau_ei);
end

