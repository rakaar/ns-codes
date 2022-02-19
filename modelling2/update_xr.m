function new_xr = update_xr(M, x_r, x_i, tau_re, tau_ir)
    new_xr = x_r + (-M*(x_r/tau_re)) + (x_i/tau_ir);
end
