function new_xi = update_xi(x_e, x_i, tau_ei, tau_ir)
    new_xi = x_i + (x_e/tau_ei) - (x_i/tau_ir);
end
