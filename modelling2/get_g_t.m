function g_t = get_g_t(spike_train, dt, t, tspan)
    % kernel_kt = shift_1(get_k_t(spike_train));
    tau_syn = 10;
	kernel_kt = [0 exp(-[0:length(spike_train)-1])./tau_syn];
    conv_vector = conv(spike_train, kernel_kt);
    g_t = conv_vector(1, t);
end
