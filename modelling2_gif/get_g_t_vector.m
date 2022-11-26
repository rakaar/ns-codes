function g_t_vector = get_g_t_vector(spike_train, period)
    tau_syn = 10;
    spike_train = reshape(spike_train, 1,length(spike_train));
	kernel_kt = [0 exp(-[0:length(spike_train)-1])./tau_syn];
    g_t_vector = conv(spike_train, kernel_kt);
    g_t_vector = g_t_vector(1, 1:period);
end