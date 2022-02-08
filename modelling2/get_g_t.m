function g_t = get_g_t(spike_train, dt, t, tspan)
    kernel_kt = shift_1(get_k_t(spike_train));
	conv_vector = conv(spike_train, kernel_kt);
    g_t = conv_vector(1, t);
end
