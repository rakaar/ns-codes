function g_t = get_g_t(spike_train, dt, t, tspan)
    % kernel_kt = shift_1(get_k_t(spike_train));
    tau_syn = 10;
    spike_train = reshape(spike_train, 1,length(spike_train));
	kernel_kt = [0 exp(-[0:length(spike_train)-1])./tau_syn];

    if t >=  11
        x = conv(spike_train(t-10:t),kernel_kt);
        g_t = x(11);
    else
         x = conv(spike_train(1:t),kernel_kt);
         g_t = x(11);
    end
end
