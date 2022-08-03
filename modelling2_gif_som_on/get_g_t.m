function g_t = get_g_t(spike_train, dt, t, tspan, kernel_kt)
    % kernel_kt = shift_1(get_k_t(spike_train));
    
    spike_train = reshape(spike_train, 1,length(spike_train));
    if t >=  11
        x = conv(spike_train(t-10:t-1),kernel_kt);
        g_t = x(11);
    else
         x = conv(spike_train(1:t),kernel_kt);
         g_t = x(11);
    end
end
