function g_t = get_g_t(spike_train)
    kernel_kt = shift_1(get_k_t(spike_train));
    g_t = conv(kernel_kt, spike_train);
end
