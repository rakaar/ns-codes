function g_t = get_g_t(spike_train, dt, t, tspan)
    kernel_kt = shift_1(get_k_t(spike_train));
	g_t = 0;
	% t is index in form of time
	for i=1:tspan
		if t-i >= 1
			g_t = g_t + kernel_kt(1, t-i)*spike_train(1, i)*dt;
		end
	end
end
