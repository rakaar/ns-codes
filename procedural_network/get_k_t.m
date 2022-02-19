function k_t = get_k_t(spike_train)
    t_spike = 0;
    k_t = zeros(1,length(spike_train));
    for i=1:length(spike_train)
        if spike_train(1,i) == 1
            t_spike = i;
        end
        k_t(1,i) = exp(-(i-t_spike)/10);
    end    
end