function spike_rate = spikes_to_spike_rate(dt, spike_rate_dt, t_simulate,physical_time_in_ms, spikes)
    % dt is unit step of spikes
    % bin size is unit step of spike rate
    % bin size/dt = number of spike 
    tspan_spike_rate = 0:spike_rate_dt:t_simulate;
    spike_rate = zeros(1, length(tspan_spike_rate));
    spike_rate_bin_size = spike_rate_dt/dt;
    spikes = reshape(spikes, 1,length(spikes));
    % starting from 0+dt, that's why index is 2 not  1
    spike_rate_bin_start_points = 2:spike_rate_bin_size:length(spikes)-(spike_rate_bin_size-1);
    for i=1:length(spike_rate_bin_start_points)
        for j=spike_rate_bin_start_points(1,i):1:spike_rate_bin_start_points(1,i)+4
            spike_rate(1,i+1) = spike_rate(1,i+1) + spikes(1,j);
        end

        spike_rate(1,i+1) = spike_rate(1,i+1)/(spike_rate_bin_size*physical_time_in_ms*0.001);
    end
    
    
    
end