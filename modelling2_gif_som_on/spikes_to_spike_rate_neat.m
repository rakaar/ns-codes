function spike_rate = spikes_to_spike_rate_neat(spikes, physical_time_in_ms, dt, spike_rate_dt)
    spikes = reshape(spikes, 1, length(spikes));
    n_bins = spike_rate_dt/dt;
    original_length_spikes = length(spikes);
    % original_length_spikes-1 to make it even
    spikes = reshape(spikes(1, 1:original_length_spikes-1), n_bins,(original_length_spikes-1)/n_bins);
    n_cols = size(spikes,2);

    spike_rate = zeros(1, (original_length_spikes-1)/n_bins);
    for c=1:n_cols
%         spike_rate(1, c) = sum(spikes(:, c))/(n_bins*physical_time_in_ms*0.001);
          spike_rate(1, c) = sum(spikes(:, c))/n_bins;
    end
end