% for a random neuron
col = 1;
for random_neuron=1:n_total_neurons
    clf
figure(random_neuron)
    hold on
        mean_spike_rate_of_random_neuron = zeros(1, spike_rate_length);
        for t=1:spike_rate_length
              mean_spike_rate_of_random_neuron(1,t) = sum(spike_rates(:, col, random_neuron,t))/n_iters;
        end
        plot(mean_spike_rate_of_random_neuron);

        % adjusting from length(tspan)-1 to length(tspan)
        epsc_all_iters_for_random_neuron = squeeze(total_input_epsc(:, col, random_neuron, :));
        actual_epsc_input_to_random_neuron_size_adjusted = zeros(n_iters, length(tspan));
        for iter=1:n_iters
            actual_epsc_input_to_random_neuron_size_adjusted(iter, 2:length(tspan)) = epsc_all_iters_for_random_neuron(iter, :);
        end
        
        % taking avg iters
        random_epsc_iters_avg = zeros(1, length(tspan));
        for t=1:length(tspan)
              random_epsc_iters_avg(1,t) = sum(actual_epsc_input_to_random_neuron_size_adjusted(:,t))/n_iters;
        end
        
        epsc_to_random_neuron_binned = spikes_to_spike_rate_neat(random_epsc_iters_avg, physical_time_in_ms, dt, spike_rate_dt);
        plot(epsc_to_random_neuron_binned*(n_bins*physical_time_in_ms*0.001))
            
    hold off
        title('random neuron epsc and psth')
        legend('psth', 'epsc')
grid
pause
clf
end