% for a random neuron
col = 2;
iter_to_see=1;
n_bins = spike_rate_dt/dt;
multiply_term = (n_bins*physical_time_in_ms*0.001);

total_input_epsc = thalamic_epsc_tensor ...
                    + recurrence_exc_self_column_epsc_tensor + recurrence_inh_self_column_epsc_tensor ...
                    + recurrence_exc_neighbour_column_epsc_tensor + recurrence_inh_neighbour_column_epsc_tensor ;


% fill the spike rates tensor
for i=1:n_iters
    for n=1:n_total_neurons
        spikes1 = spikes(i, 1, n, :);
        spikes1_reshaped = reshape(spikes1, 1,length(tspan));
        spike_rate1 = spikes_to_spike_rate_neat(spikes1_reshaped, physical_time_in_ms, dt, spike_rate_dt);
        spikes_rate1 = reshape(spike_rate1, 1,1,1,spike_rate_length);
        spike_rates(i,1,n,:) = spikes_rate1; 
    end
end

for random_neuron=1:n_total_neurons
    clf
figure(random_neuron)
    hold on
        mean_spike_rate_of_random_neuron = zeros(1, spike_rate_length);
        for t=1:spike_rate_length
              mean_spike_rate_of_random_neuron(1,t) = spike_rates(iter_to_see, col, random_neuron,t);
        end
        

        plot(mean_spike_rate_of_random_neuron*100);

        % adjusting from length(tspan)-1 to length(tspan)
        epsc_all_iters_for_random_neuron = total_input_epsc(:, col, random_neuron, :);
        actual_epsc_input_to_random_neuron_size_adjusted = zeros(n_iters, length(tspan));
        
        for iter=1:n_iters
            actual_epsc_input_to_random_neuron_size_adjusted(iter, 2:length(tspan)) = squeeze(epsc_all_iters_for_random_neuron(iter, 1, 1,:));
        end
        
        % taking avg iters
        random_epsc_iters_avg = zeros(1, length(tspan));
        for t=1:length(tspan)
              random_epsc_iters_avg(1,t) = actual_epsc_input_to_random_neuron_size_adjusted(iter_to_see,t);
        end
        
        epsc_to_random_neuron_binned = spikes_to_spike_rate_neat(random_epsc_iters_avg, physical_time_in_ms, dt, spike_rate_dt);
        plot(epsc_to_random_neuron_binned)
        threshold = squeeze(theta_tensor(iter_to_see, col, random_neuron, :));
        threshold_binned = spikes_to_spike_rate_neat(threshold, physical_time_in_ms, dt, spike_rate_dt);
        plot(threshold_binned);

        v = squeeze(voltages(iter_to_see, col, random_neuron, :));
        v_binned = spikes_to_spike_rate_neat(v, physical_time_in_ms, dt, spike_rate_dt);
        plot(v_binned);

      if col == 2
          plot(lamda_std_protochol)
      elseif col == 4
          plot(lamda_dev_protochol)
      end

    hold off
        title('random neuron epsc and psth')
        legend('psth', 'epsc', 'threshold', 'voltage', 'protochol')
grid
pause
clf
end