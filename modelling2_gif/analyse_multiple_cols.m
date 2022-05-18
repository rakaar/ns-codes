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

% -- plot total input epsc to l4 (leave I_background)
total_input_epsc = thalamic_epsc_tensor ...
                    + recurrence_exc_self_column_epsc_tensor + recurrence_inh_self_column_epsc_tensor ...
                    + recurrence_exc_neighbour_column_epsc_tensor + recurrence_inh_neighbour_column_epsc_tensor ;

[mean_input_epsc_exc_for_iters, mean_input_epsc_exc_for_neurons] = get_mean(total_input_epsc(:,:,1:n_excitatory,:), n_iters, n_excitatory, length(tspan)-1, 1);
[mean_input_epsc_inh_for_iters, mean_input_epsc_inh_for_neurons] = get_mean(total_input_epsc(:,:,n_excitatory+1:n_total_neurons,:), n_iters, n_inhibitory, length(tspan)-1, 1);
[mean_input_epsc_all_for_iters, mean_input_epsc_all_for_neurons] = get_mean(total_input_epsc, n_iters, n_total_neurons, length(tspan)-1, 1);
figure
    subplot(1,2,1)     
    hold on
    
        plot(mean_input_epsc_exc_for_neurons)
             plot(mean_input_epsc_inh_for_neurons)
             plot(mean_input_epsc_all_for_neurons)
             legend('total input to exc l4', 'total input to inh l4', 'total input to l4 all')
            title('l4 - total input epsc')
        hold off
    
    subplot(1,2,2)
    plot(mean_input_epsc_all_for_neurons)
    title('l4 total input epsc')

grid

% -- plot psth of l4 exc, inh, all
figure
    subplot(1,2,1)
    hold on
        [mean_spike_rate_exc_for_iters, mean_spike_rate_exc_for_neurons] = get_mean(spike_rates(:,:,1:n_excitatory,:), n_iters, n_excitatory, spike_rate_length,1);
        [mean_spike_rate_inh_for_iters, mean_spike_rate_inh_for_neurons] = get_mean(spike_rates(:,:,n_excitatory+1:n_total_neurons,:), n_iters, n_inhibitory, spike_rate_length,1);
        [mean_spike_rate_for_iters, mean_spike_rate_for_neurons] = get_mean(spike_rates, n_iters, n_total_neurons, spike_rate_length,1);
       
       n_bins = spike_rate_dt/dt;
       
       plot(mean_spike_rate_exc_for_neurons*(n_bins*physical_time_in_ms*0.001))
       plot(mean_spike_rate_inh_for_neurons*(n_bins*physical_time_in_ms*0.001))
       plot(mean_spike_rate_for_neurons*(n_bins*physical_time_in_ms*0.001));
       legend('psth l4 exc', 'psth l4 inh','psth l4 all')
    hold off

    subplot(1,2,2)
    plot(mean_spike_rate_for_neurons*(n_bins*physical_time_in_ms*0.001))
    title('psth l4 all')
grid

recurrence_exc_epsc = recurrence_exc_self_column_epsc_tensor + recurrence_exc_neighbour_column_epsc_tensor;
recurrence_inh_epsc = recurrence_inh_self_column_epsc_tensor + recurrence_inh_neighbour_column_epsc_tensor;
recurrence_epsc = recurrence_exc_epsc + recurrence_inh_epsc;

figure
    hold on
        [mean_epsc_exc_recurrence_for_iters, mean_epsc_exc_reccurence_for_neurons] = get_mean(recurrence_exc_epsc, n_iters, n_total_neurons, length(tspan)-1, 1);
        [mean_epsc_inh_recurrence_for_iters, mean_epsc_inh_reccurence_for_neurons] = get_mean(recurrence_inh_epsc, n_iters, n_total_neurons, length(tspan)-1, 1);
        
        plot(mean_epsc_exc_reccurence_for_neurons)
        plot(mean_epsc_inh_reccurence_for_neurons)
        legend('exc epsc to all l4', 'inh epsc to all l4')
        title('epscs-exc and inh')
    hold off
grid
figure
   subplot(1,2,1)
    hold on
        [mean_recurrence_epsc_exc_for_iters, mean_recurrence_epsc_exc_for_neurons] = get_mean(recurrence_epsc(:,:,1:n_excitatory,:), n_iters, n_excitatory, length(tspan)-1,1);
        [mean_recurrence_epsc_inh_for_iters, mean_recurrence_epsc_inh_for_neurons] = get_mean(recurrence_epsc(:,:,n_excitatory+1:n_total_neurons,:), n_iters, n_inhibitory, length(tspan)-1,1);
        [mean_recurrence_epsc_all_for_iters, mean_recurrence_epsc_all_for_neurons] = get_mean(recurrence_epsc, n_iters, n_total_neurons, length(tspan)-1, 1);
    
        plot(mean_recurrence_epsc_exc_for_neurons)
        plot(mean_recurrence_epsc_inh_for_neurons)
        plot(mean_recurrence_epsc_all_for_neurons)
        legend('recurrence epsc input to exc in l4', 'recurrence epsc to inh in l4', 'recurrence epsc to all in l4')
        title('recurrrence epscs')
    hold off

    subplot(1,2,2)
    plot(mean_recurrence_epsc_all_for_neurons)
    title('mean recurrence epsc to all l4 neurons')


grid

% --- thalamic epsc to l4 --
[mean_thalamic_epsc_for_iters, mean_thalamic_epsc_for_neurons] = get_mean(thalamic_epsc_tensor, n_iters, n_total_neurons, length(tspan)-1, 1);
figure
    plot(mean_thalamic_epsc_for_neurons)
    title('thalamic epsc to l4 all')
grid

figure
     [mean_spike_rate_for_iters, mean_spike_rate_for_neurons] = get_mean(spike_rates, n_iters, n_total_neurons, spike_rate_length,1);
        n_bins = spike_rate_dt/dt;
    hold on
        mean_input_epsc_extended = zeros(1, length(tspan));
        mean_input_epsc_extended(1, 2:length(tspan)) = mean_input_epsc_all_for_neurons;
        mean_input_epsc_binned = spikes_to_spike_rate_neat(mean_input_epsc_extended, 1, dt, spike_rate_dt);
        plot(mean_input_epsc_binned*(n_bins*physical_time_in_ms*0.001))
        plot(mean_spike_rate_for_neurons)
    hold off
    title('total input and psth')
    legend('total input epsc', 'psth l4')
grid

% for a random neuron
random_neuron = 23;
col = 1;
figure
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

% raster of l4
iter_to_see = 1;
for c=1:n_columns
    figure(c*100 + 77)
        hold on
            imagesc(squeeze(spikes(iter_to_see,c,:,:)));
            if c == 2
                plot(lamda_std_protochol/6)
            elseif c == 4
                plot(lamda_dev_protochol*20)
            else
                plot(lamda_common_protochol)
            end
        hold off
        
        title('raster of l4')
    grid

end
