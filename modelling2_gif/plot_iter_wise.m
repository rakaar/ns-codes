% for a random neuron
close all;
col = 4;
iter_to_see=1;
n_bins = spike_rate_dt/dt;

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
        spikes_neuron = squeeze(spikes(iter_to_see, col, random_neuron, :));
        plot(spikes_neuron*100,'LineStyle','--');

        % threshold
        threshold = squeeze(theta_tensor(iter_to_see, col, random_neuron, :));
        threshold_binned = spikes_to_spike_rate_neat(threshold, physical_time_in_ms, dt, spike_rate_dt);
        plot(threshold_binned, 'c');
        
        % voltage
        v = squeeze(voltages(iter_to_see, col, random_neuron, :));
        v_binned = spikes_to_spike_rate_neat(v, physical_time_in_ms, dt, spike_rate_dt);
        plot(v_binned,'r');

        % thalamic current
        thalamic_curr_to_neuron = thalamic_epsc_tensor(iter_to_see, col, random_neuron, :);
        thalamic_curr_to_neuron = squeeze(thalamic_curr_to_neuron);
        thalamic_curr_to_neuron = [0; thalamic_curr_to_neuron];
        plot(thalamic_curr_to_neuron,'g');
    
        % exc epsc own col
        exc_epsc = recurrence_exc_self_column_epsc_tensor(iter_to_see, col, random_neuron, :) + recurrence_exc_neighbour_column_epsc_tensor(iter_to_see, col, random_neuron, :);
        exc_epsc = squeeze(exc_epsc);
        exc_epsc = [0; exc_epsc];
        plot(exc_epsc,'b', 'LineWidth',2);

        % inh epsc own col
        inh_epsc = recurrence_inh_self_column_epsc_tensor(iter_to_see, col, random_neuron, :) + recurrence_inh_neighbour_column_epsc_tensor(iter_to_see, col, random_neuron, :);
        inh_epsc = squeeze(inh_epsc);
        inh_epsc = [0; inh_epsc];
        plot(inh_epsc,'k');

        % epsc from neighbour
        epsc_neigh = recurrence_exc_neighbour_column_epsc_tensor(iter_to_see, col, random_neuron, :);
        epsc_neigh = squeeze(epsc_neigh);
        epsc_neigh = [0; epsc_neigh];
        plot(epsc_neigh,'--s', 'LineWidth',2)

        % recurrence own column
        rec_epsc_own_col = exc_epsc + inh_epsc;
        plot(rec_epsc_own_col,'y','LineWidth',3);
        
         % background
        i_back = I_background_tensor(iter_to_see, col, random_neuron, :);
        i_back = squeeze(i_back);
        i_back = [0; i_back];
        plot(i_back, '--p');

        % all epsc
        all_epsc = rec_epsc_own_col + epsc_neigh + thalamic_curr_to_neuron;
        plot(all_epsc, 'm');

        % protochol
        plot(squeeze(lamda(iter_to_see,col+2,1,:)))
      
    hold off
        title('random neuron epsc*100 and psth')
        legend('spike rate','threshold','voltage','thalamic col','exc own col','inh own col','epsc neighbour','recurrence own column','i background','all epsc','protochol')
grid
pause
clf
end