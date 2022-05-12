close all;
% for a random neuron
col = 1;
iter_to_see=1;
n_bins = spike_rate_dt/dt;
multiply_term = (n_bins*physical_time_in_ms*0.001);
for random_neuron=1:n_total_neurons
    clf
figure(random_neuron)
hold on
        % spike rate
        mean_spike_rate_of_random_neuron = zeros(1, spike_rate_length);
        for t=1:spike_rate_length
              mean_spike_rate_of_random_neuron(1,t) = spike_rates(iter_to_see, col, random_neuron,t);
        end
        

        plot(mean_spike_rate_of_random_neuron*100,'LineStyle','--');

%         % adjusting from length(tspan)-1 to length(tspan)
%         epsc_all_iters_for_random_neuron = squeeze(total_input_epsc(:, col, random_neuron, :));
%         actual_epsc_input_to_random_neuron_size_adjusted = zeros(n_iters, length(tspan));
%         for iter=1:n_iters
%             actual_epsc_input_to_random_neuron_size_adjusted(iter, 2:length(tspan)) = epsc_all_iters_for_random_neuron(iter, :);
%         end
%         
%         % taking avg iters
%         random_epsc_iters_avg = zeros(1, length(tspan));
%         for t=1:length(tspan)
%               random_epsc_iters_avg(1,t) = actual_epsc_input_to_random_neuron_size_adjusted(iter_to_see,t);
%         end
%         
%         epsc_to_random_neuron_binned = spikes_to_spike_rate_neat(random_epsc_iters_avg, physical_time_in_ms, dt, spike_rate_dt);
%         plot(epsc_to_random_neuron_binned,'m')



        % threshold
        threshold = squeeze(theta_tensor(iter_to_see, col, random_neuron, :));
        threshold_binned = spikes_to_spike_rate_neat(threshold, physical_time_in_ms, dt, spike_rate_dt);
        plot(threshold_binned, 'c');
        
        % voltage
        v = squeeze(voltages(iter_to_see, col, random_neuron, :));
        v_binned = spikes_to_spike_rate_neat(v, physical_time_in_ms, dt, spike_rate_dt);
        plot(v_binned,'r');


        %  thalamic current
        thalamic_curr_to_neuron = thalamic_epsc_tensor(iter_to_see, col, random_neuron, :);
        thalamic_curr_to_neuron = squeeze(thalamic_curr_to_neuron);
        thalamic_curr_to_neuron = [0; thalamic_curr_to_neuron];
        plot(thalamic_curr_to_neuron,'g');
    
        %  exc epsc
        exc_epsc = recurrence_exc_self_column_epsc_tensor(iter_to_see, col, random_neuron, :) + recurrence_exc_neighbour_column_epsc_tensor(iter_to_see, col, random_neuron, :);
        exc_epsc = squeeze(exc_epsc);
        exc_epsc = [0; exc_epsc];
        plot(exc_epsc,'b');

        %  inh epsc
        inh_epsc = recurrence_inh_self_column_epsc_tensor(iter_to_see, col, random_neuron, :) + recurrence_inh_neighbour_column_epsc_tensor(iter_to_see, col, random_neuron, :);
        inh_epsc = squeeze(inh_epsc);
        inh_epsc = [0; inh_epsc];
        plot(inh_epsc,'k');

        %  recurrence epsc
        rec_epsc = exc_epsc + inh_epsc;
        plot(rec_epsc,'y');

        % all epsc
        all_epsc = rec_epsc + thalamic_curr_to_neuron;
        plot(all_epsc, 'm');

        % background
        i_back = I_background_tensor(iter_to_see, col, random_neuron, :);
        i_back = squeeze(i_back);
        i_back = [0; i_back];
        plot(i_back, '--p');

    hold off
        title('random neuron epsc and psth')
        legend('spike rate','threshold','voltage','thalamic epsc','exc epsc','inh epsc','rec epsc','all epsc')
grid
pause
clf
end

