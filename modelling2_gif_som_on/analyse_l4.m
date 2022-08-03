% stimulus 

figure
    plot(lamda_std_protochol)
    title('std protocol')
grid

figure
    plot(lamda_dev_protochol)
    title('dev protocol')
grid


% psth of each columns

spike_rates = zeros(n_iters, n_columns, n_total_neurons, spike_rate_length);

for i=1:n_iters
    for c=1:n_columns
        for n=1:n_total_neurons
            spikes1 = spikes(i,c,n,:);
            spikes1_reshaped = reshape(spikes1, 1,length(tspan));
            spike_rate1 = spikes_to_spike_rate_neat(spikes1_reshaped, physical_time_in_ms, dt, spike_rate_dt);
            spikes_rate1 = reshape(spike_rate1, 1,1,1,spike_rate_length);
            spike_rates(i,c,n,:) = spikes_rate1; 
        end
    end
end

% get a mean of all spikes
spike_rate_l4 = zeros(n_columns, n_total_neurons, spike_rate_length);
for i=1:spike_rate_length
   for c=1:n_columns
     for n=1:n_total_neurons
        spike_rate_l4(c,n, i) = sum(spike_rates(:,c,n,i))/n_iters;
    end
   end
end

% psth of each columns
for c=1:n_columns
    figure(c)
        spike_rate_per_column = zeros(1, spike_rate_length);
        for t=1:spike_rate_length
            spike_rate_per_column(1,t) = sum(spike_rate_l4(c,:,t))/n_total_neurons;
        end

        plot(spike_rate_per_column);
        title('psth of col', num2str(c));
    grid
end

%% stimulus 
% std stimulus
thalamic_poission_std_avg = zeros(n_thalamic, length(tspan));
for t=1:length(tspan)
    for n=1:n_thalamic
        thalamic_poission_std_avg(n,t) = sum(thalamic_poisson_spikes_std(:,n,t))/n_iters;
    end
end
thalamic_poisson_std_all_neurons_avg = zeros(1, length(tspan));
for t=1:length(tspan)
    thalamic_poisson_std_all_neurons_avg(1,t) = sum(thalamic_poission_std_avg(:,t))/n_thalamic;
end

% dev stimulus
thalamic_poission_dev_avg = zeros(n_thalamic, length(tspan));
for t=1:length(tspan)
    for n=1:n_thalamic
        thalamic_poission_dev_avg(n,t) = sum(thalamic_poisson_spikes_dev(:,n,t))/n_iters;
    end
end
thalamic_poisson_dev_all_neurons_avg = zeros(1, length(tspan));
for t=1:length(tspan)
    thalamic_poisson_dev_all_neurons_avg(1,t) = sum(thalamic_poission_dev_avg(:,t))/n_thalamic;
end

% col 2 and 4 - epsc and psth
total_input_epsc = thalamic_epsc_tensor ...
                    + recurrence_exc_self_column_epsc_tensor + recurrence_inh_self_column_epsc_tensor ...
                    + recurrence_exc_neighbour_column_epsc_tensor + recurrence_inh_neighbour_column_epsc_tensor ;

[mean_input_epsc_all_for_iters_col_2, mean_input_epsc_all_for_neurons_col_2] = get_mean(total_input_epsc, n_iters, n_total_neurons, length(tspan)-1, 2);
[mean_input_epsc_all_for_iters_col_4, mean_input_epsc_all_for_neurons_col_4] = get_mean(total_input_epsc, n_iters, n_total_neurons, length(tspan)-1, 4);

%% - col 2 and 4 - epsc and psth
visual_comfort_scale_factor = 30;

figure
    hold on
    
    plot(mean_input_epsc_all_for_neurons_col_2/visual_comfort_scale_factor);

    c = 2;
    spike_rate_per_column_std = zeros(1, spike_rate_length);
     for t=1:spike_rate_length
        spike_rate_per_column_std(1,t) = sum(spike_rate_l4(c,:,t))/n_total_neurons;
      end

        plot(spike_rate_per_column_std);
        plot(lamda_std_protochol./100);
        title('col 2 psth with input/visualcomfortscalefactor')
        legend('stim','psth','proto')
grid

% dev
figure
    hold on
     
    plot(mean_input_epsc_all_for_neurons_col_4/visual_comfort_scale_factor);
    c = 4;
    spike_rate_per_column = zeros(1, spike_rate_length);
     for t=1:spike_rate_length
        spike_rate_per_column(1,t) = sum(spike_rate_l4(c,:,t))/n_total_neurons;
      end

        plot(spike_rate_per_column);
        plot(lamda_dev_protochol./100);
        title('col 4 psth with input')
        legend('stim','psth','proto')
%% end of col 2 and 4 - psth and epsc

%%  raster plot
for c=1:n_columns
    spikes1 = spikes(:,c,:,:);
    spikes1 = squeeze(spikes1);
    figure(c*10)
        hold on
            raster = reshape(spikes1, n_iters*n_total_neurons, length(tspan));
            imagesc(raster);
            if c==2
                plot(lamda_std_protochol./10, 'g')
            elseif c==4
                plot(lamda_dev_protochol./10, 'g')
            end
            title('raster of col,', num2str(c));
            
        hold off
    grid
end

%% time constant by fitting exponential
% num_bins = spike_rate_dt/dt;
% t_offset = (pre_stimulus_time/num_bins);
% psth_during_stim = spike_rate_per_column_std(1, pre_stimulus_time/num_bins:(pre_stimulus_time + (single_stimulus_duration+gap_duration)*n_tokens)/num_bins);
% period_stimulus = (pre_stimulus_time/num_bins)-t_offset:((pre_stimulus_time + (single_stimulus_duration+gap_duration)*n_tokens)/num_bins) - t_offset;
% f = fit(period_stimulus',psth_during_stim','exp1');
% coeffs = coeffvalues(f);
% time_const = 1/coeffs(2)
% figure
%   hold on
%     plot(psth_during_stim)
%     plot(f)
%   hold off
% grid

%% inhibitory and excitatory of std and dev
% -- std --
c = 2;
    spike_rate_per_column_std_exc = zeros(1, spike_rate_length);
    spike_rate_per_column_std_inh = zeros(1, spike_rate_length);
     for t=1:spike_rate_length
        spike_rate_per_column_std_exc(1,t) = sum(spike_rate_l4(c,1:n_excitatory,t))/n_excitatory;
     end
     for t=1:spike_rate_length
        spike_rate_per_column_std_inh(1,t) = sum(spike_rate_l4(c,n_excitatory+1:n_total_neurons,t))/n_inhibitory;
      end
figure
hold on
        plot(spike_rate_per_column_std_exc);
        plot(spike_rate_per_column_std_inh)
        title('exc inh std col -2')
        legend('exc', 'inh')
hold off
grid
% --- dev ---
c = 4;
    spike_rate_per_column_dev_exc = zeros(1, spike_rate_length);
    spike_rate_per_column_dev_inh = zeros(1, spike_rate_length);
     for t=1:spike_rate_length
        spike_rate_per_column_dev_exc(1,t) = sum(spike_rate_l4(c,1:n_excitatory,t))/n_excitatory;
     end
     for t=1:spike_rate_length
        spike_rate_per_column_dev_inh(1,t) = sum(spike_rate_l4(c,n_excitatory+1:n_total_neurons,t))/n_inhibitory;
      end
figure
hold on
        plot(spike_rate_per_column_dev_exc);
        plot(spike_rate_per_column_dev_inh)
        title('exc inh std col -4')
        legend('exc', 'inh')
hold off
grid

% epsc - exc and inhi
recurrence_exc_epsc = recurrence_exc_self_column_epsc_tensor + recurrence_exc_neighbour_column_epsc_tensor;
recurrence_inh_epsc = recurrence_inh_self_column_epsc_tensor + recurrence_inh_neighbour_column_epsc_tensor;
recurrence_epsc = recurrence_exc_epsc + recurrence_inh_epsc;

figure
    hold on
        [mean_epsc_exc_recurrence_for_iters, mean_epsc_exc_reccurence_for_neurons] = get_mean(recurrence_exc_epsc, n_iters, n_total_neurons, length(tspan)-1, 2);
        [mean_epsc_inh_recurrence_for_iters, mean_epsc_inh_reccurence_for_neurons] = get_mean(recurrence_inh_epsc, n_iters, n_total_neurons, length(tspan)-1, 2);
        
        plot(mean_epsc_exc_reccurence_for_neurons)
        plot(mean_epsc_inh_reccurence_for_neurons)
        legend('exc epsc to all l4', 'inh epsc to all l4')
        title('epscs-exc and inh to col 2')
    hold off
grid

figure
    hold on
        [mean_epsc_exc_recurrence_for_iters, mean_epsc_exc_reccurence_for_neurons] = get_mean(recurrence_exc_epsc, n_iters, n_total_neurons, length(tspan)-1, 4);
        [mean_epsc_inh_recurrence_for_iters, mean_epsc_inh_reccurence_for_neurons] = get_mean(recurrence_inh_epsc, n_iters, n_total_neurons, length(tspan)-1, 4);
        
        plot(mean_epsc_exc_reccurence_for_neurons)
        plot(mean_epsc_inh_reccurence_for_neurons)
        legend('exc epsc to all l4', 'inh epsc to all l4')
        title('epscs-exc and inh to col 4')
    hold off
grid