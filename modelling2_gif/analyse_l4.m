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

figure
    hold on
    psth_stim = spikes_to_spike_rate_neat(thalamic_poisson_std_all_neurons_avg,physical_time_in_ms, dt, spike_rate_dt );
    plot(psth_stim) 
    
    c = 2;
    spike_rate_per_column = zeros(1, spike_rate_length);
     for t=1:spike_rate_length
        spike_rate_per_column(1,t) = sum(spike_rate_l4(c,:,t))/n_total_neurons;
      end

        plot(spike_rate_per_column);

%% dev
figure
    hold on
    psth_stim = spikes_to_spike_rate_neat(thalamic_poisson_dev_all_neurons_avg,physical_time_in_ms, dt, spike_rate_dt );
    plot(psth_stim) 
    
    c = 4;
    spike_rate_per_column = zeros(1, spike_rate_length);
     for t=1:spike_rate_length
        spike_rate_per_column(1,t) = sum(spike_rate_l4(c,:,t))/n_total_neurons;
      end

        plot(spike_rate_per_column);


%%  raster plot
for c=1:n_columns
    spikes1 = spikes(:,c,:,:);
    spikes1 = squeeze(spikes1);
    figure(c*10)
        raster = reshape(spikes1, n_iters*n_total_neurons, length(tspan));
        imagesc(raster);
        title('raster of col', num2str(c));
    grid
end