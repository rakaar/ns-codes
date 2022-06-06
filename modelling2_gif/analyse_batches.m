clear all;close all;
batches = [1, 5, 10, 100, 105, 110, 190, 195, 200];

for b=1:length(batches)
    close all;
    batch_num = batches(b);
    fprintf("\n batch num %d \n",batch_num);
    
col = 1;
% batch_data_path = "D:\batches_data";
batch_data_path = "D:\4_batches_data";
batches_images_path = "D:\4_batches_images\";
batch_file_name = batch_data_path + "\batch_" + num2str(batch_num) + ".mat";
folder_name = "batch_" + num2str(batch_num);
mkdir(batches_images_path,folder_name)

images_path = batches_images_path + "batch_" + num2str(batch_num) + "\";

spikes = load(batch_file_name,"spikes").spikes;
n_iters = load(batch_file_name,"n_iters").n_iters;
n_total_neurons = load(batch_file_name,"n_total_neurons").n_total_neurons;
tspan = load(batch_file_name,"tspan").tspan;
dt = load(batch_file_name,"dt").dt;

thalamic_epsc_tensor = load(batch_file_name,"thalamic_epsc_tensor").thalamic_epsc_tensor;
recurrence_exc_self_column_epsc_tensor = load(batch_file_name,"recurrence_exc_self_column_epsc_tensor").recurrence_exc_self_column_epsc_tensor;
recurrence_inh_self_column_epsc_tensor = load(batch_file_name,"recurrence_inh_self_column_epsc_tensor").recurrence_inh_self_column_epsc_tensor;
recurrence_exc_neighbour_column_epsc_tensor = load(batch_file_name,"recurrence_exc_neighbour_column_epsc_tensor").recurrence_exc_neighbour_column_epsc_tensor;
recurrence_inh_neighbour_column_epsc_tensor = load(batch_file_name,"recurrence_inh_neighbour_column_epsc_tensor").recurrence_inh_neighbour_column_epsc_tensor;
spike_rates = load(batch_file_name,"spike_rates").spike_rates;

protochol = load(batch_file_name,"protochol").protochol;
physical_time_in_ms = load(batch_file_name,"physical_time_in_ms").physical_time_in_ms;
spike_rate_dt = load(batch_file_name,"spike_rate_dt").spike_rate_dt;
spike_rate_length = load(batch_file_name,"spike_rate_length").spike_rate_length;
n_excitatory = load(batch_file_name,"n_excitatory").n_excitatory;
n_inhibitory = load(batch_file_name,"n_inhibitory").n_inhibitory;
n_columns = load(batch_file_name,"n_columns").n_columns;

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
image_name = images_path + "total-input-epsc.fig"; 
saveas(gcf, image_name);
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

image_name = images_path + "psth-type-wise.fig"; 
saveas(gcf, image_name);        
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
    
    image_name = images_path + "exc-and-inh-epsc.fig"; 
    saveas(gcf, image_name);
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

image_name = images_path + "recurrence.fig"; 
saveas(gcf, image_name);        
grid


% --- thalamic epsc to l4 --
[mean_thalamic_epsc_for_iters, mean_thalamic_epsc_for_neurons] = get_mean(thalamic_epsc_tensor, n_iters, n_total_neurons, length(tspan)-1, 1);
figure
    plot(mean_thalamic_epsc_for_neurons)
    title('thalamic epsc to l4 all')
    image_name = images_path + "thalamic-epsc.fig"; 
    saveas(gcf, image_name);
grid

figure
    [mean_spike_rate_for_iters, mean_spike_rate_for_neurons] = get_mean(spike_rates, n_iters, n_total_neurons, spike_rate_length,1);
    n_bins = spike_rate_dt/dt;
    hold on
        mean_input_epsc_extended = zeros(1, length(tspan));
        mean_input_epsc_extended(1, 2:length(tspan)) = mean_input_epsc_all_for_neurons;
        mean_input_epsc_binned = spikes_to_spike_rate_neat(mean_input_epsc_extended, 1, dt, spike_rate_dt);
        plot(mean_input_epsc_binned)
        plot(mean_spike_rate_for_neurons/0.001)
        
        plot(protochol);
    hold off
    title('total input and psth')
    image_name = images_path + "total-input-and-psth.fig"; 
    saveas(gcf, image_name);
    legend('total input epsc', 'psth l4', 'protochol')
grid


% 5 bins - psth and epsc
spike_rate_dt_5 = 5*dt;
spike_rate_length_5 = (length(tspan)-1)/(spike_rate_dt_5/dt);
spike_rates_5 = zeros(n_iters, n_columns, n_total_neurons, spike_rate_length_5);

for i=1:n_iters
    for n=1:n_total_neurons
        spikes1 = spikes(i, 1, n, :);
        spikes1_reshaped = reshape(spikes1, 1,length(tspan));
        spike_rate1 = spikes_to_spike_rate_neat(spikes1_reshaped, physical_time_in_ms, dt, spike_rate_dt_5);
        spikes_rate1 = reshape(spike_rate1, 1,1,1,spike_rate_length_5);
        spike_rates_5(i,1,n,:) = spikes_rate1;
    end
end

figure
    [mean_spike_rate_for_iters, mean_spike_rate_for_neurons] = get_mean(spike_rates_5, n_iters, n_total_neurons, spike_rate_length_5,1);
    n_bins = spike_rate_dt_5/dt;
    hold on
        mean_input_epsc_extended = zeros(1, length(tspan));
        mean_input_epsc_extended(1, 2:length(tspan)) = mean_input_epsc_all_for_neurons;
        mean_input_epsc_binned = spikes_to_spike_rate_neat(mean_input_epsc_extended, 1, dt, spike_rate_dt_5);
        plot(mean_input_epsc_binned)
        plot(mean_spike_rate_for_neurons/(n_bins*0.001))
        
        protochol_binned = spikes_to_spike_rate_neat(protochol, 1, dt, spike_rate_dt_5);
        plot(protochol_binned);
    hold off
    title('5 bins total input and psth')
    image_name = images_path + "5-bins.fig"; 
    saveas(gcf, image_name);
    legend('total input epsc', 'psth l4','protochol')
grid


for iter=1:n_iters
    figure(iter*100 + 77)
        hold on
            c = 1;
            spike_reshaped = reshape(spikes(iter,c,:,:),  n_total_neurons, length(tspan));
            imagesc(spike_reshaped);
            title('single iter raster l4')
            plot(protochol/10)
        hold off
        image_name = images_path + "raster.fig"; 
        saveas(gcf, image_name);
    grid
end




end % end of for all batches

