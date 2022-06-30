clear all;close all;

batches = [1, 5, 10, 20, 25, 30];


for b=1:length(batches)
    close all;
    batch_num = batches(b);
    fprintf("\n batch num %d \n",batch_num);
    
col = 1;
batch_data_path = "D:\2_multi_col_across_plastic";
batches_images_path = "D:\2_multi_col_images\";
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

lamda = load(batch_file_name,"lamda").lamda;
physical_time_in_ms = load(batch_file_name,"physical_time_in_ms").physical_time_in_ms;
spike_rate_dt = load(batch_file_name,"spike_rate_dt").spike_rate_dt;
spike_rate_length = load(batch_file_name,"spike_rate_length").spike_rate_length;
n_excitatory = load(batch_file_name,"n_excitatory").n_excitatory;
n_inhibitory = load(batch_file_name,"n_inhibitory").n_inhibitory;
n_columns = load(batch_file_name,"n_columns").n_columns;

% fill the spike rates tensor
for i=1:n_iters
    for c=1:n_columns
        for n=1:n_total_neurons
            spikes1 = spikes(i, c, n, :);
            spikes1_reshaped = reshape(spikes1, 1,length(tspan));
            spike_rate1 = spikes_to_spike_rate_neat(spikes1_reshaped, physical_time_in_ms, dt, spike_rate_dt);
            spikes_rate1 = reshape(spike_rate1, 1,1,1,spike_rate_length);
            spike_rates(i,c,n,:) = spikes_rate1;
        end
   end
end

% -- plot total input epsc to l4 (leave I_background)
total_input_epsc = thalamic_epsc_tensor ...
    + recurrence_exc_self_column_epsc_tensor + recurrence_inh_self_column_epsc_tensor ...
    + recurrence_exc_neighbour_column_epsc_tensor + recurrence_inh_neighbour_column_epsc_tensor ;


iter_to_see = 1; random_thal_neuron = 2;
spike_rate_dt_5 = 5*dt;
spike_rate_length_5 = (length(tspan)-1)/(spike_rate_dt_5/dt);
spike_rates_5 = zeros(n_iters, n_columns, n_total_neurons, spike_rate_length_5);

for i=1:n_iters
    for c=1:n_columns
        for n=1:n_total_neurons
            spikes1 = spikes(i, c, n, :);
            spikes1_reshaped = reshape(spikes1, 1,length(tspan));
            spike_rate1 = spikes_to_spike_rate_neat(spikes1_reshaped, physical_time_in_ms, dt, spike_rate_dt_5);
            spikes_rate1 = reshape(spike_rate1, 1,1,1,spike_rate_length_5);
            spike_rates_5(i,c,n,:) = spikes_rate1;
        end
    end
end

for c=1:n_columns
    [mean_spike_rate_for_iters, mean_spike_rate_for_neurons] = get_mean(spike_rates_5(iter_to_see,c,:,:), n_iters, n_total_neurons, spike_rate_length_5,1);
    [mean_input_epsc_all_for_iters, mean_input_epsc_all_for_neurons] = get_mean(total_input_epsc(iter_to_see,c,:,:), n_iters, n_total_neurons, length(tspan)-1, 1);

    figure
        n_bins = spike_rate_dt_5/dt;
        hold on
            mean_input_epsc_extended = zeros(1, length(tspan));
            mean_input_epsc_extended(1, 2:length(tspan)) = mean_input_epsc_all_for_neurons;
            mean_input_epsc_binned = spikes_to_spike_rate_neat(mean_input_epsc_extended, 1, dt, spike_rate_dt_5);
            plot(mean_input_epsc_binned)
            plot(mean_spike_rate_for_neurons/(n_bins*0.001))
            
            protochol = squeeze(lamda(iter_to_see,c+2,random_thal_neuron,:));
            protochol_binned = spikes_to_spike_rate_neat(protochol, 1, dt, spike_rate_dt_5);
            plot(protochol_binned);

        hold off
        title(['5 bins total input and psth','col ', num2str(c)])
        legend('total input epsc', 'psth l4','protochol/100')
    
         image_name = images_path + "-col-" + num2str(c) + "psth.fig"; 
         saveas(gcf, image_name);
    grid

    
end




iter_to_see = 1;
for col_to_see=1:n_columns
    figure
        hold on
            spikes_col = spikes(iter_to_see, col_to_see, :, :);
            imagesc(squeeze(spikes_col))
            plot(squeeze(lamda(iter_to_see,col_to_see+2,1,:)))
            title(['col ',num2str(col_to_see)])
        hold off
        image_name = images_path + "-col-" + num2str(col_to_see) + "raster.fig"; 
        saveas(gcf, image_name);
    grid
end




end % end of for all batches

