% to be run after distr del t del w
%% check rates of pre and post neurons
close all;
spike_rate_pre_neuron = zeros(batches,1);
spike_rate_post_neuron = zeros(batches,1);

for b=1:batches
    fprintf("\n batch is %d \n",b)
    batch_file_name = data_path + "\batch_" + num2str(b) + ".mat";
    spikes = load(batch_file_name, "spikes").spikes;
    spike_rate_pre_neuron(b,1) = sum(spikes(iter,pre_col,pre_neuron,:));
    spike_rate_post_neuron(b,1) = sum(spikes(iter,post_col,post_neuron,:));
end

%% pre and post
figure
    hold on
        plot(spike_rate_pre_neuron)
        plot(spike_rate_post_neuron)
    hold off
    title(['spikes-','-pre-col-',num2str(pre_col),'-post-col-',num2str(post_col),'-pre-neuron-',num2str(pre_neuron),'-post-neuron-',num2str(post_neuron)])
    legend('rate pre syn', 'rate post syn')
grid

%% rates pre and post except target
spike_rate_pre_col_except_target = zeros(batches,n_exc-1);
spike_rate_post_col_except_target = zeros(batches,n_exc-1);

for b=1:batches
    fprintf("\n batch is %d \n",b)
    batch_file_name = data_path + "\batch_" + num2str(b) + ".mat";
    spikes = load(batch_file_name, "spikes").spikes;
    
    % pre col
    n_index = 1;
    for n=1:n_exc
        if n ~= pre_neuron
            spike_rate_pre_col_except_target(b,n_index) = sum(spikes(iter,pre_col,n,:));
            n_index = n_index + 1; 
        end
    end

    % post col
    n_index = 1;
    for n=1:n_exc
        if n ~= post_neuron
            spike_rate_post_col_except_target(b,n_index) = sum(spikes(iter,post_col,n,:));
            n_index = n_index + 1; 
        end
    end
end

%% pre col and post col
figure
    hold on
        plot(spike_rate_pre_col_except_target)
        plot(spike_rate_pre_neuron, 'LineWidth',5)
    hold off

    title(['rates-pre-col-',num2str(pre_col),'--target-neuron--', num2str(pre_neuron)])
grid

figure
    hold on
        plot(spike_rate_post_col_except_target)
        plot(spike_rate_post_neuron, 'LineWidth',5)
    hold off

    title(['rates-post-col-',num2str(post_col),'--target-neuron--', num2str(post_neuron)])
grid

%% within column weights and post column weights
weights_within_column_to_post_neuron_all = zeros(batches, n_exc);
weights_across_col_to_post_neuron_all = zeros(batches,n_exc);
target_synapse_weight_over_time = zeros(batches,1);

across_col_index = post_col - 1;

for b=1:batches
    fprintf("\n batch is %d \n",b);
    batch_file_name = data_path + "\batch_" + num2str(b) + ".mat";
    network_weight_matrix = load(batch_file_name, "network_weight_matrix").network_weight_matrix;

    weights_within_column_to_post_neuron_all(b,:) = mean(squeeze(network_weight_matrix(iter,:,(post_col-1)*n_total_neurons+1:(post_col-1)*n_total_neurons+n_exc,n_total_neurons*(post_col-1) + post_neuron)), 1);
    weights_across_col_to_post_neuron_all(b,:) = mean(squeeze(network_weight_matrix(iter,:,(across_col_index-1)*n_total_neurons+1:(across_col_index-1)*n_total_neurons+n_exc,n_total_neurons*(post_col-1) + post_neuron)), 1);

    target_synapse_weight_over_time(b,1) = mean(squeeze(network_weight_matrix(iter,:,(pre_col-1)*n_total_neurons + pre_neuron, (post_col-1)*n_total_neurons + post_neuron)));
 
end

%% weights within and across and target
figure
    hold on
        plot(weights_within_column_to_post_neuron_all)
        plot(target_synapse_weight_over_time, 'LineWidth',5)
    hold off
    title('within post col')
grid


figure
    hold on
        plot(weights_across_col_to_post_neuron_all)
        plot(target_synapse_weight_over_time, 'LineWidth',5)
    hold off
    title('across col to post col')
grid