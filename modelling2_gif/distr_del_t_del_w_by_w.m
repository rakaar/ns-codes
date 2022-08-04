%% imagesc of weights
batches = 100;
data_path = "D:\new-final-working-baseline-params-som-off";
images_path = "D:\new-final-working-baseline-params-som-off-analysis\";
n_excitatory=20; n_pv = 3; n_som  = 2;
n_neurons = n_excitatory + n_pv + n_som;
n_total_neurons = load(strcat(data_path, '\', 'batch_1.mat'), "n_total_neurons").n_total_neurons;
tspan = load(strcat(data_path, '\', 'batch_1.mat'), "tspan").tspan;
num_network_neurons = 125;
batch_avg_network_weights = zeros(batches, num_network_neurons, num_network_neurons);
n_exc = 20;
iter=1;

for b=1:batches
    fprintf("\n batch is %d \n",b)
    batch_file_name = data_path + "\batch_" + num2str(b) + ".mat";
    network_weight_matrix = load(batch_file_name, "network_weight_matrix").network_weight_matrix;
    for n1=1:num_network_neurons
        for n2=1:num_network_neurons
             batch_avg_network_weights(b,n1,n2) = mean(squeeze(network_weight_matrix(iter,:,n1,n2)));
        end
    end
end

n_columns = 5;
for c1=1:n_columns
    for c2=1:n_columns
        if c1 - c2 == 0
            % within column matrix
            within_column_matrix = batch_avg_network_weights(:,(c1-1)*n_total_neurons + 1:(c1-1)*n_total_neurons + n_exc, (c1-1)*n_total_neurons + 1:(c1-1)*n_total_neurons + n_exc);
            fprintf("\n %d %d %d %d \n ",(c1-1)*n_total_neurons + 1,(c1-1)*n_total_neurons + n_exc, (c1-1)*n_total_neurons + 1,(c1-1)*n_total_neurons + n_exc)
            reshaped_within_column_matrix = reshape(within_column_matrix, batches, n_exc*n_exc);
            figure
                imagesc(transpose(reshaped_within_column_matrix))
                title(['weights batches wise', num2str(c1)])
            grid
        elseif abs(c1 - c2) == 1 || abs(c1 - c2) == 2
            % c1 pre ,c2 post
            across_column_matrix = batch_avg_network_weights(:,(c1-1)*n_total_neurons + 1:(c1-1)*n_total_neurons + n_exc, (c2-1)*n_total_neurons + 1:(c2-1)*n_total_neurons + n_exc);
            fprintf("\n %d %d %d %d \n ",(c1-1)*n_total_neurons + 1,(c1-1)*n_total_neurons + n_exc, (c2-1)*n_total_neurons + 1,(c2-1)*n_total_neurons + n_exc)
            reshaped_across_column_matrix = reshape(across_column_matrix, batches, n_exc*n_exc);
            figure
                imagesc(transpose(reshaped_across_column_matrix))
                title(['weights batches wise', num2str(c1), '-',num2str(c2)])
            grid

        else
            continue
        end
    end
end

%% too see pre post weight over batches
close all
pre_col = 2;
post_col = 4;
pre_neuron = 4;
post_neuron = 1;

pre_index_in_matrix = (pre_col-1)*n_total_neurons + pre_neuron;
post_index_in_matrix = (post_col-1)*n_total_neurons + post_neuron;

figure
    weight_over_time = squeeze(batch_avg_network_weights(:, pre_index_in_matrix, post_index_in_matrix));
    plot(weight_over_time)
grid

%% delta t, delta w/w distribution
Amp_strength = 0.015; Amp_weak = 0.021;
tau_strength = 30; tau_weak = 50;

delta_t = [];
delta_w_by_w = [];

iter = 1;
pre_col = 2;
post_col = 4;
pre_neuron = 4;
post_neuron = 1;

for b=1:batches
    fprintf("\n batch num %d \n",b);
    batch_file_name = data_path + "\batch_" + num2str(b) + ".mat";
    spikes = load(batch_file_name,"spikes").spikes;

    for i=6:length(tspan)
         % post - pre = +ve
        if spikes(iter,post_col,post_neuron,i) == 1
            for presyn_time=i-1:-1:i-20
                if presyn_time >= 1 && spikes(iter,pre_col,pre_neuron,presyn_time) == 1 && spikes(iter,pre_col,pre_neuron,i) == 0
                    time_diff = i-presyn_time;
                    delta_t = [delta_t, time_diff];
                    delta_w_by_w = [delta_w_by_w, Amp_strength*exp(-abs(time_diff)/tau_strength)];
                    break;
                end
            end
        end

        % post - pre = -ve
        if spikes(iter,pre_col,pre_neuron,i) == 1
            for postsyn_time=i-1:-1:i-20
                if postsyn_time >= 1 && spikes(iter,post_col,post_neuron,postsyn_time) == 1 && spikes(iter,post_col,post_neuron,i) == 0
                    time_diff = postsyn_time - i;
                    delta_t = [delta_t, time_diff];
                    delta_w_by_w = [delta_w_by_w, -Amp_weak*exp(-abs(time_diff)/tau_weak)];
                    break;
                end
            end
        end

    end
end

%% distr plots
figure
    histogram(delta_t)
    title(['delta-t-','-pre-col-',num2str(pre_col),'-post-col-',num2str(post_col),'-pre-neuron-',num2str(pre_neuron),'-post-neuron-',num2str(post_neuron)])
grid

figure
   histogram(delta_w_by_w)
   title(['delta-w-by-w-','-pre-col-',num2str(pre_col),'-post-col-',num2str(post_col),'-pre-neuron-',num2str(pre_neuron),'-post-neuron-',num2str(post_neuron)])
grid