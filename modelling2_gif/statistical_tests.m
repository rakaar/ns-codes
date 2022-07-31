%% get the weight matrices
close all;
som_on_data_path = "D:\29july-som-on";
som_off_data_path = "D:\29uly-morning";
batches = 100;
n_total_neurons = 25;
num_network_neurons = 125;
n_exc = 20;
iter=1;

% ------------------for som on ---------------------------
data_path = som_on_data_path; 
tspan = load(strcat(data_path, '\', 'batch_1.mat'), "tspan").tspan;
batch_avg_network_weights = zeros(batches, num_network_neurons, num_network_neurons);

for b=1:batches
    fprintf("\n batch is %d \n",b)
    batch_file_name = data_path + "\batch_" + num2str(b) + ".mat";
    network_weight_matrix = load(batch_file_name, "network_weight_matrix").network_weight_matrix;
    for n1=1:num_network_neurons
        for n2=1:num_network_neurons
             batch_avg_network_weights(b,n1,n2) = mean(squeeze(network_weight_matrix(iter,:,n1,n2)));
%              all_1ms_network_weights((b-1)*length(tspan) + 1:(b-1)*length(tspan) + length(tspan)   ,n1,n2) =  network_weight_matrix(iter,:,n1,n2);
        end
    end
end

som_on_batch_avg_network_weights = batch_avg_network_weights; % weights som on

% ------------------ for som off ---------------------------
data_path = som_off_data_path; 
tspan = load(strcat(data_path, '\', 'batch_1.mat'), "tspan").tspan;
batch_avg_network_weights = zeros(batches, num_network_neurons, num_network_neurons);

for b=1:batches
    fprintf("\n batch is %d \n",b)
    batch_file_name = data_path + "\batch_" + num2str(b) + ".mat";
    network_weight_matrix = load(batch_file_name, "network_weight_matrix").network_weight_matrix;
    for n1=1:num_network_neurons
        for n2=1:num_network_neurons
             batch_avg_network_weights(b,n1,n2) = mean(squeeze(network_weight_matrix(iter,:,n1,n2)));
%              all_1ms_network_weights((b-1)*length(tspan) + 1:(b-1)*length(tspan) + length(tspan)   ,n1,n2) =  network_weight_matrix(iter,:,n1,n2);
        end
    end
end

som_off_batch_avg_network_weights = batch_avg_network_weights; % weights som on

%% compare weight matrices
c1 = 4; c2 = 4;
som_on_across_column_matrix = som_on_batch_avg_network_weights(:,(c1-1)*n_total_neurons + 1:(c1-1)*n_total_neurons + n_exc, (c2-1)*n_total_neurons + 1:(c2-1)*n_total_neurons + n_exc);
som_on_reshaped_across_column_matrix = reshape(som_on_across_column_matrix, batches, n_exc*n_exc);


som_off_across_column_matrix = som_off_batch_avg_network_weights(:,(c1-1)*n_total_neurons + 1:(c1-1)*n_total_neurons + n_exc, (c2-1)*n_total_neurons + 1:(c2-1)*n_total_neurons + n_exc);
som_off_reshaped_across_column_matrix = reshape(som_off_across_column_matrix, batches, n_exc*n_exc);

diff_som_on_weights = diff(som_on_reshaped_across_column_matrix, 1, 1);
diff_mean_som_on = mean(diff_som_on_weights, 2);

diff_som_off_weights = diff(som_off_reshaped_across_column_matrix, 1, 1);
diff_mean_som_off = mean(diff_som_off_weights, 2);

[h, p] = ttest(diff_mean_som_on,diff_mean_som_off)

%% rates comparison
n_columns = 5;
n_excitatory = 20;
% ------- for som on -----------
psth_avg_exc = zeros(n_columns,n_excitatory, batches);
data_path = som_on_data_path;
for b=1:batches
    fprintf("\n batch num %d \n",b);

    batch_file_name = data_path + "\batch_" + num2str(b) + ".mat";
    spikes = load(batch_file_name,"spikes").spikes;
    lamda = load(batch_file_name,"lamda").lamda;
    t_simulate = load(batch_file_name,"t_simulate").t_simulate;
    
    for c=1:n_columns
        for nn=1:n_excitatory
            psth_avg_exc(c,nn,b) = sum(spikes(iter,c,nn,:))/(1700 * 0.001);
        end
    end
end
som_on_psth_avg_exc = psth_avg_exc;

% ------- for som off -----------
psth_avg_exc = zeros(n_columns,n_excitatory, batches);
data_path = som_off_data_path;
for b=1:batches
    fprintf("\n batch num %d \n",b);

    batch_file_name = data_path + "\batch_" + num2str(b) + ".mat";
    spikes = load(batch_file_name,"spikes").spikes;
    lamda = load(batch_file_name,"lamda").lamda;
    t_simulate = load(batch_file_name,"t_simulate").t_simulate;
    
    for c=1:n_columns
        for nn=1:n_excitatory
            psth_avg_exc(c,nn,b) = sum(spikes(iter,c,nn,:))/(1700 * 0.001);
        end
    end
end
som_off_psth_avg_exc = psth_avg_exc;

%% compare rates
col = 2;
som_on_col_psth_exc = squeeze(som_on_psth_avg_exc(col,:,:));
som_on_diff_col_psth_exc = diff(som_on_col_psth_exc, 1, 2);
mean_som_on_diff_col_psth_exc = mean(som_on_diff_col_psth_exc, 1);

som_off_col_psth_exc = squeeze(som_off_psth_avg_exc(col,:,:));
som_off_diff_col_psth_exc = diff(som_off_col_psth_exc, 1, 2);
mean_som_off_diff_col_psth_exc = mean(som_off_diff_col_psth_exc, 1);

[h, p] = ttest(mean_som_on_diff_col_psth_exc, mean_som_off_diff_col_psth_exc)