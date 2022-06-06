%% 1 ms bins
clear all;close all;
pre_syn = 1;
post_syn = 2;
iter = 1;
col = 1;
weight_over_time = [];
batches = 200;
% batch_data_path = "D:\batches_data";
batch_data_path = "D:\4_batches_data";

for b=1:batches
        fprintf("\n batch is %d \n", b);
        batch_file_name = batch_data_path + "\batch_" + num2str(b) + ".mat";
        batch_weight_matrix_struct = load(batch_file_name, "exc_to_exc_weight_matrix");
        batch_weight_matrix = batch_weight_matrix_struct.exc_to_exc_weight_matrix;
        weights_over_time_current_batch = squeeze(batch_weight_matrix(iter,col,:,pre_syn, post_syn));
        weight_over_time = [weight_over_time; weights_over_time_current_batch];
end

figure
    plot(weight_over_time)
    title("presyn is "+ pre_syn + "," + "postsyn is " + post_syn)
grid
%% 200 ms bins
clear all;close all;
pre_syn = 10;
post_syn = 15;
iter = 1;
col = 1;
weights_over_time = [];
batches = 200;
% batch_data_path = "D:\batches_data";
batch_data_path = "D:\4_batches_data";

for b=1:batches
        fprintf("\n batch is %d \n", b);
        batch_file_name = batch_data_path + "\batch_" + num2str(b) + ".mat";
        batch_weight_matrix_struct = load(batch_file_name, "exc_to_exc_weight_matrix");
        batch_weight_matrix = batch_weight_matrix_struct.exc_to_exc_weight_matrix;
        weights_over_time_current_batch = squeeze(batch_weight_matrix(iter,col,:,pre_syn, post_syn));
        weights_over_time_current_batch_avg = sum(weights_over_time_current_batch)/length(weights_over_time_current_batch);
        weights_over_time = [weights_over_time, weights_over_time_current_batch_avg];
        
end

figure
    plot(weights_over_time)
    title("presyn is "+ pre_syn + "," + "postsyn is " + post_syn)
grid

%% imagesc
close all;clear all;
iter = 1;
col = 1;
batches = 200;
n_excitatory = 20;
% batch_data_path = "D:\batches_data";
batch_data_path = "D:\4_batches_data";
weights_stacked_side_by_side = zeros(n_excitatory*n_excitatory, batches);
for b=1:batches
        fprintf("\n batch is %d \n", b);
        batch_file_name = batch_data_path + "\batch_" + num2str(b) + ".mat";
        batch_weight_matrix_struct = load(batch_file_name, "exc_to_exc_weight_matrix");
        batch_weight_matrix = batch_weight_matrix_struct.exc_to_exc_weight_matrix;
        batch_weight_matrix = squeeze(batch_weight_matrix);

        
        batch_weight_matrix_last = batch_weight_matrix(end, :,:);
        batch_weight_matrix_last = squeeze(batch_weight_matrix_last);
        batch_weight_matrix_reshaped = reshape(batch_weight_matrix_last, n_excitatory*n_excitatory,1);

        weights_stacked_side_by_side(:,b) = batch_weight_matrix_reshaped; 
end

figure
    imagesc(weights_stacked_side_by_side);
    title('weights at end of batch stacked side by side')
grid
