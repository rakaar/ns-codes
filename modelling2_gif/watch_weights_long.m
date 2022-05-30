%% 1 ms bins
clear all;close all;
pre_syn = 1;
post_syn = 2;
iter = 1;
col = 1;
weight_over_time = [];
batches = 200;
batch_data_path = "D:\batches_data";

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
batch_data_path = "D:\batches_data";

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
