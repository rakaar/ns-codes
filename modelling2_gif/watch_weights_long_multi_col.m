%% 1 ms bins
clear all;close all;
pre_syn = 10;
post_syn = 15;
iter = 1;
col = 3;
weight_over_time = [];
batches = 200;
batch_data_path = "D:\1_mult_col_pvsom";

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
col = 4;
weights_over_time = [];
batches = 200;
batch_data_path = "D:\1_mult_col_pvsom";

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
        title("presyn is "+ pre_syn + "," + "postsyn is " + post_syn + " col " + num2str(col))
grid

%% all weights over time in batches
clear all;close all;
iter = 1;
n_cols = 5;
batches = 200;
batch_data_path = "D:\1_mult_col_pvsom";
n_excitatory = 20;

batch_file_name = batch_data_path + "\batch_" + num2str(1) + ".mat";
batch_weight_matrix_struct = load(batch_file_name, "exc_to_exc_weight_matrix");
batch_weight_matrix = batch_weight_matrix_struct.exc_to_exc_weight_matrix;
time_span = size(batch_weight_matrix,3);

weights_over_time_avg = zeros(n_cols, n_excitatory*n_excitatory,batches);

for b=1:batches
        fprintf("\n batch is %d \n", b);

        batch_file_name = batch_data_path + "\batch_" + num2str(b) + ".mat";
        batch_weight_matrix_struct = load(batch_file_name, "exc_to_exc_weight_matrix");
        batch_weight_matrix = batch_weight_matrix_struct.exc_to_exc_weight_matrix;

    for c=1:n_cols
        for presyn=1:n_excitatory
            for postsyn=1:n_excitatory
                weights_over_time_avg(c,(presyn-1)*20 + postsyn, b) = mean(squeeze(batch_weight_matrix(iter,c,:,presyn, postsyn)));
            end
        end
    end

end

for c=1:n_cols
    figure
        imagesc(squeeze(weights_over_time_avg(c,:,:)));
        title(['col ',num2str(c)])
    grid
end

%% imagesc
close all;clear all;
iter = 1;
batches = 200;
n_excitatory = 20;
n_cols = 5;
batch_data_path = "D:\1_mult_col_pvsom";
weights_stacked_side_by_side = zeros(n_cols, n_excitatory*n_excitatory, batches);
for b=1:batches
        fprintf("\n batch is %d \n", b);
        batch_file_name = batch_data_path + "\batch_" + num2str(b) + ".mat";
        batch_weight_matrix_struct = load(batch_file_name, "exc_to_exc_weight_matrix");
        batch_weight_matrix = batch_weight_matrix_struct.exc_to_exc_weight_matrix;
        
        
        for col=1:5
            batch_weight_matrix1 = squeeze(batch_weight_matrix(iter,col,:,:,:));
            batch_weight_matrix_last = batch_weight_matrix1(end, :,:);
            batch_weight_matrix_last1 = squeeze(batch_weight_matrix_last);
            batch_weight_matrix_reshaped = reshape(batch_weight_matrix_last1, n_excitatory*n_excitatory,1);
            
            weights_stacked_side_by_side(col,:,b) = batch_weight_matrix_reshaped; 
        end
        
end

for c=1:n_cols
    figure
            imagesc(squeeze(weights_stacked_side_by_side(c,:,:)));
            title(['stack side by side weights ',num2str(c)])
    grid
end

