close all;
n_total_neurons = 25;
num_network_neurons = 125;
batches = 30;
tspan = load("D:\2_multi_col_across_plastic\batch_1.mat", "tspan").tspan;

batch_avg_network_weights = zeros(batches, num_network_neurons, num_network_neurons);
% all_1ms_network_weights = zeros((batches*(length(tspan)-1))/10, num_network_neurons, num_network_neurons);
n_exc = 20;

batch_data_path = "D:\2_multi_col_across_plastic";
images_path = "D:\2_multi_col_images\";


iter=1;

for b=1:batches
    fprintf("\n batch is %d \n",b)
    batch_file_name = batch_data_path + "\batch_" + num2str(b) + ".mat";
    network_weight_matrix = load(batch_file_name, "network_weight_matrix").network_weight_matrix;
    for n1=1:num_network_neurons
        for n2=1:num_network_neurons
             batch_avg_network_weights(b,n1,n2) = mean(network_weight_matrix(iter,:,n1,n2));
%              all_1ms_network_weights((b-1)*length(tspan) + 1:(b-1)*length(tspan) + length(tspan)   ,n1,n2) =  network_weight_matrix(iter,:,n1,n2);
        end
    end


end

% plot weight matrix
n_columns = 5;
for c1=1:n_columns
    for c2=1:n_columns
        if c1 - c2 == 0
            % within column matrix
            within_column_matrix = batch_avg_network_weights(:,(c1-1)*n_total_neurons + 1:(c1-1)*n_total_neurons + n_exc, (c1-1)*n_total_neurons + 1:(c1-1)*n_total_neurons + n_exc);
            fprintf("\n %d %d %d %d \n ",(c1-1)*n_total_neurons + 1,(c1-1)*n_total_neurons + n_exc, (c1-1)*n_total_neurons + 1,(c1-1)*n_total_neurons + n_exc)
            reshaped_within_column_matrix = reshape(within_column_matrix, batches, n_exc*n_exc);
            figure
                plot(reshaped_within_column_matrix)
                title(['weights batches wise', num2str(c1)])
                image_name = images_path + "w_batch_avg-c-" + num2str(c1) + ".fig"; 
                saveas(gcf, image_name);

            grid
        elseif abs(c1 - c2) == 1 || abs(c1 - c2) == 2
            % c1 pre ,c2 post
            across_column_matrix = batch_avg_network_weights(:,(c1-1)*n_total_neurons + 1:(c1-1)*n_total_neurons + n_exc, (c2-1)*n_total_neurons + 1:(c2-1)*n_total_neurons + n_exc);
            fprintf("\n %d %d %d %d \n ",(c1-1)*n_total_neurons + 1,(c1-1)*n_total_neurons + n_exc, (c2-1)*n_total_neurons + 1,(c2-1)*n_total_neurons + n_exc)
            reshaped_across_column_matrix = reshape(across_column_matrix, batches, n_exc*n_exc);
            figure
                plot(reshaped_across_column_matrix)
                title(['weights batches wise', num2str(c1), '-',num2str(c2)])
                image_name = images_path + "w_batch_avg-c1-" + num2str(c1) + "-c2-" + num2str(c2)+ ".fig"; 
                saveas(gcf, image_name);
            grid

        else
            continue
        end
    end
end


function new_vec = make_bins_of_weights(vec, bin_size)
    time_span = length(vec) - 1;
    new_vec = zeros(time_span/10,1);
    for i=1:bin_size:time_span-9
        new_vec(i,1) = mean(vec(1,i:i+9));
    end
end