close all
bin_size = 10;
overlap_size = 5;
n_batches = 200;
n_excitatory = 20;
batch_data_path = "D:\7_multi_col_big_clip_range";

gini_coeff_row_over_time = zeros(1,n_batches);
gini_coeff_col_over_time = zeros(1,n_batches);
gini_coeff_of_transpose_and_sum_row_over_time = zeros(1,n_batches);
gini_coeff_of_transpose_and_sum_col_over_time = zeros(1,n_batches);

iter=1;
for col=1:5

for t=1:n_batches
    fprintf("\n for batch %d \n",t);
    weights_avg_matrix = zeros(n_excitatory, n_excitatory);
    batch_file_name = batch_data_path + "\batch_" + num2str(t) + ".mat";
    batch_weight_matrix_struct = load(batch_file_name, "exc_to_exc_weight_matrix");
    batch_weight_matrix = batch_weight_matrix_struct.exc_to_exc_weight_matrix;
    batch_weight_matrix1 = squeeze(batch_weight_matrix(iter,col,:,:,:));
    batch_tspan = size(batch_weight_matrix1, 1);
    for presyn=1:n_excitatory
        for postsyn=1:n_excitatory
            weights_avg_matrix(presyn, postsyn) = sum(batch_weight_matrix1(:,presyn,postsyn))/batch_tspan;
        end
    end
    
    % over row
    gini_coeff = calculate_gini_coeff(weights_avg_matrix,0);
    gini_coeff_row_over_time(1,t) = gini_coeff;
    
    % over column
    gini_coeff = calculate_gini_coeff(weights_avg_matrix,1);
    gini_coeff_col_over_time(1,t) = gini_coeff;

    % of matrix = w + w'
    weights_avg_matrix_sum_with_transpose = weights_avg_matrix + transpose(weights_avg_matrix);
    gini_coeff_0 = calculate_gini_coeff(weights_avg_matrix_sum_with_transpose,0);
    gini_coeff_1 = calculate_gini_coeff(weights_avg_matrix_sum_with_transpose,1);

    gini_coeff_of_transpose_and_sum_row_over_time(1,t) = gini_coeff_0;
    gini_coeff_of_transpose_and_sum_col_over_time(1,t) = gini_coeff_1;
    
end

figure
    hold on
        plot(gini_coeff_row_over_time);
        plot(gini_coeff_col_over_time);
        plot(gini_coeff_of_transpose_and_sum_row_over_time);
%         plot(gini_coeff_of_transpose_and_sum_col_over_time, '--.');
    hold off
    title(['gini coeffs', num2str(col)])
%     legend('row', 'col', 'a + transpose(a) row','a + transpose(a) col')
        legend('row', 'col', 'a + transpose(a) row')
    
grid

end

