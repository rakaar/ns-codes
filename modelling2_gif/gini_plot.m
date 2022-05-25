close all
bin_size = 10;
overlap_size = 5;
bin_starting_times=1:overlap_size:length(tspan)-(bin_size-1);
binned_matrices = zeros(length(bin_starting_times), n_excitatory, n_excitatory);

gini_coeff_row_over_time = zeros(1,length(bin_starting_times));
gini_coeff_col_over_time = zeros(1,length(bin_starting_times));
gini_coeff_of_transpose_and_sum_row_over_time = zeros(1,length(bin_starting_times));
gini_coeff_of_transpose_and_sum_col_over_time = zeros(1,length(bin_starting_times));

iter=1;col=1;
bin_index = 1;
for t=bin_starting_times
    weights_avg_matrix = zeros(n_excitatory, n_excitatory);
    weights_avg_matrix(:,:) = reshape(sum(exc_to_exc_weight_matrix(iter,col,t:t+(bin_size-1),:,:)), n_excitatory, n_excitatory)/bin_size;
    binned_matrices(bin_index, :, :) = reshape(weights_avg_matrix, 1,n_excitatory,n_excitatory);
    
    % over row
    gini_coeff = calculate_gini_coeff(weights_avg_matrix,0);
    gini_coeff_row_over_time(1,bin_index) = gini_coeff;
    
    % over column
    gini_coeff = calculate_gini_coeff(weights_avg_matrix,1);
    gini_coeff_col_over_time(1,bin_index) = gini_coeff;

    % of matrix = w + w'
    weights_avg_matrix_sum_with_transpose = weights_avg_matrix + transpose(weights_avg_matrix);
    gini_coeff_0 = calculate_gini_coeff(weights_avg_matrix_sum_with_transpose,0);
    gini_coeff_1 = calculate_gini_coeff(weights_avg_matrix_sum_with_transpose,1);

    gini_coeff_of_transpose_and_sum_row_over_time(1,bin_index) = gini_coeff_0;
    gini_coeff_of_transpose_and_sum_col_over_time(1,bin_index) = gini_coeff_1;
    
    bin_index = bin_index + 1;
end

figure
    hold on
        plot(gini_coeff_row_over_time);
        plot(gini_coeff_col_over_time);
        plot(gini_coeff_of_transpose_and_sum_row_over_time);
%         plot(gini_coeff_of_transpose_and_sum_col_over_time, '--.');
    hold off
    title('gini coeffs')
%     legend('row', 'col', 'a + transpose(a) row','a + transpose(a) col')
        legend('row', 'col', 'a + transpose(a) row')
    
grid