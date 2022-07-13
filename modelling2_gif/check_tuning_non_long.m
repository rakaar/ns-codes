%% thalamic epsc given by 9 cols
% epsc_thalamic = zeros(n_iters, n_thalamic_cols, n_thalamic_neurons, length(tspan));

close all
iter=1;
bin_size = 50;
for c=1:n_thalamic_cols
    figure
%         plot(  mean(squeeze(epsc_thalamic(iter,c,:,1:t_simulate)), 1)  ); % without bin
        squeezed_epsc = squeeze(epsc_thalamic(iter,c,:,1:t_simulate));
        mean_over_neurons_epsc = mean(squeezed_epsc, 1);
        mean_epsc_reshaped = reshape(mean_over_neurons_epsc,  bin_size, t_simulate/bin_size);
        mean_epsc_binned = mean(mean_epsc_reshaped, 1);
        plot(mean_epsc_binned)
        title(['thalamic-epsc-avg-over-neurons-',num2str(c)])
    grid
end

%% input given to each n
close all;
iter=1;
% thalamic_epsc_to_neuron_thalamic_column_wise = zeros(n_iters, n_columns, n_total_neurons, n_thalamic_cols);
for c=1:n_columns
    thalamic_epsc_to_c_col = zeros(1, n_thalamic_cols);
    for c_thal=1:n_thalamic_cols
        thalamic_epsc_to_c_col(1, c_thal) = mean(squeeze(thalamic_epsc_to_neuron_thalamic_column_wise(iter,c,:,c_thal)));
    end
    
    figure
        plot(thalamic_epsc_to_c_col)
        title(['thalamic inputs to col-', num2str(c)])
    grid
end