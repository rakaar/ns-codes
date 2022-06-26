%% single batch
weights_stacked_side_by_side = zeros(num_network_neurons*num_network_neurons, length(tspan));
for t=1:length(tspan)
    weights_stacked_side_by_side(:,t) = reshape(network_weight_matrix(1,t,:,:),  num_network_neurons*num_network_neurons,1);
end

figure
    imagesc(weights_stacked_side_by_side);
    title('weights side by side all network neurons')
grid

%% within column 
n_columns=5;

within_column_weight_matrices = zeros(n_columns, length(tspan),n_excitatory, n_excitatory);
within_column_weight_matrices_diagonal_removed = zeros(n_columns, length(tspan),n_excitatory, n_excitatory-1);

for c1=1:n_columns
    for c2=1:n_columns
        if c1 == c2
            
               n1 = (c1-1)*25 + 1;
               n2 = n1 + 19;
               within_column_weights_tensor = network_weight_matrix(1,:,n1:n2,n1:n2);
               within_column_weights_tensor_squeezed = squeeze(within_column_weights_tensor);
               within_column_weights_tensor_squeezed_dia_removed = remove_diagonal_elements(within_column_weights_tensor_squeezed);
               
               within_column_weight_matrices(c1, :, :, :) = reshape(within_column_weights_tensor_squeezed,  1,length(tspan),n_excitatory,n_excitatory);
               within_column_weight_matrices_diagonal_removed(c1,:,:,:) = reshape(within_column_weights_tensor_squeezed_dia_removed,  1,length(tspan),n_excitatory, n_excitatory-1);

                weights_stacked_side_by_side = zeros(n_excitatory*n_excitatory, length(tspan));
                weights_stacked_side_by_side_diagonal_removed = zeros(n_excitatory*(n_excitatory-1), length(tspan));
                for t=1:length(tspan)
                    weights_stacked_side_by_side(:,t) = reshape(within_column_weights_tensor_squeezed(t,:,:),  n_excitatory*n_excitatory,1);
                    weights_stacked_side_by_side_diagonal_removed(:,t) = reshape(within_column_weights_tensor_squeezed_dia_removed(t,:,:), (n_excitatory-1)*n_excitatory,1);
                end

            
%            figure
%                 imagesc(weights_stacked_side_by_side);
%                 title(['col',num2str(c1), ' weight matrix'])
%             grid

           figure
                imagesc(weights_stacked_side_by_side_diagonal_removed);
                title(['col',num2str(c1), ' weight matrix'])
            grid
        else % c1 not equal to c2
            continue
        end
    end
end

%% across column 
close all
iter = 1;
pre = 15; 
col_pre = 2;
col_post = 1;

pre_n_index = (col_pre - 1)*25 + pre;
for post=1:20
    post_n_index = (col_post-1)*25 + post;
    figure
        plot(squeeze(network_weight_matrix(iter, :,pre_n_index,post_n_index)))
        title(['pre ', num2str(pre), ' post ', num2str(post)])
    grid
end

%% within column individual synapses
close all
iter = 1;
pre = 5; 

col = 4;
for post=1:20
    figure
        plot(squeeze(within_column_weight_matrices(col, :,pre,post)))
        title(['pre ', num2str(pre), ' post ', num2str(post)])
    grid
end