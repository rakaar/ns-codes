% there is fast code, ignore this
close all;
batch_data_path = "D:\1_multi_col_across_plastic";
n_excitatory = 20;
n_total_neurons = 25;
batches = 50;
n_columns = 5;
iter = 1;
tspan = load("D:\1_multi_col_across_plastic\batch_1.mat", "tspan").tspan;

for c1=1:n_columns
    for c2=1:n_columns
        if c1 == c2
            % within column
            weights_1ms = zeros(batches*length(tspan),n_excitatory,n_excitatory);
            weights_batch = zeros(batches,n_excitatory,n_excitatory);

            for b=1:batches
                fprintf("\n col is %d. batch is %d\n", c1,b)
                batch_file_name = batch_data_path + "\batch_" + num2str(b) + ".mat";
                network_weight_matrix = load(batch_file_name, "network_weight_matrix").network_weight_matrix;
                batch_weight_tensor = zeros(length(tspan), n_excitatory, n_excitatory);
                for n1=1:n_excitatory
                    for n2=1:n_excitatory
                        n1_index_in_network = (c1-1)*n_total_neurons + n1;
                        n2_index_in_network = (c1-1)*n_total_neurons + n2;
                        fprintf("\n col-%d, n1-%d, network-%d, n2-%d, network-%d \n",c1,n1,n1_index_in_network,n2,n2_index_in_network)
                        batch_weight_tensor(:,n1,n2) = reshape(network_weight_matrix(iter,:,n1_index_in_network,n2_index_in_network),   length(tspan),1);
                        weights_batch(b,n1,n2) = mean(batch_weight_tensor(:,n1,n2));
%                         fprintf("\n b-%d,n1-%d,n2-%d,mean %f \n",b,n1,n2,weights_batch(b,n1,n2))
%                         if b == 1 && c1 == 1 
%                             fprintf("\n mean is %f \n",mean(batch_weight_tensor(:,n1,n2)))
%                             pause(1)
%                         end
                        weights_1ms((b-1)*length(tspan) + 1:(b-1)*length(tspan) + length(tspan),n1,n2) = reshape(batch_weight_tensor(:,n1,n2)   ,length(tspan),1);
                    end
                end
            end

            figure
                imagesc(transpose(reshape(weights_batch,   batches,n_excitatory*n_excitatory)))
                title(['1ms-within weights-col- ', num2str(c1)]);
            grid

              figure
                plot(transpose(reshape(weights_1ms,   batches*length(tspan),n_excitatory*n_excitatory)))
                title(['batchavg-within weights-col- ', num2str(c1)]);
            grid
            
        elseif abs(c1-c2) == 1
            continue
        elseif abs(c1-c2) == 2
            continue
        else
            continue
        end
    end
end
