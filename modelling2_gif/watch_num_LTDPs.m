close all;clear all;
iter = 1;
col = 1;
batches = 200;
n_excitatory = 20;
batch_data_path = "D:\batches_data";

num_of_LTPs = 0;
num_of_LTDs_same_time_spike = 0;
num_of_LTDs_non_same_time_spike = 0;

for b=1:batches
        fprintf("\n batch is %d \n", b);
        batch_file_name = batch_data_path + "\batch_" + num2str(b) + ".mat";
        
        batch_spikes = load(batch_file_name, "spikes").spikes;
        batch_weight_matrix = load(batch_file_name, "exc_to_exc_weight_matrix").exc_to_exc_weight_matrix;
        
        batch_spikes = squeeze(batch_spikes);
        batch_weight_matrix = squeeze(batch_weight_matrix);

        
        time_period = size(batch_spikes,2);

        for t=1:time_period-1
            for presyn=1:n_excitatory
                for postsyn=1:n_excitatory
                    % LTP or LTD
                    if batch_weight_matrix(t,presyn,postsyn) < batch_weight_matrix(t+1,presyn,postsyn)
                        num_of_LTPs = num_of_LTPs + 1;
                    elseif batch_weight_matrix(t,presyn,postsyn) > batch_weight_matrix(t+1,presyn,postsyn)
                             % is it due to same time spike
                            if batch_spikes(presyn,t) == batch_spikes(postsyn,t)
                                num_of_LTDs_same_time_spike = num_of_LTDs_same_time_spike + 1;
                            else
                                num_of_LTDs_non_same_time_spike = num_of_LTDs_non_same_time_spike + 1 ;
                            end
                    end
                end
            end
        end
        
end