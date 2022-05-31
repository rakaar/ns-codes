close all;clear all;
iter = 1;
col = 1;
batches = 200;
n_excitatory = 20;
% batch_data_path = "D:\batches_data";
batch_data_path = "D:\batches_data";

num_LTP_same_time_spike = 0;
num_LTP_diff_time_spike = 0;

num_LTD_same_time_spike = 0;
num_LTD_diff_time_spike = 0;

same_time_spike_LTP = 0; 
same_time_spike_LTD = 0;


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
                    % same time spike check
                    if (batch_spikes(presyn,t) == 1) && (batch_spikes(postsyn,t) == 1)
                        if batch_weight_matrix(t,presyn,postsyn) < batch_weight_matrix(t+1,presyn,postsyn)
                            same_time_spike_LTP = same_time_spike_LTP + 1;
                        elseif batch_weight_matrix(t,presyn,postsyn) > batch_weight_matrix(t+1,presyn,postsyn)
                            same_time_spike_LTD = same_time_spike_LTD + 1;
                        end
                    end


                    % LTP or LTD
                    if batch_weight_matrix(t,presyn,postsyn) < batch_weight_matrix(t+1,presyn,postsyn)
                        % LTP
                        % same time spike or not
                        if (batch_spikes(presyn,t) == 1) && (batch_spikes(postsyn,t) == 1)
                            num_LTP_same_time_spike = num_LTP_same_time_spike + 1;
                        else
                            num_LTP_diff_time_spike = num_LTP_diff_time_spike + 1;
                        end

                    elseif batch_weight_matrix(t,presyn,postsyn) > batch_weight_matrix(t+1,presyn,postsyn)
                         if (batch_spikes(presyn,t) == 1) && (batch_spikes(postsyn,t) == 1)
                            num_LTD_same_time_spike = num_LTD_same_time_spike + 1;
                        else
                            num_LTD_diff_time_spike = num_LTD_diff_time_spike + 1;
                        end
                    end
                end
            end
        end
        
end

fprintf("\n num_LTP_same_time_spike %d \n",num_LTP_same_time_spike)
fprintf("\n num_LTP_diff_time_spike %d \n",num_LTP_diff_time_spike)
fprintf("\n num_LTD_same_time_spike %d \n",num_LTD_same_time_spike)
fprintf("\n num_LTD_diff_time_spike %d \n", num_LTD_diff_time_spike)

fprintf("\n same_time_spike_LTP %d \n", same_time_spike_LTP)
fprintf("\n same_time_spike_LTD %d \n", same_time_spike_LTD)