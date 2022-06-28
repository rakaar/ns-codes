%% num spikes
close all;
batches = 1:40;
n_excitatory = 20;
iter = 1;
inc_neuron = 4; ctrl_neuron = 3;
inc_col = 1;
num_of_spikes_inc_neuron = 0;
num_of_spikes_ctrl_neuron = 0;
num_of_spikes = zeros(n_excitatory, length(batches));
for batch_num=batches
    fprintf("\n batch num is %d \n", batch_num);
    batch_data_path = "D:\7_multi_col_big_clip_range";
    batch_file_name = batch_data_path + "\batch_" + num2str(batch_num) + ".mat";
  
    spikes = load(batch_file_name,"spikes").spikes;
    
    
    for n=1:n_excitatory
        squeezed_spikes_n = squeeze(spikes(iter,inc_col,n,:));
        spike_count = length(find(squeezed_spikes_n == 1));
        num_of_spikes(n,batch_num) = num_of_spikes(n,batch_num) + spike_count;
    end
end

figure
    plot(transpose(num_of_spikes));
    title('num of spikes')
grid

figure
    imagesc(num_of_spikes);
    title('num of spikes')
grid


%% synapse 1 -> 3, 1 -> 4
close all;
batches = 1:40;
n_excitatory = 20;
iter = 1;
inc_neuron = 4; ctrl_neuron = 3;
inc_col = 1;
num_of_spikes_inc_neuron = 0;
num_of_spikes_ctrl_neuron = 0;
pre_neuron_test = 6;

weights_to_n = zeros(n_excitatory, length(batches));
for batch_num=batches
    fprintf("\n batch num is %d \n", batch_num);
    batch_data_path = "D:\7_multi_col_big_clip_range";
    batch_file_name = batch_data_path + "\batch_" + num2str(batch_num) + ".mat";
  
    exc_to_exc_weight_matrix = load(batch_file_name,"exc_to_exc_weight_matrix").exc_to_exc_weight_matrix;
    for n=1:n_excitatory
        weights_to_n(n,batch_num) = mean(exc_to_exc_weight_matrix(iter,inc_col,:,pre_neuron_test,n));
    end
end 

figure
    plot(transpose(weights_to_n));
    title('weights to n')
grid

figure
    imagesc(weights_to_n);
    title('weights to n')
grid
%% number of times other neuron didn't have spikes where there was a spike 
close all;
batches = 1:40;
n_excitatory = 20;
iter = 1;
inc_col = 1;

num_of_times_spike_was_not_there = zeros(n_excitatory, length(batches));

for batch_num=batches
    fprintf("\n batch num is %d \n", batch_num);
    batch_data_path = "D:\7_multi_col_big_clip_range";
    batch_file_name = batch_data_path + "\batch_" + num2str(batch_num) + ".mat";
  
    exc_to_exc_weight_matrix = load(batch_file_name,"exc_to_exc_weight_matrix").exc_to_exc_weight_matrix;
    tspan = load(batch_file_name,"tspan").tspan;
    spikes = load(batch_file_name,"spikes").spikes;


    for n=1:n_excitatory
        for t=2:length(tspan)
            oldest_time = t-20;
            if oldest_time < 1
                oldest_time = 1;
            end
            if spikes(iter,inc_col,n,t) == 1
                for pre=1:n_excitatory
                    if spikes(iter,inc_col,pre,t) == 0
                        for ttt=t-1:-1:oldest_time
                            if spikes(iter,inc_col,pre,ttt) == 1
                                num_of_times_spike_was_not_there(n,batch_num) = num_of_times_spike_was_not_there(n,batch_num) + 1;
                                break
                            end
                        end
                    end
                end
            end
        end
    end
end 

figure
    plot(transpose(num_of_times_spike_was_not_there));
    title('spikes not there')
grid

figure
    imagesc(num_of_times_spike_was_not_there);
    title('spikes not there')
grid
%% ratio
close all
the_ratio = zeros(n_excitatory, length(batches));

for b=1:length(batches)
    for n=1:n_excitatory
        the_ratio(n,b) = num_of_times_spike_was_not_there(n,b)/num_of_spikes(n,b);
    end
end
figure
    plot(transpose(the_ratio));
    title('ratio')
grid