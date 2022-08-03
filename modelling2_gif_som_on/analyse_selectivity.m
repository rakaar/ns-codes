close all;

bin_size = 10;
for c=1:n_columns
    spikes_col = squeeze(spikes(1,c,:,1:t_simulate+1));
    spikes_col_mean = mean(spikes_col, 1);
    spike_rate_10 = spikes_to_spike_rate_neat(spikes_col_mean, physical_time_in_ms, dt, bin_size*dt);

    figure
        plot(spike_rate_10);
        title(['spikes-col-', num2str(c)])
    grid  
end

%% same diff style
% for c=1:n_columns
%     spikes_col = squeeze(spikes(1,c,:,1:t_simulate));
%     spikes_col_mean = mean(spikes_col, 1);
% 
% 
%     spikes_col_mean_reshape = reshape(spikes_col_mean,   bin_size, t_simulate/bin_size);
%     spikes_col_mean_reshape_avg = mean(spikes_col_mean_reshape, 1);
%     figure
%         plot(spikes_col_mean_reshape_avg);
%         title(['spikes-col-', num2str(c)])
%     grid  
% end
%% seperate responses towards a and b
token_size = 170;
spike_rate_a = zeros(n_columns, n_excitatory, n_tokens);
spike_rate_b = zeros(n_columns, n_excitatory, n_tokens);
token_num = 1;
for t=1:t_simulate
    if mod(t,token_size) == 0
        token_num = token_num + 1;
    end
    
    if lamda(1,4,1,t) == 300 % a 
        for c=1:n_columns
            for n=1:n_excitatory
                spike_rate_a(c,n,token_num) = spike_rate_a(c,n,token_num) + spikes(1,c,n,t);
            end
        end
    elseif lamda(1,4,1,t) == 50 % b
        for c=1:n_columns
            for n=1:n_excitatory
                spike_rate_b(c,n,token_num) = spike_rate_b(c,n,token_num) + spikes(1,c,n,t);
            end
        end
    end
end

for c=1:n_columns
    spike_rate_a_all_neurons = squeeze(spike_rate_a(c,:,:));
    spike_rate_b_all_neurons = squeeze(spike_rate_b(c,:,:));

    spike_rate_a_mean_over_neurons = mean(spike_rate_a_all_neurons, 1);
    spike_rate_b_mean_over_neurons = mean(spike_rate_b_all_neurons, 1);

    figure
        hold on
            plot(spike_rate_a_mean_over_neurons)
            plot(spike_rate_b_mean_over_neurons)
        hold off

        title(['spike rate a and b-col-', num2str(c)])
        legend('spike rate a', 'spike rate b')
    grid
end

col = 4;
figure
    hold on
        imagesc(squeeze(spike_rate_a(col,:,:)))
    hold off
    title(['spike rate a-col-', num2str(col)])
grid

figure
    hold on
        imagesc(squeeze(spike_rate_b(col,:,:)))
    hold off
    title(['spike rate b-col-', num2str(col)])
grid