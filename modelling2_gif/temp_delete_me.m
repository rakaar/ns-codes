%% check 1
neurons = [1 5 7 9 11 13 15 15 17 19];
neuron_len = length(neurons);

    figure
        for i=1:length(neurons)
                        subplot(neuron_len,1,i)
                        plot(squeeze(spikes(1,2,neurons(i),:)))
                        title(['spikes ', num2str(neurons(i))])

        end
    grid

%% check 2
neurons = [1 5 7 9 11 13 15  17 19];
neuron_len = length(neurons);

    figure
        for i=1:length(neurons)
                        subplot(neuron_len,1,i)
                        plot(squeeze(spikes(1,2,neurons(i),:)))
                        title(['spikes ', num2str(neurons(i))])

        end
    grid

%% check 3
neurons = [2 4 6 8 10 12 14 16 18 20 1 5 7 9 11 13 15  17 19];
neuron_len = length(neurons);

    figure
        for i=1:length(neurons)
                        subplot(neuron_len,1,i)
                        plot(squeeze(spikes(1,2,neurons(i),:)))
                        title(['spikes ', num2str(neurons(i))])

        end
    grid


%% check epsc
close all;
total_input_epsc = thalamic_epsc_tensor ...
                    + recurrence_exc_self_column_epsc_tensor + recurrence_inh_self_column_epsc_tensor ...
                    + recurrence_exc_neighbour_column_epsc_tensor + recurrence_inh_neighbour_column_epsc_tensor ;

neurons = [2 4 6 8 10 12 14 16 18 20 1 5 7 9 11 13 15  17 19];
neuron_len = length(neurons);

    figure
        for i=1:length(neurons)
                        subplot(neuron_len,1,i)
                        plot(squeeze(total_input_epsc(1,4,neurons(i),:)))
                        title(['spikes ', num2str(neurons(i))])

        end
    grid


%% check thalamic epsc
close all;

neurons = [2 4 6 8 10 12 14 16 18 20 1 5 7 9 11 13 15  17 19];
neuron_len = length(neurons);

    figure
        for i=1:length(neurons)
                        subplot(neuron_len,1,i)
                        plot(squeeze(thalamic_epsc_tensor(1,3,neurons(i),:)))
                        title(['spikes ', num2str(neurons(i))])

        end
    grid


%% thalamic epsc
close all
figure
i=1;
for c=1:9
        subplot(9,2,i)
        plot(squeeze(epsc_thalamic(1,c,1,:)))

        subplot(9,2,i+1)
        plot(squeeze(epsc_thalamic(1,c,2,:)))
        i = i + 2;
   end
grid

%% thalamic epsc
close all
figure
   for c=1:9
        subplot(9,2,c)
        plot(squeeze(epsc_thalamic(1,c,2,:)))
   end
grid