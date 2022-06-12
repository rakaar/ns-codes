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

    figureg
        for i=1:length(neurons)
                        subplot(neuron_len,1,i)
                        plot(squeeze(spikes(1,2,neurons(i),:)))
                        title(['spikes ', num2str(neurons(i))])

        end
    grid



