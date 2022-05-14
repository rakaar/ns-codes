close all;
iter_to_see = 1;
first_spike_arr = zeros(1,n_total_neurons);
for n=1:n_total_neurons
    x = spikes(iter_to_see, col, n, :);
    x = squeeze(x);
    z = find(x == 1);
    first_spike_arr(1,n) = z(1);
end

figure
    stem(first_spike_arr);
grid