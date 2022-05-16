close all;
iter_to_see = 5;
col = 1;
for n=1:n_total_neurons
    xe_vec = squeeze(xe(iter_to_see,col,n,:));
    figure(n)
        plot(xe_vec);
        title('xe of neuron')
    grid
end