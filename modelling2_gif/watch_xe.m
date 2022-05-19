close all;
iter_to_see = 5;
col = 2;
for n=1:n_thalamic
    xe_vec = squeeze(xe_thalamic_std(n,:));
    figure(n)
        plot(xe_vec);
        title('xe of neuron')
    grid
end