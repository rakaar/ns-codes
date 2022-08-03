close all
iter_to_see = 2;
col = 1;
for n=1:n_total_neurons
    figure(n)
        g = get_g_t_vector(squeeze(spikes(iter_to_see,col,n,:)), length(tspan));
        x = squeeze(xe(iter_to_see,col,n,:));
        x = x';
        current = g.*x;
        hold on
            plot(current);
            stem(squeeze(spikes(iter_to_see,col,n,:))/10);
        hold off
        
    grid
end