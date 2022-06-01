close all;
%% presyn fix
close all;
pre_syn = 5;
iter = 1;
col = 1;
for post_syn=1:n_excitatory
    figure(post_syn)
    hold on
    weight_over_time = squeeze(exc_to_exc_weight_matrix(iter,col,:,pre_syn, post_syn));
    plot(weight_over_time,'--p')
    title("presyn is "+ pre_syn + "," + "postsyn is " + post_syn)
    %         plot(protochol);
    hold off
    grid
end

%% postsyn fix
close all;
post_syn = 13;
iter = 1;
col = 1;
for pre_syn=1:n_excitatory
    figure(pre_syn)
    hold on
    weight_over_time = squeeze(exc_to_exc_weight_matrix(iter,col,:,pre_syn, post_syn));
    plot(weight_over_time,'--p')
    title("presyn is "+ pre_syn + "," + "postsyn is " + post_syn)
    %         plot(protochol);
    hold off
    grid
end

%% watch number of ltps and ltds over time
close all;
iter = 1; col = 1;
figure
    hold on
    plot(protochol)
    plot(squeeze(num_of_LTPs(iter,col,:)))
    title('num of ltps over time')
    legend('protochol', 'ltps')
    hold off
grid

figure
    hold on
    plot(protochol)
    plot(squeeze(num_of_LTDs(iter,col,:)))
    title('num of ltds over time')
    legend('protochol', 'ltds')
    hold off
grid

%% watch synapse and spikes
presyn = 5;
postsyn = 15;
iter=1;c=1;
figure
    hold on
        plot(squeeze(spikes(iter,c,presyn,:))*100, '--p')
        plot(squeeze(spikes(iter,c,postsyn,:))*100)
        plot(squeeze(exc_to_exc_weight_matrix(iter,c,:,presyn,postsyn)));
        plot(squeeze(exc_to_exc_weight_matrix(iter,c,:,postsyn,presyn)));
        l1 = num2str(presyn) + ' ' + 'spikes';
        l2 = num2str(postsyn) + ' ' + 'spikes';
        l3 = num2str(presyn) + '->' + num2str(postsyn);
        l4 = num2str(postsyn) + '->' + num2str(presyn);
        legend(l1,  l2, l3, l4)
    hold off
    
grid

%% imagesc of stacked side by side
close all;
weights_stacked_side_by_side = zeros(n_excitatory*n_excitatory,length(tspan));
for t=1:length(tspan)
    weight_at_t = exc_to_exc_weight_matrix(1,1,t,:,:);
    weight_at_t_reshaped = reshape(weight_at_t, n_excitatory*n_excitatory,1);
    weights_stacked_side_by_side(:,t) = weight_at_t_reshaped;
end

figure
    imagesc(weights_stacked_side_by_side);
    title('weights stacked side by side')
grid