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

close all;
post_syn = 15;
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
