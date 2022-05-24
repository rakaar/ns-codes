close all;
%% presyn fix
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