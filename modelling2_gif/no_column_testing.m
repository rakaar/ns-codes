
 exc_index = find(shuffled_neuron_types == 0);
nt = exc_index(randi(length(exc_index)));
figure
    total_input_current = thalamus_to_all_l4_epsc + rec_exc_epsc + rec_pv_epsc + rec_som_epsc + background_epsc;
    hold on
    stem(10*spikes(nt,:))
    plot(voltages(nt,:))
    plot(thetas(nt,:))
    plot(-3*squeeze(thalamic_cols_spike_rates(1,1,:)))

    plot(thalamus_to_all_l4_epsc(nt,:), 'LineWidth',2)
    plot(rec_exc_epsc(nt,:),'LineWidth',2)
    plot(rec_pv_epsc(nt,:),'LineWidth',2)
    plot(rec_som_epsc(nt,:),'LineWidth',2)
    plot(total_input_current(nt,:),'LineWidth',4,'LineStyle','--')
    
 
    legend('spikes', 'voltage', 'theta','stimulus', 'thalamus-epsc', 'exc','pv','som','total epsc')
    hold off
    title(num2str(nt))
grid

%% number of spikes for A and B
a_rate = zeros(80,4);
b_rate = zeros(80,4);

for e=1:length(exc_index) 
    k_counter=1;
    for k=11:170:521
        a_rate(e,k_counter) = sum(spikes(exc_index(e), k:k+49));
        k_counter = k_counter + 1;
    end

    k_counter=1;
    for k=111:170:621
        b_rate(e,k_counter) = sum(spikes(exc_index(e), k:k+49));
        k_counter = k_counter + 1;
    end
end

figure
    imagesc(a_rate)
    title('a counter')
grid

figure
    imagesc(b_rate)
    title('b counter')
grid

figure
    a_mean = mean(a_rate,1);
    b_mean = mean(b_rate,1);

    hold on
        plot(a_mean)
        plot(b_mean)
    hold off
    legend('a mean', 'b mean')
grid
% figure
%     hold on
%         plot(thetas.')
% %         xline(token_start_times(1), 'LineWidth',0.5)
% %         xline(token_start_times(2), 'LineWidth',0.5)
% %         xline(token_start_times(3), 'LineWidth',0.5)
%     hold off
% grid

% figure
%     imagesc(spikes(find(shuffled_neuron_types == 2),:))
% grid
% figure
%     imagesc(rec_exc_epsc(find(shuffled_neuron_types == 0),:))
% grid

% figure
%     plot(squeeze(g_t_test(87,:,:)).')
% grid

% figure
%     imagesc(rec_som_epsc(find(shuffled_neuron_types == 0), :));
% grid
% 
% nt = exc_index(randi(length(exc_index)));
% figure
%     imagesc(squeeze(weight_matrix(find(shuffled_neuron_types == 1),nt, :)))
% grid

% x1 = load( 'som_epsc_1.mat')
% x1 = x1.rec_som_epsc;
% x0 = load('som_epsc_0.mat').rec_som_epsc;
% xd = x1 - x0;
% imagesc(xd)
% plot(xd.')


%% 
% s = load('batch_3.mat','shuffled_neuron_types').shuffled_neuron_types;
% w1 = load('batch_1.mat', 'weight_matrix').weight_matrix;
% w2 = load('batch_2.mat', 'weight_matrix').weight_matrix;
% w3 = load('batch_3.mat', 'weight_matrix').weight_matrix;
% 
% 
% figure
% hold on
% plot(squeeze(w1(70,87,:)))
% plot(squeeze(w2(70,93,:)))
% plot(squeeze(w3(70,93,:)))
% hold off
% grid