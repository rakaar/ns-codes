
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
