%% check non column
figure
    nt = 65;
    hold on
        plot(thalamus_to_all_l4_epsc(nt,:))
        stem(spikes(nt,:))
    hold off
grid

nt = 20;
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
    plot(total_input_current(nt,:),'LineWidth',2,'LineStyle','-')
    
 
    legend('spikes', 'voltage', 'theta','stimulus', 'thalamus-epsc', 'exc','pv','som','total epsc')
    hold off
grid

figure
    stem(shuffled_neuron_types)
grid

figure
    plot(squeeze(weight_matrix(70,42,:)))
grid
figure
    plot(squeeze(weight_matrix(20,70,:)))
    hold on
    plot(squeeze(weight_matrix(70,20,:)))
    stem(3.1*spikes(20,:), '-o')
    stem(3*spikes(70,:), '-+')
    legend('20-70','70-20', '20','70')
grid

figure
    plot(voltages.')
grid

figure
    plot(squeeze(thalamic_epsc(:,10,:)).')
    legend('1','2','3','4','5','6','7','8')
grid

%% 
close all
a = 227;
b = 232 ;

figure
hold on
    plot(squeeze(network_weight_matrix(1,:,a,b)))
    plot(squeeze(network_weight_matrix(1,:,b,a)))
    l1=find(squeeze(spikes(1,10,2,:))==1);
    l2=find(squeeze(spikes(1,10,7,:))==1);
    
    plot(l1,5.5*ones(1,length(l1)),'*')
    plot(l2,5.3*ones(1,length(l2)),'*')
    legend('a to b', 'b to a', 'spikes a', 'spikes b')
hold off
grid

%% last from each batch
dir = "D:\som-on-21cols-10-12-BA-data";
w_a_to_b = [];
w_b_to_a = [];
for i=1:500
    disp(i)
    file = strcat(dir,'\batch_', num2str(i),'.mat');
    w_matrix = load(file, 'network_weight_matrix').network_weight_matrix;
    w_a_to_b = [w_a_to_b, w_matrix(1,end,a,b)];
    w_b_to_a = [w_b_to_a, w_matrix(1,end,b,a)];
end

figure
    hold on
        plot(w_a_to_b)
        plot(w_b_to_a)
        legend('a to b', 'b to a')
    hold off
grid