%% basic vars
batches = 50;
data_path = "D:\28june-test";
images_path = "D:\28june-test-analysis\";
n_columns = 5;
n_excitatory=20; n_pv = 3; n_som  = 2;
n_neurons = n_excitatory + n_pv + n_som;
iter=1;

%% rate of exc, pv, som
close all;
psth_all = zeros(batches, n_columns, n_neurons);

psth_avg_exc = zeros(n_columns,n_excitatory, batches);
psth_avg_pv = zeros(n_columns,n_pv, batches);
psth_avg_som = zeros(n_columns,n_som, batches);

psth_a = zeros(batches, n_columns, n_neurons);
psth_b = zeros(batches, n_columns, n_neurons);

for b=1:batches
    fprintf("\n batch num %d \n",b);

    batch_file_name = data_path + "\batch_" + num2str(b) + ".mat";
    spikes = load(batch_file_name,"spikes").spikes;
    lamda = load(batch_file_name,"lamda").lamda;
    t_simulate = load(batch_file_name,"t_simulate").t_simulate;
    
    for c=1:n_columns
        t_a = 0;
        t_b = 0;
        for t=1:t_simulate
            if lamda(1,4,1,t) == 300
                t_a = t_a + 1;
                for N=1:25
                    psth_a(b,c,N) = psth_a(b,c,N) + spikes(1,c,N,t);
                end
            end

            if lamda(1,4,1,t) == 50
                t_b = t_b + 1;
                for N=1:25
                    psth_b(b,c,N) = psth_b(b,c,N) + spikes(1,c,N,t);
                end
            end
        end

        psth_a(b,c,:) = psth_a(b,c,:)/(t_a*0.001);
        psth_b(b,c,:) = psth_b(b,c,:)/(t_b*0.001);
             
        

        for N=1:25
            psth_all(b,c,N) = sum(spikes(iter,c,N,:))/(1700 * 0.001);
        end

        for nn=1:n_excitatory
            psth_avg_exc(c,nn,b) = sum(spikes(iter,c,nn,:))/(1700 * 0.001);
        end

        for nn=n_excitatory+1:n_excitatory+n_pv
            psth_avg_pv(c,nn-n_excitatory,b) = sum(spikes(iter,c,nn,:))/(1700 * 0.001);
        end

        for nn=n_excitatory+n_pv+1:n_excitatory+n_pv+n_som
            psth_avg_som(c,nn-(n_excitatory+n_pv),b) = sum(spikes(iter,c,nn,:))/(1700 * 0.001);
        end
    end
end

for c=1:n_columns
    figure
        plot(transpose(squeeze(psth_avg_exc(c,:,:))));
        title(['exc-col ',num2str(c)])
        image_name = images_path + "col-" + num2str(c) + "psth_exc_imagesc.fig"; 
        saveas(gcf, image_name);
    grid

    figure
        plot(transpose(squeeze(psth_avg_pv(c,:,:))));
        title(['pv-col ',num2str(c)])
        image_name = images_path + "col-" + num2str(c) + "psth_pv_imagesc.fig"; 
        saveas(gcf, image_name);
    grid

    figure
        plot(transpose(squeeze(psth_avg_som(c,:,:))));
        title(['som-col ',num2str(c)])
        image_name = images_path + "col-" + num2str(c) + "psth_som_imagesc.fig"; 
        saveas(gcf, image_name);
    grid

    fprintf("\n writing col %d \n",c)
end

%% rate for a stimulus, b stimulus
close all;
n_columns = 5;
for c=1:n_columns
    col_psth = squeeze(psth_all(:,c,1:20));
    figure
         plot(col_psth)
         title(['all neurons psth batch wise avg-col-', num2str(c)])
         image_name = images_path + "psth-c-" + num2str(c) + ".fig"; 
         saveas(gcf, image_name);
    grid
end

for c=1:5
    col_psth_a = squeeze(psth_a(:,c,1:20));
    col_psth_a_mean = mean(col_psth_a,2);

    col_psth_b = squeeze(psth_b(:,c,1:20));
    col_psth_b_mean = mean(col_psth_b,2);

    figure
        hold on
            plot(col_psth_a_mean)
            plot(col_psth_b_mean)
            legend('psth for 1st stim','psth for 2nd stim')
            title(['psth - a,b col ',num2str(c)])
            image_name = images_path + "psth-a-b-col-" + num2str(c) + ".fig"; 
            saveas(gcf, image_name);
        hold off
    grid
end

%% rate in bins
close all;
bin_size = 100;


t_simulate = load(strcat(data_path, '\', 'batch_1.mat'), "t_simulate").t_simulate;
n_columns = load(strcat(data_path, '\', 'batch_1.mat'), "n_columns").n_columns;
n_total_neurons = load(strcat(data_path, '\', 'batch_1.mat'), "n_total_neurons").n_total_neurons;
n_bins = t_simulate/bin_size;
psth_binned = zeros(n_columns, n_total_neurons, batches*n_bins);
iter=1;
for b=1:batches
    batch_file_name = data_path + "\batch_" + num2str(b) + ".mat";
    spikes = load(batch_file_name,"spikes").spikes;
    
    for c=1:n_columns
        for n=1:n_total_neurons
            spikes_cn = squeeze(spikes(iter,c,n,1:t_simulate));
            spikes_cn_reshaped = reshape(spikes_cn,   n_bins, bin_size);
            spikes_cn_reshaped_avg = sum(spikes_cn_reshaped, 2);
            psth_binned(c,n,(b-1)*n_bins + 1:(b-1)*n_bins + n_bins) = reshape(spikes_cn_reshaped_avg,  1,1,n_bins);
        end
    end
end


for c=1:n_columns
    rate_c = squeeze(psth_binned(c,:,:));
    rate_c_avg = sum(rate_c, 1);
    figure
        plot(rate_c_avg);
        title(['rate 100ms binned-col-', num2str(c)])
        image_name = images_path + "rate-100ms-avg-col-" + num2str(c) + ".fig"; 
        saveas(gcf, image_name);

    grid

    figure
        plot(rate_c.');
        title(['rate 100ms binned-col-', num2str(c)])
        image_name = images_path + "rate-100ms-all-neurons-col-" + num2str(c) + ".fig"; 
        saveas(gcf, image_name);
    grid
end

%% for each token AB 
close all
t_simulate = load(strcat(data_path, '\', 'batch_1.mat'), "t_simulate").t_simulate;
n_columns = load(strcat(data_path, '\', 'batch_1.mat'), "n_columns").n_columns;
n_exc = 20;
iter=1;
rate_for_ab = zeros(n_columns, n_exc, (t_simulate/170)*batches);
for b=1:batches
    fprintf("\n batch is %d \n",b)
    batch_file_name = data_path + "\batch_" + num2str(b) + ".mat";
    spikes = load(batch_file_name,"spikes").spikes;

    for c=1:n_columns
        for n=1:n_exc
            spikes_cn = squeeze(spikes(iter,c,n,1:t_simulate));
            spikes_cn_reshaped = reshape(spikes_cn,   10, 170);
            spikes_cn_reshaped_avg = sum(spikes_cn_reshaped, 2)/(170*0.001);
            rate_for_ab(c,n,(b-1)*10 + 1:(b-1)*10 + 10) = reshape(spikes_cn_reshaped_avg,  1,1,10);
        end
    end
end



for c=1:n_columns
    figure
        plot(transpose(squeeze(rate_for_ab(c,:,:))))
        title(['rate-170ms-c-',num2str(c)])
        image_name = images_path + "-rate-170ms-bin-c-" + num2str(c) + ".fig"; 
        saveas(gcf, image_name);
    grid
end

for c=1:n_columns
    figure
        mean_rate_over_neurons = mean(squeeze(rate_for_ab(c,:,:)), 1);
        plot(mean_rate_over_neurons)
        title(['rate-170ms-c-',num2str(c)])
        image_name = images_path + "-rate-avg-over-neurons-170ms-bin-c-" + num2str(c) + ".fig"; 
        saveas(gcf, image_name);
    grid
end

%% weights
close all;
tspan = load(strcat(data_path, '\', 'batch_1.mat'), "tspan").tspan;
num_network_neurons = 125;
batch_avg_network_weights = zeros(batches, num_network_neurons, num_network_neurons);
n_exc = 20;



iter=1;

for b=1:batches
    fprintf("\n batch is %d \n",b)
    batch_file_name = data_path + "\batch_" + num2str(b) + ".mat";
    network_weight_matrix = load(batch_file_name, "network_weight_matrix").network_weight_matrix;
    for n1=1:num_network_neurons
        for n2=1:num_network_neurons
             batch_avg_network_weights(b,n1,n2) = mean(squeeze(network_weight_matrix(iter,:,n1,n2)));
%              all_1ms_network_weights((b-1)*length(tspan) + 1:(b-1)*length(tspan) + length(tspan)   ,n1,n2) =  network_weight_matrix(iter,:,n1,n2);
        end
    end


end

% plot weight matrix
n_columns = 5;
for c1=1:n_columns
    for c2=1:n_columns
        if c1 - c2 == 0
            % within column matrix
            within_column_matrix = batch_avg_network_weights(:,(c1-1)*n_total_neurons + 1:(c1-1)*n_total_neurons + n_exc, (c1-1)*n_total_neurons + 1:(c1-1)*n_total_neurons + n_exc);
            fprintf("\n %d %d %d %d \n ",(c1-1)*n_total_neurons + 1,(c1-1)*n_total_neurons + n_exc, (c1-1)*n_total_neurons + 1,(c1-1)*n_total_neurons + n_exc)
            reshaped_within_column_matrix = reshape(within_column_matrix, batches, n_exc*n_exc);
            figure
                plot(reshaped_within_column_matrix)
                title(['weights batches wise', num2str(c1)])
                image_name = images_path + "w_batch_avg-c-" + num2str(c1) + ".fig"; 
                saveas(gcf, image_name);

            grid
        elseif abs(c1 - c2) == 1 || abs(c1 - c2) == 2
            % c1 pre ,c2 post
            across_column_matrix = batch_avg_network_weights(:,(c1-1)*n_total_neurons + 1:(c1-1)*n_total_neurons + n_exc, (c2-1)*n_total_neurons + 1:(c2-1)*n_total_neurons + n_exc);
            fprintf("\n %d %d %d %d \n ",(c1-1)*n_total_neurons + 1,(c1-1)*n_total_neurons + n_exc, (c2-1)*n_total_neurons + 1,(c2-1)*n_total_neurons + n_exc)
            reshaped_across_column_matrix = reshape(across_column_matrix, batches, n_exc*n_exc);
            figure
                plot(reshaped_across_column_matrix)
                title(['weights batches wise', num2str(c1), '-',num2str(c2)])
                image_name = images_path + "w_batch_avg-c1-" + num2str(c1) + "-c2-" + num2str(c2)+ ".fig"; 
                saveas(gcf, image_name);
            grid

        else
            continue
        end
    end
end


