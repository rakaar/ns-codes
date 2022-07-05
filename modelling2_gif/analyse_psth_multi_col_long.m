
clear all; close all;


n_columns = 5;
n_excitatory=20; n_pv = 3; n_som  = 2;
n_neurons = n_excitatory + n_pv + n_som;
batches = 100;
iter=1;

batch_data_path = "D:\5_multi_col_across_ap";
images_path = "D:\5_multi_col_across_ap_images\";


psth_all = zeros(batches, n_columns, n_neurons);

psth_avg_exc = zeros(n_columns,n_excitatory, batches);
psth_avg_pv = zeros(n_columns,n_pv, batches);
psth_avg_som = zeros(n_columns,n_som, batches);

psth_a = zeros(batches, n_columns, n_neurons);
psth_b = zeros(batches, n_columns, n_neurons);

for b=1:batches
    fprintf("\n batch num %d \n",b);

    batch_file_name = batch_data_path + "\batch_" + num2str(b) + ".mat";
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

%% rate for non long code
% close all
% bin_size = 20;
% iter = 1;
% rate = zeros(n_columns, n_excitatory, (length(tspan)-1)/bin_size);
% for c=1:n_columns
%     for n=1:n_excitatory
%         spikes_n = reshape(spikes(iter,c,n,1:length(tspan)-1),  bin_size,(length(tspan)-1)/bin_size);
%         spikes_avg = mean(spikes_n(:, 1:(length(tspan)-1)/bin_size))/0.001;
%         rate(c,n,:) = spikes_avg;
%     end
% end
% 
% for c=1:n_columns
%     figure
%         plot(transpose(squeeze(rate(c,:,:))))
%         title(['rates-col-',num2str(c)])
%     grid
% end


%% analyse spikes binned
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

%% binned rate

bin_size = 100;
batch_data_path = "D:\5_multi_col_across_ap";
images_path = "D:\5_multi_col_across_ap_images\";

t_simulate = load("D:\5_multi_col_across_ap\batch_1.mat", "t_simulate").t_simulate;
n_columns = load("D:\5_multi_col_across_ap\batch_1.mat", "n_columns").n_columns;
n_total_neurons = load("D:\5_multi_col_across_ap\batch_1.mat", "n_total_neurons").n_total_neurons;
n_bins = t_simulate/bin_size;
psth_binned = zeros(n_columns, n_total_neurons, batches*n_bins);
iter=1;
for b=1:batches
    batch_file_name = batch_data_path + "\batch_" + num2str(b) + ".mat";
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

%% for one AB 
close all
batch_data_path = "D:\5_multi_col_across_ap";
images_path = "D:\5_multi_col_across_ap_images\";

t_simulate = load("D:\5_multi_col_across_ap\batch_1.mat", "t_simulate").t_simulate;
n_columns = load("D:\5_multi_col_across_ap\batch_1.mat", "n_columns").n_columns;
n_exc = 20;
batches = 100;
iter=1;
rate_for_ab = zeros(n_columns, n_exc, (t_simulate/170)*batches);
for b=1:batches
    fprintf("\n batch is %d \n",b)
    batch_file_name = batch_data_path + "\batch_" + num2str(b) + ".mat";
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
        title(['rate-ab-c-',num2str(c)])
        image_name = images_path + "-rate-ab-c-" + num2str(c) + ".fig"; 
        saveas(gcf, image_name);
    grid
end


for c=1:n_columns
    figure
        plot(transpose(squeeze(rate_for_ab(c,:,:))))
        title(['rate-ab-c-',num2str(c)])
        image_name = images_path + "-rate-ab-c-" + num2str(c) + ".fig"; 
        saveas(gcf, image_name);
    grid
end

for c=1:n_columns
    figure
        mean_rate_over_neurons = mean(squeeze(rate_for_ab(c,:,:)), 1);
        plot(mean_rate_over_neurons)
        title(['rate-ab-c-',num2str(c)])
        image_name = images_path + "-rate-avg-over-neurons-ab-c-" + num2str(c) + ".fig"; 
        saveas(gcf, image_name);
    grid
end


