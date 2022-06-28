clear all; close all;


n_columns = 5;
n_excitatory=20; n_pv = 3; n_som  = 2;
n_neurons = n_excitatory + n_pv + n_som;
batches = 50;
iter=1;

batch_data_path = "D:\1_multi_col_across_plastic";
images_path = "D:\1_multi_col_across_col_plastic_images\";


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
            if lamda(1,c+2,1,t) == 300
                t_a = t_a + 1;
                for N=1:25
                    psth_a(b,c,N) = psth_a(b,c,N) + spikes(1,c,N,t);
                end
            end

            if lamda(1,c+2,1,t) == 50
                t_b = t_b + 1;
                for N=1:25
                    psth_b(b,c,N) = psth_b(b,c,N) + spikes(1,c,N,t);
                end
            end
        end

        psth_a(b,c,:) = psth_a(b,c,:)/(t_a*0.001);
        psth_b(b,c,:) = psth_b(b,c,:)/(t_b*0.001);
             
        

        for N=1:25
            psth_all(b,c,N) = mean(spikes(iter,c,N,:))/0.001;
        end

        for nn=1:n_excitatory
            psth_avg_exc(c,nn,b) = mean(spikes(iter,c,nn,:))/0.001;
        end

        for nn=n_excitatory+1:n_excitatory+n_pv
            psth_avg_pv(c,nn-n_excitatory,b) = mean(spikes(iter,c,nn,:))/0.001;
        end

        for nn=n_excitatory+n_pv+1:n_excitatory+n_pv+n_som
            psth_avg_som(c,nn-(n_excitatory+n_pv),b) = mean(spikes(iter,c,nn,:))/0.001;
        end
    end
end

for c=1:n_columns
    figure
        plot(transpose(squeeze(psth_avg_exc(c,:,:))));
        title(['col ',num2str(c)])
        image_name = images_path + "col-" + num2str(c) + "psth_exc_imagesc.fig"; 
        saveas(gcf, image_name);
    grid

    figure
        plot(transpose(squeeze(psth_avg_pv(c,:,:))));
        title(['col ',num2str(c)])
        image_name = images_path + "col-" + num2str(c) + "psth_pv_imagesc.fig"; 
        saveas(gcf, image_name);
    grid

    figure
        plot(transpose(squeeze(psth_avg_som(c,:,:))));
        title(['col ',num2str(c)])
        image_name = images_path + "col-" + num2str(c) + "psth_som_imagesc.fig"; 
        saveas(gcf, image_name);
    grid

    fprintf("\n writing col %d \n",c)
end

%% rate for non long code
close all
iter = 1;
rate = zeros(n_columns, n_excitatory, (length(tspan)-1)/10);
for c=1:n_columns
    for n=1:n_excitatory
        spikes_n = reshape(spikes(iter,c,n,1:length(tspan)-1),  10,(length(tspan)-1)/10);
        spikes_avg = mean(spikes_n(:, 1:(length(tspan)-1)/10))/0.001;
        rate(c,n,:) = spikes_avg;
    end
end

for c=1:n_columns
    figure
        imagesc(squeeze(rate(c,:,:)))
        title(['rates-col-',num2str(c)])
    grid
end

