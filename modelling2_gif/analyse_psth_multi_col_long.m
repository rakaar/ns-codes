clear all; close all;


n_columns = 5;
n_excitatory=20; n_pv = 3; n_som  = 2;
batches = 200;
iter=1;

batch_data_path = "D:\7_multi_col_big_clip_range";
images_path = "D:\7_multi_col_big_clip_images\";




psth_avg_exc = zeros(n_columns,n_excitatory, batches);
psth_avg_pv = zeros(n_columns,n_pv, batches);
psth_avg_som = zeros(n_columns,n_som, batches);

for b=1:batches
    fprintf("\n batch num %d \n",b);

    batch_file_name = batch_data_path + "\batch_" + num2str(b) + ".mat";
    spikes = load(batch_file_name,"spikes").spikes;
    for c=1:n_columns
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
        imagesc(squeeze(psth_avg_exc(c,:,:)));
        title(['col ',num2str(c)])
        image_name = images_path + "col-" + num2str(c) + "psth_exc_imagesc.fig"; 
        saveas(gcf, image_name);
    grid

    figure
        imagesc(squeeze(psth_avg_pv(c,:,:)));
        title(['col ',num2str(c)])
        image_name = images_path + "col-" + num2str(c) + "psth_pv_imagesc.fig"; 
        saveas(gcf, image_name);
    grid

    figure
        imagesc(squeeze(psth_avg_som(c,:,:)));
        title(['col ',num2str(c)])
        image_name = images_path + "col-" + num2str(c) + "psth_som_imagesc.fig"; 
        saveas(gcf, image_name);
    grid

    fprintf("\n writing col %d \n",c)
end