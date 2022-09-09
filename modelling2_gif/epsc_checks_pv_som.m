close all
pv_neurons_index=21:23;
iter_to_see = 1; col = 2;

for p=pv_neurons_index
    figure
        hold on
            plot(squeeze(recurrence_inh_pv_epsc_tensor(iter_to_see, col, p, :)))
            
            plot(squeeze(recurrence_inh_som_epsc_tensor(iter_to_see, col, p, :)))
            
            inh_epsc = recurrence_inh_self_column_epsc_tensor(iter_to_see, col, p, :) + recurrence_inh_neighbour_column_epsc_tensor(iter_to_see, col, p, :);
            plot(squeeze(inh_epsc));
        hold off
        title(['pv-n-',num2str(p),'-negative epscs-','som-red-factor-',num2str(som_reduction_factor)])
        legend('pv','som','negative')
        
    grid
end

%% negative current avg
close all;
iter_to_see = 1;
images_path = "D:\on-modified\";
bin_size = 50;
for col=1:n_columns
    figure
        inh_epsc = recurrence_inh_self_column_epsc_tensor(iter_to_see, col, 1:20, :) + recurrence_inh_neighbour_column_epsc_tensor(iter_to_see, col, 1:20, :);
        inh_epsc_squeezed = squeeze(inh_epsc);
        mean_over_all_exc_neurons_inh_epsc = mean(inh_epsc_squeezed, 1); 
        reshaped_inh_epsc = reshape(mean_over_all_exc_neurons_inh_epsc, bin_size, length(mean_over_all_exc_neurons_inh_epsc)/bin_size);
        binned_inh_epsc = mean(reshaped_inh_epsc,1);
        plot(binned_inh_epsc);
        title(['avg-epsc-20-neurons-som-red-fac-',num2str(som_reduction_factor),'-col-',num2str(col)])
        image_name = images_path + "col-" + num2str(col) + "---" + num2str(som_reduction_factor) + "---" + 'bin-size-' + num2str(bin_size) + '--' +  "-neg_epsc.fig"; 
        saveas(gcf, image_name);
    grid
end

%% pv current avg

close all;
iter_to_see = 1;
images_path = "D:\on-modified\";
bin_size = 50;
for col=1:n_columns
    figure
        pv_epsc = recurrence_inh_pv_epsc_tensor(iter_to_see, col, 1:20, :);
        pv_epsc_squeezed = squeeze(pv_epsc);
        mean_over_all_exc_neurons_pv_epsc = mean(pv_epsc_squeezed, 1); 
        reshaped_pv_epsc = reshape(mean_over_all_exc_neurons_pv_epsc, bin_size, length(mean_over_all_exc_neurons_pv_epsc)/bin_size);
        binned_pv_epsc = mean(reshaped_pv_epsc,1);
        plot(binned_pv_epsc);
        title(['pv-epsc-20-neurons-som-red-fac-',num2str(som_reduction_factor),'-col-',num2str(col)])
        image_name = images_path + "col-" + num2str(col) + "---" + num2str(som_reduction_factor) + "---" + 'bin-size-' + num2str(bin_size) + '--' +  "-pv-epsc.fig"; 
        saveas(gcf, image_name);
    grid
end

%% som current average

close all;
iter_to_see = 1;
images_path = "D:\on-modified\";
bin_size = 50;
for col=1:n_columns
    figure
        som_epsc = recurrence_inh_som_epsc_tensor(iter_to_see, col, 1:20, :);
        som_epsc_squeezed = squeeze(som_epsc);
        mean_over_all_exc_neurons_som_epsc = mean(som_epsc_squeezed, 1); 
        reshaped_som_epsc = reshape(mean_over_all_exc_neurons_som_epsc, bin_size, length(mean_over_all_exc_neurons_som_epsc)/bin_size);
        binned_som_epsc = mean(reshaped_som_epsc,1);
        plot(binned_som_epsc);
        title(['som-epsc-20-neurons-som-red-fac-',num2str(som_reduction_factor),'-col-',num2str(col)])
        image_name = images_path + "col-" + num2str(col) + "---" + num2str(som_reduction_factor) + "---" + 'bin-size-' + num2str(bin_size) + '--' +  "-som-epsc.fig"; 
        saveas(gcf, image_name);
    grid
end

%% analyse negative epsc inc/dec

close all;
iter_to_see = 1;
images_path = "D:\on-modified\";
bin_size = 25;
som_on_rec_self_col_inh_tensor = load('batch_1_1_.mat', 'recurrence_inh_self_column_epsc_tensor').recurrence_inh_self_column_epsc_tensor;
som_on_rec_neigh_col_inh_tensor = load('batch_1_1_.mat', 'recurrence_inh_neighbour_column_epsc_tensor').recurrence_inh_neighbour_column_epsc_tensor;

som_off_rec_self_col_inh_tensor = load('batch_1_0_.mat', 'recurrence_inh_self_column_epsc_tensor').recurrence_inh_self_column_epsc_tensor;
som_off_rec_neigh_col_inh_tensor = load('batch_1_0_.mat', 'recurrence_inh_neighbour_column_epsc_tensor').recurrence_inh_neighbour_column_epsc_tensor; 

for col=1:n_columns
    figure

        inh_epsc = som_on_rec_self_col_inh_tensor(iter_to_see, col, 1:20, :) + som_on_rec_neigh_col_inh_tensor(iter_to_see, col, 1:20, :);
        inh_epsc_squeezed = squeeze(inh_epsc);
        mean_over_all_exc_neurons_inh_epsc = mean(inh_epsc_squeezed, 1); 
        reshaped_inh_epsc = reshape(mean_over_all_exc_neurons_inh_epsc, bin_size, length(mean_over_all_exc_neurons_inh_epsc)/bin_size);
        binned_inh_epsc_som_on = mean(reshaped_inh_epsc,1);
        

        inh_epsc = som_off_rec_self_col_inh_tensor(iter_to_see, col, 1:20, :) + som_off_rec_neigh_col_inh_tensor(iter_to_see, col, 1:20, :);
        inh_epsc_squeezed = squeeze(inh_epsc);
        mean_over_all_exc_neurons_inh_epsc = mean(inh_epsc_squeezed, 1); 
        reshaped_inh_epsc = reshape(mean_over_all_exc_neurons_inh_epsc, bin_size, length(mean_over_all_exc_neurons_inh_epsc)/bin_size);
        binned_inh_epsc_som_off = mean(reshaped_inh_epsc,1);
        hold on
            plot(binned_inh_epsc_som_off-binned_inh_epsc_som_on);
            plot( zeros(length(binned_inh_epsc_som_off),1) )
        hold off
        title(['avg-epsc-20-neurons-som-red-fac-',num2str(som_reduction_factor),'-col-',num2str(col)])
        image_name = images_path + "diffference-in-epsc-col-" + num2str(col) + "---" + num2str(som_reduction_factor) + "---" + 'bin-size-' + num2str(bin_size) + '--' +  "-neg_epsc.fig"; 
        saveas(gcf, image_name);
    grid
end

%% in 9 cols - som on vs off case
close all
som_on_file_path = "D:\Downloads\ns-codes\modelling2_gif\batch_1_1_.mat";
som_off_file_path = "D:\Downloads\ns-codes\modelling2_gif\batch_1_0_.mat";

som_on_avg_pv_epsc = load(som_on_file_path,'recurrence_inh_pv_epsc_tensor').recurrence_inh_pv_epsc_tensor;
som_on_avg_som_epsc = load(som_on_file_path,'recurrence_inh_som_epsc_tensor').recurrence_inh_som_epsc_tensor;
som_on_avg_inh_self_epsc = load(som_on_file_path,'recurrence_inh_self_column_epsc_tensor').recurrence_inh_self_column_epsc_tensor;
som_on_avg_inh_across_epsc = load(som_on_file_path,'recurrence_inh_neighbour_column_epsc_tensor').recurrence_inh_neighbour_column_epsc_tensor;

som_off_avg_pv_epsc = load(som_off_file_path,'recurrence_inh_pv_epsc_tensor').recurrence_inh_pv_epsc_tensor;
som_off_avg_som_epsc = load(som_off_file_path,'recurrence_inh_som_epsc_tensor').recurrence_inh_som_epsc_tensor;
som_off_avg_inh_self_epsc = load(som_off_file_path,'recurrence_inh_self_column_epsc_tensor').recurrence_inh_self_column_epsc_tensor;
som_off_avg_inh_across_epsc = load(som_off_file_path,'recurrence_inh_neighbour_column_epsc_tensor').recurrence_inh_neighbour_column_epsc_tensor;

n_columns = load(som_off_file_path,'n_columns').n_columns;

iter = 1;

% figure
%     hold on
%            plot(squeeze(mean(som_on_avg_pv_epsc(iter,col,1:20,:), 3)))
%            plot(squeeze(mean(som_off_avg_pv_epsc(iter,col,1:20,:), 3)))
%     hold off
%     legend('on', 'off')
%     title('pv epsc - on vs off')
% grid
% 
% figure
%     hold on
%          plot(squeeze(mean(som_on_avg_som_epsc(iter,col,1:20,:), 3)))
%          plot(squeeze(mean(som_off_avg_som_epsc(iter,col,1:20,:), 3)))
%     hold off
%     legend('on', 'off')
%     title('som epsc - on vs off')
% grid
% 
% figure
%     hold on
%         plot(squeeze(mean(som_on_avg_inh_self_epsc(iter,col,1:20,:), 3)) + squeeze(mean(som_on_avg_inh_across_epsc(iter,col,1:20,:), 3)))
%         plot(squeeze(mean(som_off_avg_inh_self_epsc(iter,col,1:20,:), 3)) + squeeze(mean(som_off_avg_inh_across_epsc(iter,col,1:20,:), 3)))
%     hold off
%     legend('on', 'off')
%     title('inh epsc - on vs off')
% grid

images_path = "D:\rough-off-on";

for col=1:n_columns
    figure
        on_inh_epsc = squeeze(mean(som_on_avg_inh_self_epsc(iter,col,1:20,:), 3)) + squeeze(mean(som_on_avg_inh_across_epsc(iter,col,1:20,:), 3));
        off_inh_epsc = squeeze(mean(som_off_avg_inh_self_epsc(iter,col,1:20,:), 3)) + squeeze(mean(som_off_avg_inh_across_epsc(iter,col,1:20,:), 3));
        hold on
            plot(get_binned_vec(off_inh_epsc - on_inh_epsc, 25))
            yline(0)
        hold off
        title(['off - on(should be more negative)-col-6-fac', num2str(col)])
        image_path = strcat(images_path,'\',num2str(col),'-off-on-6-fac');
        saveas(gcf, image_path);
    grid
end

%% spikes for 2nd stimulus in on and off - (off - on) (should be more positive)

close all
som_on_file_path = "D:\Downloads\ns-codes\modelling2_gif\batch_1_1_.mat";
som_off_file_path = "D:\Downloads\ns-codes\modelling2_gif\batch_1_0_.mat";

som_on_spikes = load(som_on_file_path,'spikes').spikes;
som_off_spikes = load(som_off_file_path,'spikes').spikes;

second_stim_start = 110:170:790;

n_columns = load(som_off_file_path,'n_columns').n_columns;
n_tokens = load(som_off_file_path,'n_tokens').n_tokens;
iter = 1;

som_on_spikes_stimB = zeros(n_columns, n_tokens);
som_off_spikes_stimB = zeros(n_columns, n_tokens);

for c=1:n_columns
    on_col_mean_spikes = squeeze(mean(som_on_spikes(iter,c,:,:), 3));
    off_col_mean_spikes = squeeze(mean(som_off_spikes(iter,c,:,:), 3));
    for tok=1:n_tokens
        som_on_spikes_stimB(c,tok) = sum(on_col_mean_spikes(second_stim_start(tok):second_stim_start(tok) + 50));
        som_off_spikes_stimB(c,tok) = sum(off_col_mean_spikes(second_stim_start(tok):second_stim_start(tok) + 50));
    end
end

images_path = "D:\rough-spikes-off-on";
for c=1:n_columns
    figure
        hold on
            stem(som_off_spikes_stimB(c,:) - som_on_spikes_stimB(c,:));
            yline(0)
        hold off
        
        title(['spikes off vs on-', num2str(c)])
        image = strcat(images_path,'\','col-',num2str(c),'-spikes-off-on(should be more +ve)');
        saveas(gcf, image);
    grid
end