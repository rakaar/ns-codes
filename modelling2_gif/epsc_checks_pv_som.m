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
images_path = "D:\20-ration\";
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
images_path = "D:\reduce-imbalance\";
bin_size = 100;
for col=1:n_columns
    figure
        pv_epsc = recurrence_inh_pv_epsc_tensor(iter_to_see, col, 1:20, :);
        pv_epsc_squeezed = squeeze(pv_epsc);
        mean_over_all_exc_neurons_pv_epsc = mean(pv_epsc_squeezed, 1); 
        reshaped_pv_epsc = reshape(mean_over_all_exc_neurons_pv_epsc, bin_size, length(mean_over_all_exc_neurons_pv_epsc)/bin_size);
        binned_pv_epsc = mean(reshaped_pv_epsc,1);
        plot(binned_pv_epsc);
        title(['old-pv-epsc-20-neurons-som-red-fac-',num2str(som_reduction_factor),'-col-',num2str(col)])
        image_name = images_path + "-new-col-" + num2str(col) + "---" + num2str(som_reduction_factor) + "---" + 'bin-size-' + num2str(bin_size) + '--' +  "-pv-epsc.fig"; 
        saveas(gcf, image_name);
    grid
end
