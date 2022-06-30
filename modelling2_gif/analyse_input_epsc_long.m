%% input recurrent current
close all;
n_columns = 5;
batches = 30;
n_excitatory = 20;
iter=1;
batch_data_path = "D:\2_multi_col_across_plastic";
images_path = "D:\2_multi_col_images\";

% total_input_epsc = thalamic_epsc_tensor ...
%                     + recurrence_exc_self_column_epsc_tensor + recurrence_inh_self_column_epsc_tensor ...
%                     + recurrence_exc_neighbour_column_epsc_tensor + recurrence_inh_neighbour_column_epsc_tensor ;

thalamic_input_epsc_batch_avg = zeros(batches,n_columns,n_excitatory);
recurrence_exc_self_column_epsc_tensor_batch_avg = zeros(batches,n_columns,n_excitatory);
recurrence_inh_self_column_epsc_tensor_batch_avg = zeros(batches,n_columns,n_excitatory);
recurrence_exc_neighbour_column_epsc_tensor_batch_avg = zeros(batches,n_columns,n_excitatory);
recurrence_inh_neighbour_column_epsc_tensor_batch_avg = zeros(batches,n_columns,n_excitatory);

for b=1:batches
    batch_file_name = batch_data_path + "\batch_" + num2str(b) + ".mat";
    fprintf("\n batch num is %d \n", b)
    thalamic_epsc_tensor = load(batch_file_name,"thalamic_epsc_tensor").thalamic_epsc_tensor;
    recurrence_exc_self_column_epsc_tensor = load(batch_file_name,"recurrence_exc_self_column_epsc_tensor").recurrence_exc_self_column_epsc_tensor;
    recurrence_inh_self_column_epsc_tensor = load(batch_file_name,"recurrence_inh_self_column_epsc_tensor").recurrence_inh_self_column_epsc_tensor;
    recurrence_exc_neighbour_column_epsc_tensor = load(batch_file_name,"recurrence_exc_neighbour_column_epsc_tensor").recurrence_exc_neighbour_column_epsc_tensor;
    recurrence_inh_neighbour_column_epsc_tensor = load(batch_file_name,"recurrence_inh_neighbour_column_epsc_tensor").recurrence_inh_neighbour_column_epsc_tensor;
    
    for c=1:n_columns
        for n=1:n_excitatory
            thalamic_input_epsc_batch_avg(b,c,n) = mean(squeeze(thalamic_epsc_tensor(iter,c,n,:)));
            recurrence_exc_self_column_epsc_tensor_batch_avg(b,c,n) = mean(squeeze(recurrence_exc_self_column_epsc_tensor(iter,c,n,:)));
            recurrence_inh_self_column_epsc_tensor_batch_avg(b,c,n) = mean(squeeze(recurrence_inh_self_column_epsc_tensor(iter,c,n,:)));
            recurrence_exc_neighbour_column_epsc_tensor_batch_avg(b,c,n) = mean(squeeze(recurrence_exc_neighbour_column_epsc_tensor(iter,c,n,:)));
            recurrence_inh_neighbour_column_epsc_tensor_batch_avg(b,c,n) = mean(squeeze(recurrence_inh_neighbour_column_epsc_tensor(iter,c,n,:)));
        end
    end

end

% for c=1:n_columns
%         figure
%            thalamic_epsc_col = squeeze(thalamic_input_epsc_batch_avg(:,c,:)); 
%            imagesc(thalamic_epsc_col)
%            title(['thalamic epsc-col-',num2str(c)])
%             image_name = images_path + "epsc-thalamic-col-" + num2str(c) + ".fig"; 
%             saveas(gcf, image_name);
%         grid
% 
%         figure
%                recurrence_exc_self_column_epsc_tensor_col = squeeze(recurrence_exc_self_column_epsc_tensor_batch_avg(:,c,:)); 
%                imagesc(recurrence_exc_self_column_epsc_tensor_col)
%                title(['recurrence exc self-col-',num2str(c)])
%                image_name = images_path + "epsc-rec-exc-self-col-" + num2str(c) + ".fig"; 
%                 saveas(gcf, image_name);
%         grid
%         
%         figure
%                recurrence_inh_self_column_epsc_tensor_col = squeeze(recurrence_inh_self_column_epsc_tensor_batch_avg(:,c,:)); 
%                imagesc(recurrence_inh_self_column_epsc_tensor_col)
%                title(['recurrence inh self-col-',num2str(c)])
%                image_name = images_path + "epsc-rec-inh-self-col-" + num2str(c) + ".fig"; 
%                 saveas(gcf, image_name);
%         grid
%         
%         figure
%                recurrence_exc_neighbour_column_epsc_tensor_col = squeeze(recurrence_exc_neighbour_column_epsc_tensor_batch_avg(:,c,:)); 
%                imagesc(recurrence_exc_neighbour_column_epsc_tensor_col)
%                title(['recurrence-exc-neigbour-col-',num2str(c)])
%                image_name = images_path + "epsc-rec-exc-neigh-col-" + num2str(c) + ".fig"; 
%                saveas(gcf, image_name);
%         grid
%         
%         figure
%                recurrence_inh_neighbour_column_epsc_tensor_col = squeeze(recurrence_inh_neighbour_column_epsc_tensor_batch_avg(:,c,:)); 
%                imagesc(recurrence_inh_neighbour_column_epsc_tensor_col)
%                title(['recurrence-inh-neighbour-col-',num2str(c)])
%                image_name = images_path + "epsc-rec-inh-neigh-col-" + num2str(c) + ".fig"; 
%                 saveas(gcf, image_name);
%         grid
% 
% 
% end

for c=1:n_columns
    thalamic_epsc_col = squeeze(thalamic_input_epsc_batch_avg(:,c,:));
    thalamic_neurons_avg = mean(thalamic_epsc_col, 2);

    recurrence_exc = recurrence_exc_self_column_epsc_tensor_batch_avg + recurrence_exc_neighbour_column_epsc_tensor_batch_avg;
    recurrence_exc_col = squeeze(recurrence_exc(:,c,:));
    recurrence_exc_avg  = mean(recurrence_exc_col,2);

    recurrence_inh = recurrence_inh_self_column_epsc_tensor_batch_avg + recurrence_inh_neighbour_column_epsc_tensor_batch_avg;
    recurrence_inh_col = squeeze(recurrence_inh(:,c,:));
    recurrence_inh_avg = mean(recurrence_inh_col, 2);

    figure
        hold on
            plot(thalamic_neurons_avg);
            plot(recurrence_exc_avg);
            plot(recurrence_inh_avg);
        hold off
        title(['epscs-col-', num2str(c)])
        legend('feedforward', 'recurrence exc','recurrence inh' )
    grid
end