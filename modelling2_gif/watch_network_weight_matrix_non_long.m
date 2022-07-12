close all;
n_columns = 5;
n_exc = 20;
n_total = 25;

for c1=1:n_columns
    for c2=1:n_columns

        pre_start_exc = (c1-1)*n_total + 1;
        pre_end_exc = (c1-1)*n_total + n_exc;
        
        post_start_exc = (c2-1)*n_total + 1;
        post_end_exc = (c2-1)*n_total + n_exc;

        if abs(c1 - c2) == 0 || abs(c1 - c2) == 1 || abs(c1 - c2) == 2
            % within column 
            weight_matrix = network_weight_matrix(1, :,pre_start_exc:pre_end_exc, post_start_exc:post_end_exc);
            squeezed_weight_matrix = squeeze(weight_matrix);
            time_period = size(squeezed_weight_matrix, 1);
            reshaped_weight_matrix = reshape(squeezed_weight_matrix, time_period, n_exc*n_exc);
            
            figure
                plot(reshaped_weight_matrix);
                title(['weights-pre-',num2str(c1),'-post-',num2str(c2)])
            grid
       
        else
            continue
        
        end
    end
end