close all;
% every 5 batches
n_columns = 21;
n_steps = 50;
ab_col_spike_rates = zeros(n_steps,n_columns);
starting_bbb = 251:5:500;
data_path = "D:\som-off-21cols-data";

for bbb=1:length(starting_bbb)
    fprintf("\n bbb is %d \n", bbb)

    cols_spike_rates_subatch = zeros(5,n_columns);

    b2_starting = starting_bbb(bbb):starting_bbb(bbb)+4;
    for b2=1:length(b2_starting)
        previous_batch_file = strcat(data_path,'\batch_', num2str(b2_starting(b2)), '.mat');
        spikes = load(previous_batch_file, 'spikes').spikes;

        for ccc=1:n_columns
            cols_spike_rates_subatch(b2, ccc) = mean( mean(squeeze(spikes(1,ccc,:,:)), 1) );
        end
    end % end of b2

    for ccc=1:n_columns
         ab_col_spike_rates(bbb, ccc) = mean(cols_spike_rates_subatch(:,ccc));
    end
    

end % end of bbb


%% generate imagesc
figure
    imagesc(ab_col_spike_rates.')
title('ab')
load('BA_tuning_over_time.mat')
figure
    imagesc(ba_col_spike_rates.')
title('ba')