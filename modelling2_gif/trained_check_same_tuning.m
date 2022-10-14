close all;
% every 5 batches
n_columns = 21;
n_steps = 1;
ba_from_last = zeros(n_steps,n_columns);
starting_bbb = 495;
% starting_bbb = 1:5;
data_path = "D:\som-off-21cols-10-12-BA-data";
fname = 'c21_10_12_off_BA_trained_BA_final.mat';

for bbb=1:length(starting_bbb)
    fprintf("\n bbb is %d \n", bbb)

    cols_spike_rates_subatch = zeros(5,n_columns);

    b2_starting = starting_bbb(bbb)+1:starting_bbb(bbb)+5;
    for b2=1:length(b2_starting)
        previous_batch_file = strcat(data_path,'\batch_', num2str(b2_starting(b2)), '.mat');
        spikes = load(previous_batch_file, 'spikes').spikes;

        for ccc=1:n_columns
            cols_spike_rates_subatch(b2, ccc) = mean( mean(squeeze(spikes(1,ccc,:,:)), 1) );
        end
    end % end of b2

    for ccc=1:n_columns
         ba_from_last(bbb, ccc) = mean(cols_spike_rates_subatch(:,ccc));
    end
    

end % end of bbb
save(fname, 'ba_from_last')