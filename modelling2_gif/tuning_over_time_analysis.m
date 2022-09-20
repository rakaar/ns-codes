close all;
% every 5 batches
n_columns = 21;
n_steps = 50;
ba_from_middle = zeros(n_steps,n_columns);
starting_bbb = 251:5:500;
data_path = "D:\21cols-som-off-BA-training-from-middle";

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
         ba_from_middle(bbb, ccc) = mean(cols_spike_rates_subatch(:,ccc));
    end
    

end % end of bbb
save('ba_from_middle_tuning_over_time.mat', 'ba_from_middle')

%% 

for i=1:10:50
    figure
        hold on
            plot(initial_ab_col_spike_rates, 'LineStyle','--','color','b','LineWidth',3)
            plot(ab_col_spike_rates(i,:), 'b','LineWidth',4)
            plot(initial_ba_col_spike_rates,'LineStyle','--','color','r','LineWidth',3)
            plot(ba_col_spike_rates(i,:),'r','LineWidth',4)
            plot(ba_from_middle(i,:), 'LineWidth',4)
        hold off
        legend('inital ab','ab','initial ba','ba','BA from Middle')
    grid
end