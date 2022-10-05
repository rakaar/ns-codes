close all;
% every 5 batches
n_columns = 21;
n_steps = 5;
ba_from_last = zeros(n_steps,n_columns);
starting_bbb = 251:50:500;
% starting_bbb = 1:5;
data_path = "D:\som-on-21cols-10-12-BA-data";

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
         ba_from_last(bbb, ccc) = mean(cols_spike_rates_subatch(:,ccc));
    end
    

end % end of bbb
save('ba_last.mat', 'ba_from_last')

%% 
close all

batch_set=1:10:50;
for i=1:length(batch_set)
    figure
        hold on
            plot(initial_ab_col_spike_rates, 'LineStyle','--','color','b','LineWidth',3)
            plot(ab_col_spike_rates(batch_set(i),:), 'b','LineWidth',4)
            plot(initial_ba_col_spike_rates,'LineStyle','--','color','r','LineWidth',3)
            plot(ba_col_spike_rates(batch_set(i),:),'r','LineWidth',4)
            plot(ba_from_middle(batch_set(i),:), 'LineWidth',4)
            plot(ba_given_ab_col_spike_rates(i,:),'LineWidth',5)
        hold off
        legend('inital ab','ab','initial ba','ba','BA from Middle', 'BA given AB')
    grid
    disp(i)
end

%% rough
close all
for i=1:5
figure
    hold on
    plot(initial_ab_col_spike_rates(5,:),'Color','b','LineStyle','--')
    plot(ba_given_ab_col_spike_rates(i,:), 'Color','r','LineStyle','--')

    plot(ba_from_initial(5,:),'Color','b','LineWidth',2)
    plot(ba_from_last(i,:),'Color','r','LineWidth',2)
    hold off
    legend('AB initial', 'AB later','BA initial','BA later')
    title('for BA tuned')
grid
end
%% test
close all
figure
    hold on
        plot(a_initial,'Color','b','LineStyle','--')
        plot(a_final, 'Color','r','LineStyle','--')

        plot(b_initial, 'Color','b','LineWidth',2)
        plot(b_final, 'Color','r','LineWidth',2)
    hold off
    legend('A initial', 'A final','B initial','B final')
    title('BA trained')
grid

