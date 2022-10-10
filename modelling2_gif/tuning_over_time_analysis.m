close all;
% every 5 batches
n_columns = 21;
n_steps = 1;
ba_from_last = zeros(n_steps,n_columns);
starting_bbb = 5;
% starting_bbb = 1:5;
data_path = "D:\som-on-21cols-10-12-BA-data";

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
save('AB_initial.mat', 'ba_from_last')

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
close all;
figure
    hold on
        plot(ab_i,'Color','b','LineStyle','--')
        plot(ab_f, 'Color','r','LineStyle','--')

        plot(ba_i, 'Color','b','LineWidth',2)
        plot(ba_f, 'Color','r','LineWidth',2)
    hold off
    legend('AB initial', 'AB final','BA initial','BA final')
    title('9 13 on BA trained')
grid

%% ab trained    matrix

AB_trained_matrix = [...
    c21_10_12_on_ABi;...
    c21_10_12_on_ABf;...
    c21_10_12_on_BAi;...
    c21_10_12_on_BAf;...

    c21_10_12_off_ABi;...
    c21_10_12_off_ABf;...
    c21_10_12_off_BAi;...
    c21_10_12_off_BAf;...

    c21_9_13_on_ABi;...
    c21_9_13_on_ABf;...
    c21_9_13_on_BAi;...
    c21_9_13_on_BAf;...

    c21_9_13_off_ABi;...
    c21_9_13_off_ABf;...
    c21_9_13_off_BAi;...
    c21_9_13_off_BAf;...
    ];
