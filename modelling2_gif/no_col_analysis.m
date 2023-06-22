%%  response to AB at initial and final stages
neuron_types = load("E:\RK_modelling\data\batch_1.mat").shuffled_neuron_types;
data_path = "E:\RK_modelling\data";
% for initial tuning
max_batches = 300;

initial_exc_spike_rate =  zeros(80,1);
for b=5:9
    data_file = strcat(data_path, '\', 'batch_', num2str(b), '.mat');
    all_spikes = load(data_file).spikes;
    exc_spikes = all_spikes(find(neuron_types == 0),:);
    exc_mean_spikes = mean(exc_spikes,2); % 80 x 1
    initial_exc_spike_rate = initial_exc_spike_rate + exc_mean_spikes;
end
initial_exc_spike_rate = initial_exc_spike_rate/5;
initial_exc_spike_rate = initial_exc_spike_rate/(0.001);

final_exc_spike_rate =  zeros(80,1);
for b=246:250
    data_file = strcat(data_path, '\', 'batch_', num2str(b), '.mat');
    all_spikes = load(data_file).spikes;
    exc_spikes = all_spikes(find(neuron_types == 0),:);
    exc_mean_spikes = mean(exc_spikes,2); % 80 x 1
    final_exc_spike_rate = final_exc_spike_rate + exc_mean_spikes;
end
final_exc_spike_rate = final_exc_spike_rate/5;
final_exc_spike_rate = final_exc_spike_rate/(0.001);

figure
    hold on
        plot(initial_exc_spike_rate)
        plot(final_exc_spike_rate)
    hold off
    legend('initial', 'final')
grid

figure
    plot(final_exc_spike_rate - initial_exc_spike_rate)
grid

%% weight distribution initial and final
initial_weight_matrix = load("E:\RK_modelling\data\batch_1.mat").weight_matrix;
final_weight_matrix = load("E:\RK_modelling\data\batch_250.mat").weight_matrix;
neuron_types = load("E:\RK_modelling\data\batch_1.mat").shuffled_neuron_types;

exc_initial_weights = initial_weight_matrix(find(neuron_types == 0), find(neuron_types == 0),:);
exc_initial_weights_mean = mean(exc_initial_weights, 3);
exc_initial_weights_mean_reshaped = reshape(exc_initial_weights_mean, 6400,1);

exc_final_weights = final_weight_matrix(find(neuron_types == 0), find(neuron_types == 0), :);
exc_final_weights_mean = mean(exc_final_weights, 3);
exc_final_weights_mean_reshaped = reshape(exc_final_weights_mean, 6400,1);

figure
    hist(exc_initial_weights_mean_reshaped)
grid

figure
    hist(exc_final_weights_mean_reshaped)
 grid

 %% eigenvalues ?
 initial_weight_matrix = load("E:\RK_modelling\data\batch_1.mat").weight_matrix;
final_weight_matrix = load("E:\RK_modelling\data\batch_250.mat").weight_matrix;
neuron_types = load("E:\RK_modelling\data\batch_1.mat").shuffled_neuron_types;

exc_initial_weights = initial_weight_matrix(find(neuron_types == 0), find(neuron_types == 0),:);
exc_initial_weights_mean = mean(exc_initial_weights, 3);
exc_initial_weights_mean_squeeze = squeeze(exc_initial_weights_mean);
eig_initial = eig(exc_initial_weights_mean_squeeze);
real_eig_initial = real(eig_initial);
img_eig_initial = imag(eig_initial);

exc_final_weights = final_weight_matrix(find(neuron_types == 0), find(neuron_types == 0), :);
exc_final_weights_mean = mean(exc_final_weights, 3);
exc_final_weights_mean_squeeze = squeeze(exc_final_weights_mean);
eig_final = eig(exc_final_weights_mean_squeeze);
real_eig_final = real(eig_final);
img_eig_final = imag(eig_final);

figure
    hold on
        plot(real_eig_initial)
        plot(img_eig_initial)

        plot(real_eig_final)
        plot(img_eig_final)
    hold off
    legend('initial-real', 'initial-img', 'final-real','final-img')
grid
%% spike rate over time
spike_rate_over_time = zeros(80,max_batches);
neuron_types = load("E:\RK_modelling\data\batch_1.mat").shuffled_neuron_types;
data_path = "E:\RK_modelling\data";

for b=1:max_batches
    disp(b)
    data_file = strcat(data_path, '\', 'batch_', num2str(b), '.mat');
    all_spikes = load(data_file).spikes;
    exc_spikes = all_spikes(find(neuron_types == 0),:);
    mean_exc_spikes = mean(exc_spikes, 2);
    spike_rate_over_time(:,b) = mean_exc_spikes;
end

figure
    plot(spike_rate_over_time.')
grid

%% set of neurons in 1 ensemble
last_batch_spike_rates = spike_rate_over_time(:,max_batches);

set1_indices = find(last_batch_spike_rates >= 0.075 & last_batch_spike_rates <= 0.085 );
set2_indices = find(last_batch_spike_rates >= 0.05 & last_batch_spike_rates <= 0.056);
final_weight_matrix = load("E:\RK_modelling\data\batch_" + num2str(max_batches) +  ".mat").weight_matrix;
final_weight_matrix = final_weight_matrix(find(neuron_types == 0), find(neuron_types == 0), :);
set1_weights_matrix = squeeze(final_weight_matrix(set1_indices, set1_indices, end));
set2_weights_matrix = squeeze(final_weight_matrix(set2_indices, set2_indices, end));

figure
    imagesc(set1_weights_matrix)
    title('set1')
grid

figure
    imagesc(set2_weights_matrix)
    title('set2')
grid
%% 
w1 = load("E:\RK_modelling\data\batch_1.mat").weight_matrix;
w500 = load("E:\RK_modelling\data\batch_" + num2str(max_batches) +  ".mat").weight_matrix;

neuron_types = load("E:\RK_modelling\data\batch_1.mat").shuffled_neuron_types;

w_begin = squeeze(w1(find(neuron_types == 0),find(neuron_types == 0),1));
w_end = squeeze(w500(find(neuron_types == 0),find(neuron_types == 0),end));

figure
    imagesc(w_begin)
grid

figure
    imagesc(w_end)
grid

%% 
path = "E:\RK_modelling\data\batch_";

initial_AB = zeros(100,5);
final_AB = zeros(100,5);

for b=6:10
    spikes = load(strcat(path, num2str(b), '.mat')).spikes;
    initial_AB(:,b-5) = mean(spikes,2);
end

for b=max_batches-4:max_batches
    spikes = load(strcat(path, num2str(b), '.mat')).spikes;
    final_AB(:,b-(max_batches - 5)) = mean(spikes,2);
end

%% 
% initial_A = load('initial_A.mat').initial_A; 
% final_A = load('final_A.mat').final_A; 

% initial_B = load('initial_B.mat').initial_B; 
% final_B = load('final_B.mat').final_B; 

% initial_BA = load('initial_BA.mat').initial_BA; 
% final_BA = load('final_BA.mat').final_BA; 
initial_A = load('initial_response_to_A.mat').response_to_A;
final_A = load('final_response_to_A.mat').response_to_A;

initial_B = load('initial_response_to_B.mat').response_to_B;
final_B = load('final_response_to_B.mat').response_to_B;

initial_BA = load('initial_response_to_BA.mat').response_to_BA;
final_BA = load('final_response_to_BA.mat').response_to_BA;

initial_A = mean(initial_A,2);
final_A = mean(final_A,2);

initial_B = mean(initial_B,2);
final_B = mean(final_B,2);

initial_BA = mean(initial_BA,2);
final_BA = mean(final_BA,2);

initial_AB = mean(initial_AB,2);
final_AB = mean(final_AB,2);

all_tunings = zeros(100,8);
all_tunings(:,1) = initial_A;
all_tunings(:,2) = final_A;

all_tunings(:,3) = initial_B;
all_tunings(:,4) = final_B;

all_tunings(:,5) = initial_AB;
all_tunings(:,6) = final_AB;

all_tunings(:,7) = initial_BA;
all_tunings(:,8) = final_BA;

initial_response_to_AB = initial_AB;
final_response_to_AB = final_AB;

save('initial_response_to_AB.mat', 'initial_response_to_AB')
save('final_response_to_AB.mat', 'final_response_to_AB')

figure
    imagesc(all_tunings)
    title(strcat('A-','B-','AB-','BA-', '-i,f x 4'))
grid

figure
    plot(all_tunings.')
    title(strcat('A-','B-','AB-','BA-', '-i,f x 4'))
grid

%% statistical analysis
% - A
disp(['On A stimulus on AB trained matrix: ttest and ranksun test btn initial and final'])
initial_res = load('initial_response_to_A.mat').response_to_A;
final_res = load('final_response_to_A.mat').response_to_A;
ttest1 = zeros(100,1);
ranksum1 = zeros(100,1);
for n=1:100
    ttest1(n) = ttest(initial_res(n,:), final_res(n,:));
    [~,ranksum1(n)] = ranksum(initial_res(n,:), final_res(n,:));
end

disp((['Ttest: ' num2str(sum(ttest1)) '/' num2str(length(ttest1)) ' have sig difference' ]))
disp((['Ranksum: ' num2str(sum(ranksum1)) '/' num2str(length(ranksum1)) ' have sig difference' ]))
disp('-----------------------------------------------------------')
% - B
disp(['On B stimulus on AB trained matrix: ttest and ranksun test btn initial and final'])
initial_res = load('initial_response_to_B.mat').response_to_B;
final_res = load('final_response_to_B.mat').response_to_B;
ttest1 = zeros(100,1);
ranksum1 = zeros(100,1);
for n=1:100
    ttest1(n) = ttest(initial_res(n,:), final_res(n,:));
    [~,ranksum1(n)] = ranksum(initial_res(n,:), final_res(n,:));
end

disp((['Ttest: ' num2str(sum(ttest1)) '/' num2str(length(ttest1)) ' have sig difference' ]))
disp((['Ranksum: ' num2str(sum(ranksum1)) '/' num2str(length(ranksum1)) ' have sig difference' ]))
disp('-----------------------------------------------------------')
% - AB
path = "E:\RK_modelling\data\batch_";

initial_AB = zeros(100,5);
final_AB = zeros(100,5);

for b=6:10
    spikes = load(strcat(path, num2str(b), '.mat')).spikes;
    initial_AB(:,b-5) = mean(spikes,2);
end

for b=max_batches-4:max_batches
    spikes = load(strcat(path, num2str(b), '.mat')).spikes;
    final_AB(:,b-(max_batches - 5)) = mean(spikes,2);
end

disp(['On AB stimulus on AB trained matrix: ttest and ranksun test btn initial and final'])
initial_res = initial_AB;
final_res = final_AB;
ttest1 = zeros(100,1);
ranksum1 = zeros(100,1);
for n=1:100
    ttest1(n) = ttest(initial_res(n,:), final_res(n,:));
    [~,ranksum1(n)] = ranksum(initial_res(n,:), final_res(n,:));
end
disp((['Ttest: ' num2str(nansum(ttest1)) '/' num2str(length(ttest1)) ' have sig difference' ]))
disp((['Ranksum: ' num2str(sum(ranksum1)) '/' num2str(length(ranksum1)) ' have sig difference' ]))
disp('-----------------------------------------------------------')

% - BA
disp(['On BA stimulus on AB trained matrix: ttest and ranksun test btn initial and final'])
initial_res = load('initial_response_to_BA.mat').response_to_BA;
final_res = load('final_response_to_BA.mat').response_to_BA;
ttest1 = zeros(100,1);
ranksum1 = zeros(100,1);
for n=1:100
    ttest1(n) = ttest(initial_res(n,:), final_res(n,:));
    [~,ranksum1(n)] = ranksum(initial_res(n,:), final_res(n,:));
end

disp((['Ttest: ' num2str(sum(ttest1)) '/' num2str(length(ttest1)) ' have sig difference' ]))
disp((['Ranksum: ' num2str(sum(ranksum1)) '/' num2str(length(ranksum1)) ' have sig difference' ]))
disp('-----------------------------------------------------------')


disp(['-----------Between AB and BA : right------------'])
ttest1 = zeros(100,1);
ranksum1 = zeros(100,1);
for n=1:100
    ttest1(n) = ttest2(final_AB(n,:), final_BA(n,:), 'tail', 'right');
    [~,ranksum1(n)] = ranksum(final_AB(n,:), final_BA(n,:), 'tail', 'right');
end
disp(['Ttest: ' num2str(nansum(ttest1)) '/' num2str(length(ttest1)) ' have sig difference' ])
disp(['Ranksum: ' num2str(nansum(ranksum1)) '/' num2str(length(ranksum1)) ' have sig difference' ])
%% csi analysis

% a_csi = (final_A - initial_A)./(final_A + initial_A);
% b_csi = (final_B - initial_B)./(final_B + initial_B);
% ab_csi = (final_AB - initial_AB)./(final_AB + initial_AB);
% ba_csi = (final_BA - initial_BA)./(final_BA + initial_BA);

% figure
%     scatter(a_csi, b_csi, '*')

% grid

% figure
%     scatter(ab_csi, ba_csi, '*')

% grid

% a_and_b =[a_csi b_csi];

% id = kmeans(a_and_b, 2);


% ab_and_ba =[ab_csi ba_csi];

% id2 = kmeans(ab_and_ba, 2);

% figure
%     hold on
%         scatter(a_csi(find(id == 1)), b_csi(find(id == 1)), '*r')
%         scatter(a_csi(find(id == 2)), b_csi(find(id == 2)), '*b')
%     hold off
% grid

% figure
%     hold on
%         scatter(ab_csi(find(id2 == 1)), ba_csi(find(id2 == 1)), '*r')
%         scatter(ab_csi(find(id2 == 2)), ba_csi(find(id2 == 2)), '*b')
%     hold off
% grid

% %% csi diff
% initial_b_wrt_a = (initial_B - initial_A)./(initial_B + initial_A);
% final_b_wrt_a = (final_B - final_A)./(final_B + final_A);

% batch_1_path = "E:\RK_modelling\data\batch_1.mat";
% shuffled_neuron_types = load(batch_1_path, 'shuffled_neuron_types').shuffled_neuron_types;

% figure
%     hold on
%         plot( initial_b_wrt_a(  find(shuffled_neuron_types == 0)  )  )
%         plot( final_b_wrt_a(    find(shuffled_neuron_types == 0)  )  )
%     hold off
%     legend('initial b', 'final b')
% grid