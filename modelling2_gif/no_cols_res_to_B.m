clear all
tic
response_to_B = zeros(100,5);
shuffled_neuron_types = load("D:\no_col_2000tokens_more_ee_var\batch_1.mat", 'shuffled_neuron_types').shuffled_neuron_types;
weight_matrix = load("D:\no_col_2000tokens_more_ee_var\batch_496.mat").weight_matrix;
for rep=1:5
n_exc = 80;
n_pv = 12; 
n_som = 8;
n_total_neurons = n_exc + n_pv + n_som;

n_thalamus = 100;

% time step
n_tokens = 4; 
single_stimulus_duration = 50;
pre_token_silence = 10;
post_token_silence = 10;
token_start_times = zeros(n_tokens,1);

% plastic weights
Amp_strength = 0.015; Amp_weak = 0.021;
tau_strength = 30; tau_weak = 50;

% total time
physical_time_in_ms = 1; %dt time step
dt = 1;  % 0.2 dt = 20 ms, so 0 .01 = 1 ms
t_simulate = n_tokens*(single_stimulus_duration + pre_token_silence + post_token_silence);
tspan = 1:dt:t_simulate;




% neuron_types = [2*ones(n_som,1); zeros(n_exc,1); 1*ones(n_pv,1)];
% shuffled_indices = randperm(n_total_neurons);
% shuffled_neuron_types = neuron_types(shuffled_indices);
% IMPORT




lamda_1 = 300;
lamda_2 = 100;
lamda_3 = 50;
lamda_4 = 25;

% thalamus tuning
n_thalamic_cols = 8;
n_thalamic_neurons = 20;
thalamic_cols_spike_rates =  zeros(n_thalamic_cols,n_thalamic_neurons,length(tspan));
for tok=1:n_tokens
    ind =  (tok-1)*(single_stimulus_duration + pre_token_silence + post_token_silence) + 1;
    
    %%% first stim
    token_first_half_start_time = ind+pre_token_silence;
    token_first_half_end_time = token_first_half_start_time+single_stimulus_duration-1;
    thalamic_cols_spike_rates(1,:,token_first_half_start_time:token_first_half_end_time) = lamda_4;
    thalamic_cols_spike_rates(2,:,token_first_half_start_time:token_first_half_end_time) = lamda_4;
    thalamic_cols_spike_rates(3,:,token_first_half_start_time:token_first_half_end_time) = lamda_3;
    thalamic_cols_spike_rates(4,:,token_first_half_start_time:token_first_half_end_time) = lamda_2;
    thalamic_cols_spike_rates(5,:,token_first_half_start_time:token_first_half_end_time) = lamda_1;
    thalamic_cols_spike_rates(6,:,token_first_half_start_time:token_first_half_end_time) = lamda_2;
    thalamic_cols_spike_rates(7,:,token_first_half_start_time:token_first_half_end_time) = lamda_3;
    thalamic_cols_spike_rates(8,:,token_first_half_start_time:token_first_half_end_time) = lamda_4;


    %%% silence
    token_gap_duration_start = ind+pre_token_silence+single_stimulus_duration;
    token_gap_duration_end = token_gap_duration_start+post_token_silence-1;
    thalamic_cols_spike_rates(:,:,token_gap_duration_start:token_gap_duration_end) = 0;

   
    %%% register token start time
    token_start_times(tok) = ind;
end

% epsc produced by thalamic
xr_thalamic = ones(n_thalamic_cols, n_thalamic_neurons,length(tspan));
xe_thalamic = zeros(n_thalamic_cols, n_thalamic_neurons,length(tspan));
xi_thalamic = zeros(n_thalamic_cols, n_thalamic_neurons,length(tspan));

tau_re_thalamic = 0.3; tau_ir_thalamic = 300; tau_ei_thalamic = 50;
tau_re = 0.6; tau_ir = 700; tau_ei = 15;

thalamic_epsc = zeros(n_thalamic_cols, n_thalamic_neurons, length(tspan));

for c=1:n_thalamic_cols
    for n=1:n_thalamic_neurons
        spike_rate = squeeze(thalamic_cols_spike_rates(c,n,:));
        poisson_spikes = poisson_generator(spike_rate, 1);
        g_t = 0.25*get_g_t_vector(poisson_spikes, length(tspan));
        
        for t=2:length(tspan)
            current_xr_thalamic = xr_thalamic(c,n,t-1);
            current_xe_thalamic = xe_thalamic(c,n,t-1);
            current_xi_thalamic = xi_thalamic(c,n,t-1);
        
            M = 0;
            if poisson_spikes(t) == 1
                M = 1;
            end

            new_xr_thalamic = update_xr(M, current_xr_thalamic, current_xi_thalamic, tau_re_thalamic, tau_ir_thalamic);
            new_xe_thalamic = update_xe(M, current_xr_thalamic, current_xe_thalamic, tau_re_thalamic, tau_ei_thalamic);
            new_xi_thalamic = update_xi(current_xe_thalamic, current_xi_thalamic, tau_ei_thalamic, tau_ir_thalamic);
        
            new_xr_thalamic = limit_x(new_xr_thalamic);
            new_xe_thalamic = limit_x(new_xe_thalamic);
            new_xi_thalamic = limit_x(new_xi_thalamic);

            xr_thalamic(c,n,t) = new_xr_thalamic;
            xe_thalamic(c,n,t) = new_xe_thalamic;
            xi_thalamic(c,n,t) = new_xi_thalamic;

            if ismember(t,token_start_times)
                    xr_thalamic(:,:,t) = 1;
                    xe_thalamic(:,:,t) = 0;
                    xi_thalamic(:,:,t) = 0;
            end

        end % end of t

        epsc = reshape(g_t, length(tspan),1).*squeeze(xe_thalamic(c,n,:));
        thalamic_epsc(c,n,:) = epsc;
        

    end % end of n
end % end of c

%%% thalamic weights
thalamus_to_exc_epsc = zeros(n_exc, length(tspan));
thalamus_to_pv_epsc = zeros(n_pv, length(tspan));
thalamus_to_som_epsc = zeros(n_som, length(tspan));

% weights will be same, epsc will change
% rand_weights_thalamus_to_exc = zeros(n_exc,5); 
% max_weight_thalamus_to_l4 = 250;
% rand_weights_thalamus_to_exc(:,1) = make_rand_vector(max_weight_thalamus_to_l4/4,10,[n_exc 1]);
% rand_weights_thalamus_to_exc(:,2) = make_rand_vector(max_weight_thalamus_to_l4/2,10,[n_exc 1]);
% rand_weights_thalamus_to_exc(:,3) = make_rand_vector(max_weight_thalamus_to_l4,10,[n_exc 1]);
% rand_weights_thalamus_to_exc(:,4) = make_rand_vector(max_weight_thalamus_to_l4/2,10,[n_exc 1]);
% rand_weights_thalamus_to_exc(:,5) = make_rand_vector(max_weight_thalamus_to_l4/4,10,[n_exc 1]);
% JUST IMPORT IT FROM PREVIOUS BATCH/FIRST BATCH - weights from thalamus to
% exc are constant over time, but vary over L4 neurons.
rand_weights_thalamus_to_exc = load("D:\no_col_2000tokens_more_ee_var\batch_1.mat",'rand_weights_thalamus_to_exc').rand_weights_thalamus_to_exc;

weight_thalamus_to_pv_weak = 77.5; weight_thalamus_to_pv_moderate = 155; weight_thalamus_to_pv_strong = 310;
weights_thalamus_to_pv = [weight_thalamus_to_pv_weak,weight_thalamus_to_pv_moderate, weight_thalamus_to_pv_strong, weight_thalamus_to_pv_moderate, weight_thalamus_to_pv_weak];

weight_thalamus_to_som_weak = 65; weight_thalamus_to_som_moderate = 130; weight_thalamus_to_som_strong = 260;
weights_thalamus_to_som = [weight_thalamus_to_som_weak,weight_thalamus_to_som_moderate, weight_thalamus_to_som_strong, weight_thalamus_to_som_moderate, weight_thalamus_to_som_weak];

tuning_types = zeros(4,5); % 4 cols, for each col 5 thal cols % all cols imaginary due to shuffle
tuning_types(1,:) = [1 2 3 4 5];
tuning_types(2,:) = [2 3 4 5 6];
tuning_types(3,:) = [3 4 5 6 7];
tuning_types(4,:) = [4 5 6 7 8];


for i=1:4
    % thalamus to exc
    for j=1:20
        for k=1:5
            thalamus_to_exc_epsc((i-1)*20 + j, :) = reshape(thalamus_to_exc_epsc((i-1)*20 + j, :),  1,length(tspan)) + reshape(rand_weights_thalamus_to_exc((i-1)*20 + j,k).*squeeze(thalamic_epsc(tuning_types(i,k),randi([1 20]),:)),  1,length(tspan));
        end
    end

    % thalamus to pv 
    for j=1:3
        for k=1:5
            thalamus_to_pv_epsc((i-1)*3 + j,:) = reshape(thalamus_to_pv_epsc((i-1)*3 + j, :),  1,length(tspan)) + reshape(  weights_thalamus_to_pv(k)*squeeze( thalamic_epsc(tuning_types(i,k), randi([1 20]),: ) ), 1,length(tspan)  );
        end
    end

    % thalamus to som
    for j=1:2
        for k=1:5
            thalamus_to_som_epsc((i-1)*2 + j,:) = reshape(thalamus_to_som_epsc((i-1)*2 + j, :),  1,length(tspan)) + reshape(  weights_thalamus_to_som(k)*squeeze( thalamic_epsc(tuning_types(i,k), randi([1 20]),: ) ), 1,length(tspan)  );
        end
    end

end % end of i

% thalamus_to_exc_shuffle_indices = randperm(n_exc); % IMPORT
thalamus_to_exc_shuffle_indices = load("D:\no_col_2000tokens_more_ee_var\batch_1.mat", 'thalamus_to_exc_shuffle_indices').thalamus_to_exc_shuffle_indices;
thalamus_to_exc_epsc = thalamus_to_exc_epsc(thalamus_to_exc_shuffle_indices,:);

% thalamus_to_pv_shuffle_indices = randperm(n_pv); % IMPORT
thalamus_to_pv_shuffle_indices = load("D:\no_col_2000tokens_more_ee_var\batch_1.mat",'thalamus_to_pv_shuffle_indices').thalamus_to_pv_shuffle_indices;
thalamus_to_pv_epsc = thalamus_to_pv_epsc(thalamus_to_pv_shuffle_indices,:);

% thalamus_to_som_shuffle_indices = randperm(n_som); % IMPORT
thalamus_to_som_shuffle_indices = load("D:\no_col_2000tokens_more_ee_var\batch_1.mat", 'thalamus_to_som_shuffle_indices').thalamus_to_som_shuffle_indices;
thalamus_to_som_epsc = thalamus_to_som_epsc(thalamus_to_som_shuffle_indices,:);

thalamus_to_all_l4_epsc = zeros(n_total_neurons, length(tspan));
exc_counter = 1;
pv_counter = 1;
som_counter = 1;
for n=1:n_total_neurons
    if shuffled_neuron_types(n) == 0 % exc
        thalamus_to_all_l4_epsc(n,:) = thalamus_to_exc_epsc(exc_counter,:);
        exc_counter = exc_counter + 1;
    elseif shuffled_neuron_types(n) == 1 % pv 
        thalamus_to_all_l4_epsc(n,:) = thalamus_to_pv_epsc(pv_counter,:);
        pv_counter = pv_counter + 1;
    else % som
        thalamus_to_all_l4_epsc(n,:) = thalamus_to_som_epsc(som_counter,:);
        som_counter = som_counter + 1;
    end
end

% kernel for g(t)
tau_syn = 10;
kernel_kt = [0 exp(-[0:t_simulate])./tau_syn];

tau_syn_longer = 45;
kernel_kt_longer = [0 exp(-[0:t_simulate])./tau_syn_longer];

rec_exc_epsc = zeros(n_total_neurons, length(tspan));
rec_pv_epsc = zeros(n_total_neurons, length(tspan));
rec_som_epsc = zeros(n_total_neurons, length(tspan));

spikes = zeros(n_total_neurons, length(tspan));

voltages = zeros(n_total_neurons,length(tspan));
i1s = zeros(n_total_neurons,length(tspan));
i2s = zeros(n_total_neurons,length(tspan));
thetas = zeros(n_total_neurons,length(tspan));

xr = zeros(n_total_neurons, length(tspan));
xe = zeros(n_total_neurons, length(tspan));
xi = zeros(n_total_neurons, length(tspan));

% initialize
v0 = -70;  
xr(:, 1:5) = 1; 
voltages(:, 1:5) = v0; % 
i1s(:, 1:5) = 0.01;
i2s(:, 1:5) = 0.001;
thetas(:, 1:5) = -50.0;
% inhibitory synapses are not depressing
xe(shuffled_neuron_types ~= 0, :) = 1;

% sponataneous current into l4 neurons
background_epsc = zeros(n_total_neurons, length(tspan));
for n=1:n_total_neurons
        background_lamda = 0.75*ones(1, length(tspan));
        background_spikes = poisson_generator(background_lamda, 1);
        background_epsc1 = 5*get_g_t_vector_background(background_spikes, length(tspan));
        background_epsc_reshaped = reshape(background_epsc1, 1,length(tspan));
        background_epsc(n,:) = background_epsc_reshaped;
end

% test
g_t_test = zeros(n_total_neurons, n_total_neurons, length(tspan));

for t=6:length(tspan)
    for n=1:n_total_neurons
       thalamic_epsc_for_n = thalamus_to_all_l4_epsc(n,t);
       
       rec_epsc_for_n = 0;
       
       for pre=1:n_total_neurons
            if pre == n
                continue
            end

            pre_spike_train = spikes(pre,:);
            if shuffled_neuron_types(pre) == 0
                pre_g_t = get_g_t(pre_spike_train, dt, t-5, tspan,kernel_kt);
            else
                pre_g_t = get_g_t(pre_spike_train, dt, t-5, tspan,kernel_kt_longer);
            end
            g_t_test(n,pre,t) = g_t_test(n,pre,t) + pre_g_t; 
            pre_xe = xe(pre, t-5);
            pre_weight = weight_matrix(pre,n,t-5);
            pre_epsc = pre_g_t*pre_xe*pre_weight;
          
            rec_epsc_for_n = rec_epsc_for_n + pre_epsc;
            

            if shuffled_neuron_types(pre) == 0
                rec_exc_epsc(n,t) = rec_exc_epsc(n,t) + pre_epsc;
            elseif shuffled_neuron_types(pre) == 1
                rec_pv_epsc(n,t) = rec_pv_epsc(n,t) + pre_epsc;
            else 
                rec_som_epsc(n,t) = rec_som_epsc(n,t) + pre_epsc;
            end

       end % end of pre

        
    total_epsc_for_n = thalamic_epsc_for_n + rec_epsc_for_n;

    spike_vec = spikes(n,:);
    latest_spike_time = -1;
    for s=t-1:-1:t-5
        if s >= 1 && spike_vec(s) == 1
            latest_spike_time = s;
            break;
        end
    end

    [voltages(n,t), i1s(n,t), i2s(n,t), thetas(n,t), spikes(n,t)] = calculate_new_state_dynamic_threshold_rule(voltages(n,t-1), i1s(n,t-1), i2s(n,t-1), thetas(n,t-1), total_epsc_for_n, background_epsc(n,t),dt,t,latest_spike_time);
          
    if spikes(c,t-1) == 1
		spikes(c,t) = 0;
    end

     M = 0;
     if spikes(n,t) == 1
         M = 1;
     end

    current_xr = xr(n,t-1);
    current_xe = xe(n,t-1);
    current_xi = xi(n,t-1);

    if shuffled_neuron_types(n) == 0
        new_xr = update_xr(M, current_xr, current_xi, tau_re, tau_ir);
        new_xe = update_xe(M, current_xr, current_xe, tau_re, tau_ei);

        new_xr = limit_x(new_xr);
        new_xe = limit_x(new_xe);
        new_xi = 1 - (new_xr + new_xe);

        xr(n,t) = new_xr;
        xe(n,t) = new_xe;
        xi(n,t) = new_xi;
    end


                
    end % end of n

    

    % if t is initializing of a token, reset values

        if ismember(t,token_start_times)
            voltages(:,t) = v0;  
            xr(:, t) = 1;
            xe(:, t) = 0;
            xi(:, t) = 0;
            i1s(:, t) = 0.01;
            i2s(:, t) = 0.001;
            thetas(:,t) = -50.0;
        end

end % end of t

response_to_B(:,rep) = mean(spikes,2);
disp(rep)
end % end of rep
toc