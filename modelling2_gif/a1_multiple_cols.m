clear all;
close all;

tic
n_iters = 1;

% basic variables;
n_columns = 5;
n_excitatory = 20;
n_pv = 3; n_som = 2;
n_inhibitory = n_pv + n_som;
n_total_neurons = n_excitatory + n_inhibitory;
n_thalamic_neurons = 25;
n_thalamic_cols = 9;
    
% time step
n_tokens = 5;
pre_stimulus_time = 0; post_stimulus_time = 0;
single_stimulus_duration = 50; gap_duration = 50;
pre_token_silence = 10;
post_token_silence = 10;

token_start_times = zeros(n_tokens,1);

physical_time_in_ms = 1; %dt time step
dt = 1;  % 0.2 dt = 20 ms, so 0 .01 = 1 ms
t_simulate = n_tokens*(2*single_stimulus_duration + gap_duration + pre_token_silence + post_token_silence);
tspan = 0:dt:t_simulate;

% making bins
spike_rate_dt = 1*dt;
spike_rate_length = (length(tspan)-1)/(spike_rate_dt/dt);


% connection strength
% within column
inc_inh_to_exc_factor = 2.5;

J_ee_0 = 135;
J_pv_e_0 = 1.8750;
J_som_e_0 = 1.8750*3;

J_e_pv = -100*inc_inh_to_exc_factor;
J_pv_pv = -9.3750*inc_inh_to_exc_factor;
J_som_pv = 0;

J_e_som = -50*inc_inh_to_exc_factor;
J_som_som = 0;
J_pv_som = -9.3750*inc_inh_to_exc_factor;

% other column
J_ee_1 = 70;
J_pv_e_1 = 0.0131;
J_som_e_1 = 0.0131;

J_e_som_1 = -10;

J_ee_2 = 35;
J_pv_e_2 = 0.0056;
J_som_e_2 = 0.0056;

J_e_som_2 = -2;

% synaptic weight matrix - exc to exc - row: presyn, col: postsyn
exc_to_exc_weight_matrix = zeros(n_iters, n_columns, length(tspan),n_excitatory, n_excitatory);
minimum_weight_exc_to_exc = 5;
maximum_weight_exc_to_exc = 500;
across_columns_exc_to_exc_weight_matrix = zeros(n_iters, n_columns, n_columns, length(tspan), n_excitatory, n_excitatory);


% analysis of weights
num_of_LTPs = zeros(n_iters, n_columns, length(tspan));
num_of_LTDs = zeros(n_iters, n_columns, length(tspan));

% plasticity parameters
Amp_strength = 0.015; Amp_weak = 0.021;
tau_strength = 30; tau_weak = 50;

% kernel for g(t)
tau_syn = 10;
kernel_kt = [0 exp(-[0:t_simulate])./tau_syn];
    

% voltages and terms from it are 3d tensors
voltages = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
i1_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
i2_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
theta_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
spikes = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
spike_rates = zeros(n_iters, n_columns, n_total_neurons, spike_rate_length);

lamda = zeros(n_iters, n_thalamic_cols, n_thalamic_neurons, length(tspan));
thalamic_spikes = zeros(n_iters, n_thalamic_cols, n_thalamic_neurons, length(tspan));
epsc_thalamic = zeros(n_iters, n_thalamic_cols, n_thalamic_neurons, length(tspan));

I_background_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan)-1);
thalamic_epsc_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan)-1);
recurrence_exc_self_column_epsc_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan)-1);
recurrence_inh_self_column_epsc_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan)-1);
recurrence_exc_neighbour_column_epsc_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan)-1);
recurrence_inh_neighbour_column_epsc_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan)-1);

% synaptic resources
% --- l4 ---
xr = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
xe = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
xi = zeros(n_iters, n_columns, n_total_neurons, length(tspan));

% ---- thalamic ----
xr_thalamic = zeros(n_iters, n_thalamic_cols, n_thalamic_neurons, length(tspan));
xe_thalamic = zeros(n_iters, n_thalamic_cols, n_thalamic_neurons, length(tspan));
xi_thalamic = zeros(n_iters, n_thalamic_cols, n_thalamic_neurons, length(tspan));

% mapping from input thalamic neurons to a1 column neurons
lamda_i = 5;
lamda_a = 300;
lamda_b = 50;
lamda_m = 100;
lamda_s = 0;

lamda(:,:,:,1:pre_stimulus_time) = lamda_s;

for iter=1:n_iters
    for tok=1:n_tokens
        ind =  (tok-1)*(2*single_stimulus_duration + 1*gap_duration + pre_token_silence + post_token_silence) + 1;
        

        % ---- first half of token
        
        % pretoken silence
        lamda(iter,:,:,ind:ind+pre_token_silence-1) = lamda_s;

        % stimulus
        token_first_half_start_time = ind+pre_token_silence;
        token_first_half_end_time = token_first_half_start_time+single_stimulus_duration-1;

        lamda(iter,1,:,token_first_half_start_time:token_first_half_end_time) = lamda_i;
        lamda(iter,2,:,token_first_half_start_time:token_first_half_end_time) = lamda_b;
        lamda(iter,3,:,token_first_half_start_time:token_first_half_end_time) = lamda_m;
        lamda(iter,4,:,token_first_half_start_time:token_first_half_end_time) = lamda_a;
        lamda(iter,5,:,token_first_half_start_time:token_first_half_end_time) = lamda_m;
        lamda(iter,6,:,token_first_half_start_time:token_first_half_end_time) = lamda_b;
        lamda(iter,7:9,:,token_first_half_start_time:ind+pre_token_silence+single_stimulus_duration-1) = lamda_i;

        % posttoken silence
        token_gap_duration_start = ind+pre_token_silence+single_stimulus_duration;
        token_gap_duration_end = token_gap_duration_start+gap_duration-1;
        lamda(iter,:,:,token_gap_duration_start:token_gap_duration_end) = 0;

        % ---- second half of token
        % stimulus
        token_second_half_start_time = token_gap_duration_end;
        token_second_half_end_time = token_second_half_start_time + single_stimulus_duration - 1;
        lamda(iter,1,:,token_second_half_start_time:token_second_half_end_time) = lamda_i;
        lamda(iter,2,:,token_second_half_start_time:token_second_half_end_time) = lamda_i;
        lamda(iter,3,:,token_second_half_start_time:token_second_half_end_time) = lamda_i;
        lamda(iter,4,:,token_second_half_start_time:token_second_half_end_time) = lamda_b;
        lamda(iter,5,:,token_second_half_start_time:token_second_half_end_time) = lamda_m;
        lamda(iter,6,:,token_second_half_start_time:token_second_half_end_time) = lamda_a;
        lamda(iter,7,:,token_second_half_start_time:token_second_half_end_time) = lamda_m;
        lamda(iter,8,:,token_second_half_start_time:token_second_half_end_time) = lamda_b;
        lamda(iter,9,:,token_second_half_start_time:token_second_half_end_time) = lamda_i;
        % silence - saving time computationally
        
        
        lamda(iter,:,:,token_second_half_end_time:token_second_half_end_time+post_token_silence-1) = 0;
        
        token_start_times(tok,1) = token_second_half_end_time+post_token_silence;
    end
end

% thalamic connections based on neuron number
n_input_thalamic = 5;
thalamic_connections = zeros(n_total_neurons,n_input_thalamic);
% since there are 2 neurons in each thalamic column
% using coin flip outcomes/binary representation to create combinations
% coin_flip_outcomes = dec2bin(1:n_total_neurons);
% for n=1:n_total_neurons
%     outcome_str = coin_flip_outcomes(n,:);
%     outcome_split = split(outcome_str,'');
%     for z=2:6
%         str_connection_num = outcome_split{z};
%         connection_num = str2num(str_connection_num);
%         thalamic_connections(n,z-1) = connection_num + 1;
%     end
% 
% end

for n=1:n_total_neurons
    thalamic_connections(n,:) = randperm(n_thalamic_neurons,n_input_thalamic);
end

% weight_thalamic_to_exc_l4 = 550;
% weight_thalamic_to_pv_l4 = 750;
% weight_thalamic_to_som_l4 = 750;

weight_thalamic_to_exc_l4_above_col = 880;
weight_thalamic_to_exc_l4_side_col_1 = 440;
weight_thalamic_to_exc_l4_side_col_2 = 220;
weight_thalamic_to_exc_l4_arr = [weight_thalamic_to_exc_l4_side_col_2, weight_thalamic_to_exc_l4_side_col_1,weight_thalamic_to_exc_l4_above_col,weight_thalamic_to_exc_l4_side_col_1,weight_thalamic_to_exc_l4_side_col_2];

weight_thalamic_to_pv_l4_above_col = 1200;
weight_thalamic_to_pv_l4_side_col_1 = 600;
weight_thalamic_to_pv_l4_side_col_2 = 300;
weight_thalamic_to_pv_l4_arr = [weight_thalamic_to_pv_l4_side_col_2,weight_thalamic_to_pv_l4_side_col_1,weight_thalamic_to_pv_l4_above_col,weight_thalamic_to_pv_l4_side_col_1,weight_thalamic_to_pv_l4_side_col_2];

weight_thalamic_to_som_l4_above_col = 1200;
weight_thalamic_to_som_l4_side_col_1 = 600;
weight_thalamic_to_som_l4_side_col_2 = 300;
weight_thalamic_to_som_l4_arr = [weight_thalamic_to_som_l4_side_col_2,weight_thalamic_to_som_l4_side_col_1, weight_thalamic_to_som_l4_above_col, weight_thalamic_to_som_l4_side_col_1, weight_thalamic_to_som_l4_side_col_2];
%% time constant for synaptic resources
tau_re = 0.6; tau_ir = 700; tau_ei = 15;
tau_re_thalamic = 0.6; tau_ir_thalamic = 300; tau_ei_thalamic = 35;
    
% initialize
v0 = -70;  
xr(:, :, :, 1:5) = 1; 
xr_thalamic(:,:,:,1) = 1;
voltages(:, :, :, 1:5) = v0; % 
i1_tensor(:, :, :, 1:5) = 0.01;
i2_tensor(:, :, :, 1:5) = 0.001;
theta_tensor(:, :, :, 1:5) = -50.0;

J_ee_0_initial = 100;
exc_to_exc_weight_matrix(:, :, 1:5, :,:) = J_ee_0_initial;
connected_columns_arr = [  [1,2];[1 3];[2 3];[2 4];[3 4];[3 5];[4 5]; [2 1];[3 1];[3 2];[4 2];[4 3];[5 3];[5 4] ];
num_connected_pairs = size(connected_columns_arr,1);
for cc=1:num_connected_pairs
    c1 = connected_columns_arr(cc,1);
    c2 = connected_columns_arr(cc,2);
    if  abs(c2 - c1) == 1
        across_columns_exc_to_exc_weight_matrix(:,c1,c2,:,:,:) = J_ee_1;
    elseif abs(c2 - c1) == 2
        across_columns_exc_to_exc_weight_matrix(:,c1,c2,:,:,:) = J_ee_2;
    end
end

% sponataneous current into l4 neurons
background_epsc = zeros(n_iters,n_columns,n_total_neurons, length(tspan));
for iter=1:n_iters
    for col=1:n_columns
        for n=1:n_total_neurons
                background_lamda = 0.75*ones(1, length(tspan));
                background_spikes = poisson_generator(background_lamda, 1);
                background_epsc1 = 500*get_g_t_vector_background(background_spikes, length(tspan));
                background_epsc_reshaped = reshape(background_epsc1, 1,1,1,length(tspan));
                background_epsc(iter,col,n,:) = background_epsc_reshaped;
        end
    end
end


% inhibitory synapses are non-depressing
xe(:,:,n_excitatory+1:n_total_neurons,:) = 1;

for iter=1:n_iters
    
%     fprintf("------iter numm %d -----", iter);

    % thalamic
    for thal_col=1:n_thalamic_cols
        for thal_n=1:n_thalamic_neurons
            lamda_thal = squeeze(lamda(iter,thal_col,thal_n,:));
            lamda_thal = transpose(lamda_thal);
           thalamic_spikes(iter,thal_col,thal_n,:) = reshape(poisson_generator(lamda_thal,dt),  1,1,length(tspan));
        end
    end
    
    % depressing synapse thalamus -> l4
    for thal_col=1:n_thalamic_cols
        for thal_n=1:n_thalamic_neurons
            for t=2:length(tspan)
                
                

                current_xr_thalamic = xr_thalamic(iter,thal_col,thal_n, t-1);
                current_xe_thalamic = xe_thalamic(iter,thal_col,thal_n, t-1);
                current_xi_thalamic = xi_thalamic(iter,thal_col,thal_n, t-1);
                
                M = 0;
                if thalamic_spikes(iter,thal_col,thal_n,t) == 1
                    M = 1;
                end

                new_xr_thalamic = update_xr(M, current_xr_thalamic, current_xi_thalamic, tau_re_thalamic, tau_ir_thalamic);
                new_xe_thalamic = update_xe(M, current_xr_thalamic, current_xe_thalamic, tau_re_thalamic, tau_ei_thalamic);
                new_xi_thalamic = update_xi(current_xe_thalamic, current_xi_thalamic, tau_ei_thalamic, tau_ir_thalamic);

                if new_xr_thalamic > 1
                    new_xr_thalamic = 1;
                end
                if new_xr_thalamic < 0
                    new_xr_thalamic = 0;
                end

                if new_xe_thalamic > 1
                    new_xe_thalamic = 1;
                end
                if new_xe_thalamic < 0
                    new_xe_thalamic = 0;
                end

                if new_xi_thalamic > 1
                    new_xi_thalamic = 1;
                end
                if new_xi_thalamic < 0 
                    new_xi_thalamic = 0;
                end

                xr_thalamic(iter,thal_col,thal_n,t) = new_xr_thalamic;
                xe_thalamic(iter,thal_col,thal_n,t) = new_xe_thalamic;
                xi_thalamic(iter,thal_col,thal_n,t) = new_xi_thalamic;


                % initialize thalamic params 
                if ismember(t,token_start_times)
                    xr_thalamic(:,:,:,t) = 1;
                    xe_thalamic(:,:,:,t) = 0;
                    xi_thalamic(:,:,:,t) = 0;
                end
            
            end
        end
    end
    
    for thal_col=1:n_thalamic_cols
        for thal_n=1:n_thalamic_neurons
            epsc_thalamic(iter,thal_col,thal_n,:) = reshape(get_g_t_vector(thalamic_spikes(iter,thal_col,thal_n,:), length(tspan)) .* reshape(xe_thalamic(iter,thal_col,thal_n,:), 1, length(tspan)),  1,1,length(tspan));
        end
    end

    

    % simulation
    for i=6:length(tspan)

% 	fprintf("i = %d\n", i);
	for c=1:n_columns
	            
%         fprintf("\n +++++iter numm %d, column %d +++++++\n", iter, c);
            
		for n=1:n_total_neurons
					
			% voltage sum from excitatory neighbouring columns, will be useful for inhib and exc neurons
			% TODO - change column input neurons, if n is exc and if n is pv, som
            epsc_ex_neuron_back_c2 = 0;
            
			if c-2 >= 1
                 if n <= n_excitatory
                    for j=1:n_excitatory
					    spike_train_exc = spikes(iter,c-2,j,:);
                        g_t = get_g_t(spike_train_exc, dt, i-5, tspan,kernel_kt);
					    epsc_ex_neuron_back_c2 = epsc_ex_neuron_back_c2 + g_t*xe(iter, c-2,j,i-5)*across_columns_exc_to_exc_weight_matrix(iter,c-2,c,i,j,n); 
                    end
                 else % for inhibitory neurons
                        for j=1:n_excitatory
					        spike_train_exc = spikes(iter,c-2,j,:);
                            g_t = get_g_t(spike_train_exc, dt, i-5, tspan,kernel_kt);
					        epsc_ex_neuron_back_c2 = epsc_ex_neuron_back_c2 + g_t*xe(iter, c-2,j,i-5); 
                        end
                 end
			end

			epsc_ex_neuron_back_c1 = 0;
			if c-1 >= 1
                if n <= n_excitatory
                        for j=1:n_excitatory
					        spike_train_exc = spikes(iter,c-1,j,:);
					        g_t = get_g_t(spike_train_exc, dt, i-5, tspan,kernel_kt);
					        epsc_ex_neuron_back_c1 = epsc_ex_neuron_back_c1 + g_t*xe(iter,c-1,j,i-5)*across_columns_exc_to_exc_weight_matrix(iter,c-1,c,i,j,n); 
	                    end
                else
                        for j=1:n_excitatory
					            spike_train_exc = spikes(iter,c-1,j,:);
					            g_t = get_g_t(spike_train_exc, dt, i-5, tspan,kernel_kt);
					            epsc_ex_neuron_back_c1 = epsc_ex_neuron_back_c1 + g_t*xe(iter,c-1,j,i-5); 
			            end
                
                end
				
			end

			epsc_ex_neuron_front_c1 = 0;
			if c+1 <= n_columns
                if n <= n_excitatory
                        for j=1:n_excitatory
				            spike_train_exc = spikes(iter,c+1,j,:);
					        g_t = get_g_t(spike_train_exc, dt, i-5, tspan,kernel_kt);
					        epsc_ex_neuron_front_c1 = epsc_ex_neuron_front_c1 + g_t*xe(iter,c+1,j,i-5)*across_columns_exc_to_exc_weight_matrix(iter,c+1,c,i,j,n);
	                    end
                else
                        for j=1:n_excitatory
				            spike_train_exc = spikes(iter,c+1,j,:);
					        g_t = get_g_t(spike_train_exc, dt, i-5, tspan,kernel_kt);
					        epsc_ex_neuron_front_c1 = epsc_ex_neuron_front_c1 + g_t*xe(iter,c+1,j,i-5);
				        end
                end
			end	


			epsc_ex_neuron_front_c2 = 0;
			if c+2 <= n_columns
                if n <= n_excitatory
                        for j=1:n_excitatory
					        spike_train_exc = spikes(iter,c+2,j,:);
					        g_t = get_g_t(spike_train_exc, dt, i-5, tspan,kernel_kt);
					        epsc_ex_neuron_front_c2 = epsc_ex_neuron_front_c2 + g_t*xe(iter,c+2,j,i-5)*across_columns_exc_to_exc_weight_matrix(iter,c+2,c,i,j,n);
	                    end
                else
                        for j=1:n_excitatory
					        spike_train_exc = spikes(iter,c+2,j,:);
					        g_t = get_g_t(spike_train_exc, dt, i-5, tspan,kernel_kt);
					        epsc_ex_neuron_front_c2 = epsc_ex_neuron_front_c2 + g_t*xe(iter,c+2,j,i-5);
				        end
                end
			end	


			epsc_ex_own_column = 0;
			if n > n_excitatory % if n is inhibitory neuron
                for j=1:n_excitatory
                    if j == n
                        continue;
                    end
                    spike_train_exc = spikes(iter,c,j,:);
                    g_t = get_g_t(spike_train_exc, dt, i-5, tspan,kernel_kt);
                    x_e_presyn_neuron = xe(iter,c,j,i-5);
                    epsc_ex_own_column	= epsc_ex_own_column + g_t*x_e_presyn_neuron;
                end
            else
                for j=1:n_excitatory
                    if j == n
                        continue;
                    end
                    spike_train_exc = spikes(iter,c,j,:);
                    g_t = get_g_t(spike_train_exc, dt, i-5, tspan,kernel_kt);
                    x_e_presyn_neuron = xe(iter,c,j,i-5);
                    epsc_ex_own_column	= epsc_ex_own_column + g_t*x_e_presyn_neuron*exc_to_exc_weight_matrix(iter,c,i-5,j,n);
                end
            end

			% epsc from inhibitory neurons
			epsc_pv_own_column = 0; 
			for j=n_excitatory+1:n_excitatory + n_pv
				if j == n
					continue;
                end
                spike_train_inh = spikes(iter,c,j,:);
				g_t = get_g_t(spike_train_inh, dt, i-5, tspan,kernel_kt);
                epsc_pv_own_column = epsc_pv_own_column + g_t*xe(iter,c,j,i-5);
            end

            epsc_som_own_column = 0;
			for j=n_excitatory+n_pv+1:n_total_neurons
				if j == n
					continue;
                end
                spike_train_inh = spikes(iter,c,j,:);
				g_t = get_g_t(spike_train_inh, dt, i-5, tspan,kernel_kt);
                epsc_som_own_column = epsc_som_own_column + g_t*xe(iter,c,j,i-5);
            end

            epsc_som_back_c2 = 0;
            if c-2 >= 1
				for j=n_excitatory+n_pv+1:n_total_neurons
					spike_train_exc = spikes(iter,c-2,j,:);
                    g_t = get_g_t(spike_train_exc, dt, i-5, tspan,kernel_kt);
                    epsc_som_back_c2 = epsc_som_back_c2 + g_t*xe(iter, c-2,j,i-5);
                end
            end

            epsc_som_back_c1 = 0;
            if c-1 >= 1
				for j=n_excitatory+n_pv+1:n_total_neurons
					spike_train_exc = spikes(iter,c-1,j,:);
                    g_t = get_g_t(spike_train_exc, dt, i-5, tspan,kernel_kt);
                    epsc_som_back_c1 = epsc_som_back_c1 + g_t*xe(iter, c-1,j,i-5);
                end
            end

            epsc_som_front_c1 = 0;
            if c+1 <= n_columns
				for j=n_excitatory+n_pv+1:n_total_neurons
				    spike_train_exc = spikes(iter,c+1,j,:);
					g_t = get_g_t(spike_train_exc, dt, i-5, tspan,kernel_kt);
                    epsc_som_front_c1 = epsc_som_front_c1 + g_t*xe(iter,c+1,j,i-5);
				end
            end

            epsc_som_front_c2 = 0;
            if c+2 <= n_columns
				for j=n_excitatory+n_pv+1:n_total_neurons
					spike_train_exc = spikes(iter,c+2,j,:);
					g_t = get_g_t(spike_train_exc, dt, i-5, tspan,kernel_kt);
                    epsc_som_front_c2 = epsc_som_front_c2 + g_t*xe(iter,c+2,j,i-5);
				end
			end
		
          % epsc from thalamic neurons
          epsc_from_thalamic = 0;
          cols_giving_input = c:c+4;
          if n <= n_excitatory
               for col_index=1:length(cols_giving_input)
                    neuron_num = thalamic_connections(n,col_index);
                    epsc_from_thalamic = epsc_from_thalamic + epsc_thalamic(iter,cols_giving_input(col_index),neuron_num,i)*weight_thalamic_to_exc_l4_arr(col_index);
               end
          elseif n > n_excitatory && n <= n_excitatory + n_pv
             for col_index=1:length(cols_giving_input)
                    neuron_num = thalamic_connections(n,col_index);
                    epsc_from_thalamic = epsc_from_thalamic + epsc_thalamic(iter,cols_giving_input(col_index),neuron_num,i)*weight_thalamic_to_pv_l4_arr(col_index);
             end
          else
              for col_index=1:length(cols_giving_input)
                    neuron_num = thalamic_connections(n,col_index);
                    epsc_from_thalamic = epsc_from_thalamic + epsc_thalamic(iter,cols_giving_input(col_index),neuron_num,i)*weight_thalamic_to_som_l4_arr(col_index);
             end
          end

          thalamic_epsc_tensor(iter,c,n,i-5) = epsc_from_thalamic;
          
                    
		  if n <= n_excitatory % excitatory neuron
      			total_epsc = epsc_ex_neuron_back_c2 + ...
                    epsc_ex_neuron_back_c1 + ...
                    epsc_ex_neuron_front_c1 + ...
                    epsc_ex_neuron_front_c2 + ...
                    epsc_ex_own_column  + ... % weight already included in weight matrix
                    epsc_som_back_c2 * J_e_som_2 + ...
                    epsc_som_back_c1 * J_e_som_1 + ...
                    epsc_som_front_c1 * J_e_som_1 + ...
                    epsc_som_front_c2 * J_e_som_2 + ...
                    epsc_pv_own_column * J_e_pv + ...
                    epsc_som_own_column * J_e_som;
                recurrence_exc_self_column_epsc_tensor(iter,c,n,i-1) = epsc_ex_own_column;
                recurrence_exc_neighbour_column_epsc_tensor(iter,c,n,i-1) = epsc_ex_neuron_back_c2 + ...
                                                                            epsc_ex_neuron_back_c1 + ...
                                                                            epsc_ex_neuron_front_c1 + ...
                                                                            epsc_ex_neuron_front_c2 ;

                recurrence_inh_self_column_epsc_tensor(iter,c,n,i-1) = epsc_pv_own_column*J_e_pv + epsc_som_own_column*J_e_som;
                recurrence_inh_neighbour_column_epsc_tensor(iter,c,n,i-1) = epsc_som_back_c2 * J_e_som_2 + ...
                                                                            epsc_som_back_c1 * J_e_som_1 + ...
                                                                            epsc_som_front_c1 * J_e_som_1 + ...
                                                                            epsc_som_front_c2 * J_e_som_2 ;

              elseif n > n_excitatory && n <= n_excitatory + n_pv % pv neuron
                total_epsc = epsc_ex_neuron_back_c2 * J_pv_e_2 + ...
                    epsc_ex_neuron_back_c1 * J_pv_e_1 + ...
                    epsc_ex_neuron_front_c1 * J_pv_e_1 + ...
                    epsc_ex_neuron_front_c2 * J_pv_e_2 + ...
                    epsc_ex_own_column * J_pv_e_0 + ...
                    epsc_pv_own_column * J_pv_pv +...
                    epsc_som_own_column * J_pv_som;

                recurrence_exc_self_column_epsc_tensor(iter,c,n,i-1) = epsc_ex_own_column * J_pv_e_0;
                recurrence_exc_neighbour_column_epsc_tensor(iter,c,n,i-1) = epsc_ex_neuron_back_c2 * J_pv_e_2 + ...
                                                                            epsc_ex_neuron_back_c1 * J_pv_e_1 + ...
                                                                            epsc_ex_neuron_front_c1 * J_pv_e_1 + ...
                                                                            epsc_ex_neuron_front_c2 * J_pv_e_2;
                recurrence_inh_self_column_epsc_tensor(iter,c,n,i-1) = epsc_pv_own_column * J_pv_pv + epsc_som_own_column * J_pv_som;
                recurrence_inh_neighbour_column_epsc_tensor(iter,c,n,i-1) = 0;
              
              else % som neuron

      			total_epsc = epsc_ex_neuron_back_c2 * J_som_e_2 + ...
                    epsc_ex_neuron_back_c1 * J_som_e_1 + ...
                    epsc_ex_neuron_front_c1 * J_som_e_1 + ...
                    epsc_ex_neuron_front_c2 * J_som_e_2 + ...
                    epsc_ex_own_column * J_som_e_0 + ...
                    epsc_pv_own_column * J_som_pv +...
                    epsc_som_own_column * J_som_som;
                recurrence_exc_self_column_epsc_tensor(iter,c,n,i-1) = epsc_ex_own_column * J_som_e_0;
                recurrence_exc_neighbour_column_epsc_tensor(iter,c,n,i-1) = epsc_ex_neuron_back_c2 * J_som_e_2 + ...
                                                                            epsc_ex_neuron_back_c1 * J_som_e_1 + ...
                                                                            epsc_ex_neuron_front_c1 * J_som_e_1 + ...
                                                                            epsc_ex_neuron_front_c2 * J_som_e_2;
                recurrence_inh_self_column_epsc_tensor(iter,c,n,i-1) = epsc_pv_own_column * J_som_pv + epsc_som_own_column * J_som_som;
                recurrence_inh_neighbour_column_epsc_tensor(iter,c,n,i-1) = 0;
          end

          total_epsc = total_epsc + epsc_from_thalamic; % recurrence + thalamic
            % clip test - to see whether the later spike(s) is due to
            % params or really disihibition
            % uncomment and see if spikes comes or not
 
%                 if total_epsc < 0
%                     total_epsc = 0;
%                 end
           
			
                % I_background = rand * (1);
                % to see only effect of thalamic in feedforward
                
                
%                 r = normrnd(5,15);

                 I_background = background_epsc(iter,col,n,i);   
                 I_background_tensor(iter, c,n,i) = I_background;
            
            % calculate voltage using the function
% 			[voltages(iter,c, n, i), u_values(iter,c, n, i)] = calculate_v_u(v_current, u_current, dt, neuron_params_rb_ss, total_epsc, I_background );
		   
            spike_vec = squeeze(spikes(iter,c,n,:));
            latest_spike_time = -1;
            for s=i-1:-1:i-5
                if s >= 1 && spike_vec(s) == 1
                    latest_spike_time = s;
                    break;
                end
            end
            [voltages(iter,c, n, i), i1_tensor(iter,c, n, i), i2_tensor(iter,c, n, i), theta_tensor(iter,c, n, i), spikes(iter,c,n,i)] = calculate_new_state_dynamic_threshold_rule(voltages(iter,c, n, i-1), i1_tensor(iter,c, n, i-1), i2_tensor(iter,c, n, i-1), theta_tensor(iter,c, n, i-1), total_epsc, I_background,dt,i,latest_spike_time);
          
            if spikes(iter,c,n,i-1) == 1
				spikes(iter,c,n,i) = 0;
 			end

			M = 0;
			if spikes(iter,c,n,i) == 1
				M = 1;
            end

        
		
            
            
            %	fprintf("voltage returned from function is %f \n", voltages(c,n,i));
    
            % update synaptic resources
            
            
            current_xr = xr(iter,c, n, i-1);
            current_xe = xe(iter,c, n, i-1);
            current_xi = xi(iter,c, n, i-1);
            
            % only excitatory synpases are depressing, 
            % inhibitory synapses are non depressing
            if n <= n_excitatory
                new_xr = update_xr(M, current_xr, current_xi, tau_re, tau_ir);
                new_xe = update_xe(M, current_xr, current_xe, tau_re, tau_ei);
    
                % clipping for computational faults
                if new_xr < 0
                    new_xr = 0;
                end
                if new_xr > 1
                    new_xr = 1;
                end

                if new_xe < 0
                    new_xe = 0;
                end
                if new_xe > 1
                    new_xe = 1;
                end

                xr(iter,c,n,i) = new_xr;
                xe(iter,c, n, i) = new_xe;
                xi(iter,c, n, i) = 1 - (new_xr + new_xe);
            end


        end
    
   
%     fprintf("xr %f, xe %f, xi %f\n", xr(c,n,i),xe(c,n,i), xi(c,n,i));
   % pause(0.4);
    
    end

	% ------updating weights using STDP
        
                % if change has occured in LTP or LTD for a synapse, 
                % Don't undo it in the next iteration
        for col_stdp=1:n_columns   
            either_LTP_or_LTD_occured = zeros(n_excitatory, n_excitatory);
            for N=1:n_excitatory
                % N -> postyn : LTD
                for postsyn_neuron=1:n_excitatory
                    if spikes(iter,col_stdp,N,i) == 1 % if there is a spike
                        found_spike_in_window_LTD = 0;

                        for postsyn_spike_time=i-1:-1:i-20
                            if postsyn_spike_time >= 1 && spikes(iter,col_stdp,postsyn_neuron,postsyn_spike_time) == 1 && spikes(iter,col_stdp,postsyn_neuron,i) == 0 
                                exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron) = exc_to_exc_weight_matrix(iter,col_stdp,i-1,N,postsyn_neuron)*(1 - Amp_weak*exp(-abs(i-postsyn_spike_time)/tau_weak));
                                % clipping weights
                                if exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron) < minimum_weight_exc_to_exc
                                    exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron) = minimum_weight_exc_to_exc;
                                end
                                if exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron) > maximum_weight_exc_to_exc
                                    exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron) = maximum_weight_exc_to_exc;
                                end
                                num_of_LTDs(iter,col_stdp,i) = num_of_LTDs(iter,col_stdp,i) + 1;
                                
%                                 if N == 5 && (postsyn_neuron == 13 || postsyn_neuron == 12) && col_stdp == 2
%                                     fprintf("\n postsyn_neuron is %d \n", postsyn_neuron)
%                                     fprintf("\n LTD - !!! - old %f, new %f \n ",exc_to_exc_weight_matrix(iter,col_stdp,i-1,N,postsyn_neuron),exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron))
%                                     fprintf("\n due to same time ??? %d %d \n", spikes(iter,col_stdp,N,i), spikes(iter,col_stdp,postsyn_neuron,i))
%                                     pause(1)
%                                 end
                                either_LTP_or_LTD_occured(N,postsyn_neuron) = 1;
                                found_spike_in_window_LTD = 1;
                                break;
                            end
                        end

                        if found_spike_in_window_LTD == 0 && either_LTP_or_LTD_occured(N,postsyn_neuron) == 0
                            exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron) = exc_to_exc_weight_matrix(iter,col_stdp,i-1,N,postsyn_neuron);
%                             if N == 5 && (postsyn_neuron == 13 || postsyn_neuron == 12) && col_stdp == 2
%                                     fprintf("\n postsyn_neuron is %d \n", postsyn_neuron)
%                                     fprintf("\n LTD - @@@ - old %f, new %f \n ",exc_to_exc_weight_matrix(iter,col_stdp,i-1,N,postsyn_neuron),exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron))
%                                     pause(1)
%                             end
                        end
                    else % if there is no spike

                        if either_LTP_or_LTD_occured(N,postsyn_neuron) == 0
                               exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron) = exc_to_exc_weight_matrix(iter,col_stdp,i-1,N,postsyn_neuron);
%                                 if N == 5 && (postsyn_neuron == 13 || postsyn_neuron == 12) && col_stdp == 2
%                                     fprintf("\n postsyn_neuron is %d \n", postsyn_neuron)
%                                     fprintf("\n LTD - $$$ - old %f, new %f \n ",exc_to_exc_weight_matrix(iter,col_stdp,i-1,N,postsyn_neuron),exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron))
%                                     pause(1)
%                                 end
                        end
                        
                    end
                end


                % presyn -> N : LTP
                for presyn_neuron=1:n_excitatory
                    if spikes(iter,col_stdp,N,i) == 1 % if there is a spike
                        found_spike_in_window_LTP = 0;

                        for presyn_spike_time=i-1:-1:i-20
                            if presyn_spike_time >= 1 && spikes(iter,col_stdp,presyn_neuron,presyn_spike_time) == 1 && spikes(iter,col_stdp,presyn_neuron,i) == 0
                                    exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N) = exc_to_exc_weight_matrix(iter,col_stdp,i-1,presyn_neuron,N)*(1 + Amp_strength*exp(-abs(i-presyn_spike_time)/tau_strength));
                                    % clipping weights 
                                    if exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N) < minimum_weight_exc_to_exc
                                        exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N) = minimum_weight_exc_to_exc;
                                    end
                                    if exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N) > maximum_weight_exc_to_exc
                                        exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N) = maximum_weight_exc_to_exc;
                                    end

                                    num_of_LTPs(iter,col_stdp,i) = num_of_LTPs(iter,col_stdp,i) + 1; 
%                                     if presyn_neuron == 5 && (N == 13 || N == 12) && col_stdp == 2
%                                         fprintf("\n postsyn_neuron is %d \n", N)
%                                         fprintf("\n LTP - !!! - old %f, new %f \n ",exc_to_exc_weight_matrix(iter,col_stdp,i-1,presyn_neuron,N),exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N))
%                                         pause(1)
%                                     end
                                
                                    
                                
                                either_LTP_or_LTD_occured(presyn_neuron,N) = 1;
                                found_spike_in_window_LTP = 1;
                                break
                            end
                        end

                        if found_spike_in_window_LTP == 0 && either_LTP_or_LTD_occured(presyn_neuron,N) == 0
                            exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N) = exc_to_exc_weight_matrix(iter,col_stdp,i-1,presyn_neuron,N);
%                             if presyn_neuron == 5 && (N == 13 || N == 12) && col_stdp == 2
%                                     fprintf("\n postsyn_neuron is %d \n", N)
%                                     fprintf("\n LTP - @@@ - old %f, new %f \n ",exc_to_exc_weight_matrix(iter,col_stdp,i-1,presyn_neuron,N),exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N))
%                                     pause(1)
%                             end
                        end
                    else % if there is no spike
                        if either_LTP_or_LTD_occured(presyn_neuron,N) == 0
                            exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N) = exc_to_exc_weight_matrix(iter,col_stdp,i-1,presyn_neuron,N);
%                             if presyn_neuron == 5 && (N == 13 || N == 12)  && col_stdp == 2
%                                     fprintf("\n postsyn_neuron is %d \n", N)
%                                     fprintf("\n LTP - $$$ - old %f, new %f \n ",exc_to_exc_weight_matrix(iter,col_stdp,i-1,presyn_neuron,N),exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N))
%                                     pause(1)
%                             end
                        end
                     end
                end

                

            end % end of for all excitatory neurons

        end % end of stdp within column for all columns

        
        % re-initialize at the end of token
        if ismember(i,token_start_times)
%             disp("*******************token resest***************************")
            voltages(:,:,:,i) = v0;  
            xr(:, :, :, i) = 1;
            xe(:, :, :, i) = 0;
            xi(:, :, :, i) = 0;
            i1_tensor(:, :, :, i) = 0.01;
            i2_tensor(:, :, :, i) = 0.001;
            theta_tensor(:, :, :,i) = -50.0;
        end
%	break % for testing only one iteration
    
    
    either_LTP_or_LTD_occured = zeros(n_columns, n_columns, n_excitatory, n_excitatory);% pre, post
    % plasticity across columns
    for cc=1:num_connected_pairs
        c1 = connected_columns_arr(cc,1);
        c2 = connected_columns_arr(cc,2);
       
        for N=1:n_excitatory
           % N -> postsyn : LTD
           for postsyn_neuron=1:n_excitatory
                if spikes(iter,c1,N,i) == 1
                    found_spike_in_window_LTD = 0;
                    for postsyn_spike_time=i-1:-1:i-20
                        if postsyn_spike_time >= 1 && spikes(iter,c2,postsyn_neuron,i) == 0 && spikes(iter,c2,postsyn_neuron,postsyn_spike_time) == 1
                            across_columns_exc_to_exc_weight_matrix(iter,c1,c2,i,N,postsyn_neuron) = across_columns_exc_to_exc_weight_matrix(iter,c1,c2,i-1,N,postsyn_neuron)*(1 - Amp_weak*exp(-abs(i-postsyn_spike_time)/tau_weak)); 
                            either_LTP_or_LTD_occured(c1,c2,N,postsyn_neuron) = 1;
                            found_spike_in_window_LTD = 1;
                        end
                    end

                    if found_spike_in_window_LTD == 0 && either_LTP_or_LTD_occured(c1,c2,N,postsyn_neuron) == 0
                        across_columns_exc_to_exc_weight_matrix(iter,c1,c2,i,N,postsyn_neuron) =  across_columns_exc_to_exc_weight_matrix(iter,c1,c2,i-1,N,postsyn_neuron);
                    end
                else % if no spike
                    across_columns_exc_to_exc_weight_matrix(iter,c1,c2,i,N,postsyn_neuron) =  across_columns_exc_to_exc_weight_matrix(iter,c1,c2,i-1,N,postsyn_neuron);
                end
           end

            % presyn -> N : LTP
           for presyn_neuron=1:n_excitatory
                if spikes(iter,c1,N,i) == 1

                    found_spike_in_window_LTP = 0;
                    for presyn_spike_time=i-1:-1:i-20
                        if presyn_spike_time >= 1 && spikes(iter,c2,presyn_neuron,i) == 0 && spikes(iter,c2,presyn_neuron,presyn_spike_time) == 1
                            across_columns_exc_to_exc_weight_matrix(iter,c2,c1,i,presyn_neuron,N) = across_columns_exc_to_exc_weight_matrix(iter,c2,c1,i-1,presyn_neuron,N)*(1 + Amp_strength*exp(-abs(i-presyn_spike_time)/tau_strength));
                            found_spike_in_window_LTP = 1;
                            either_LTP_or_LTD_occured(c2,c1,presyn_neuron,N) = 1;
                        end
                    end

                    if found_spike_in_window_LTP == 0 && either_LTP_or_LTD_occured(c2,c1,presyn_neuron,N) == 0
                        across_columns_exc_to_exc_weight_matrix(iter,c2,c1,i,presyn_neuron,N) = across_columns_exc_to_exc_weight_matrix(iter,c2,c1,i-1,presyn_neuron,N);
                    end
                else % if no spike
                    across_columns_exc_to_exc_weight_matrix(iter,c2,c1,i,presyn_neuron,N) = across_columns_exc_to_exc_weight_matrix(iter,c2,c1,i-1,presyn_neuron,N);
                end
           end

          
            
        end % end of excitatory neurons in c1
    end % end of connected pairs
    

    end % end of a single time step

end % end of an iter

save('batch_1.mat')
toc