clear all;
close all;


n_iters = 1;

% basic variables;
n_columns = 5;
n_excitatory = 20; 
n_inhibitory = 5; 
n_total_neurons = n_excitatory + n_inhibitory;
n_thalamic = 9; num_of_input_giving_thalamic = 4;
    
% time step
physical_time_in_ms = 1; %dt time step 
dt = 1;  % 0.2 dt = 20 ms, so 0.01 = 1 ms 
t_simulate = 2780; 
tspan = 0:dt:t_simulate;

% making bins of 100ms = 20*dt and calculating spike rate
spike_rate_dt = 20*dt;
spike_rate_length = (length(tspan)-1)/(spike_rate_dt/dt);


% connection strength
weight_reducing_l4 = 0.3; % for now all weights reduced by factor of 0.2
J_ee_0 = 6*weight_reducing_l4; 
J_ie_0 = 0.5*weight_reducing_l4;
J_ei = -4*weight_reducing_l4; 
J_ii = -0.5*weight_reducing_l4;
J_ee_1 = 0.045*weight_reducing_l4; 
J_ie_1 = 0.0035*weight_reducing_l4; 
J_ee_2 = 0.015*weight_reducing_l4; 
J_ie_2 = 0.0015*weight_reducing_l4;

% voltages and terms from it are 3d tensors
voltages = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
i1_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
i2_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
theta_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
spikes = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
spike_rates = zeros(n_iters, n_columns, n_total_neurons, spike_rate_length);
recurrence_epsc_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan)-1);
epsc_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan)-1);
epsc_with_background_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan)-1);
I_background_tensor = zeros(n_iters, length(tspan));
feedforward_epsc_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan)-1);

% synaptic resources
xr = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
xe = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
xi = zeros(n_iters, n_columns, n_total_neurons, length(tspan));

% mapping from input thalamic neurons to a1 column neurons
all_combinations = nchoosek(1:n_thalamic, num_of_input_giving_thalamic);
mapping_matrix_thalamic_to_a1 = all_combinations(1:n_total_neurons,:);


thalamic_poisson_spikes_std = zeros(n_iters, n_thalamic, length(tspan));
thalamic_poisson_spikes_dev = zeros(n_iters, n_thalamic, length(tspan));
thalamic_poisson_spikes_common = zeros(n_iters, n_thalamic, length(tspan));

% protochol basic params
pre_stimulus_time = 500; post_stimulus_time = 480; 
single_stimulus_duration = 100; gap_duration = 20;
n_tokens = 15;
lamda_i = 4;

prob = 0.9;

freq_stim_in_std = 200;
unfreq_stim_in_std = 20;

freq_stim_in_dev = unfreq_stim_in_std;
unfreq_stim_in_dev = freq_stim_in_std;

lamda_pre_stimulus = lamda_i*ones(1, pre_stimulus_time);
lamda_post_stimulus = lamda_i*ones(1, post_stimulus_time);

% --- standard protochol -----
lamda_std = zeros(1, n_tokens*(single_stimulus_duration + gap_duration));

for tok=1:n_tokens
   
    if randi([1,100]) <= prob*100
        stim = freq_stim_in_std;
    else
        stim = unfreq_stim_in_std;
    end

    ind = (tok-1)*(single_stimulus_duration+gap_duration) + 1;
    
    for i=ind:ind+99
        lamda_std(1,i) = stim;
    end
    for i=ind+100:ind+119
        lamda_std(1,i) = lamda_i;
    end
end

% merging pre-stim, stim, post sim
lamda_std_protochol = [];
lamda_std_protochol = [lamda_std_protochol, lamda_pre_stimulus];
lamda_std_protochol = [lamda_std_protochol, lamda_std];
lamda_std_protochol = [lamda_std_protochol, lamda_post_stimulus];
lamda_std_protochol = [lamda_std_protochol, lamda_i];

% --- deviant protochol -----
lamda_dev = zeros(1, n_tokens*(single_stimulus_duration + gap_duration));

for tok=1:n_tokens
   
        if randi([1,100]) <= prob*100
            stim = freq_stim_in_dev;
        else
            stim = unfreq_stim_in_dev;
        end
    
        ind = (tok-1)*(single_stimulus_duration+gap_duration) + 1;
        
        for i=ind:ind+99
            lamda_dev(1,i) = stim;
        end
        for i=ind+100:ind+119
            lamda_dev(1,i) = lamda_i;
        end
    end
    
    % merging pre-stim, stim, post sim
    lamda_dev_protochol = [];
    lamda_dev_protochol = [lamda_dev_protochol, lamda_pre_stimulus];
    lamda_dev_protochol = [lamda_dev_protochol, lamda_dev];
    lamda_dev_protochol = [lamda_dev_protochol, lamda_post_stimulus];
    lamda_dev_protochol = [lamda_dev_protochol, lamda_i];


% ---- common protochol ---- thalamic poisson input other than std and dev
lamda_common = lamda_i*ones(1, length(tspan));

% calculating epsc of each thalamic neuron
weight_thalamic_to_a1 = 80; xe_thalamic = 1;
epsc_thalamic_std = zeros(n_iters,n_thalamic, length(tspan));
epsc_thalamic_dev = zeros(n_iters,n_thalamic, length(tspan));
epsc_thalamic_common = zeros(n_iters,n_thalamic, length(tspan));

%% time constant for synaptic resources
tau_re = 0.9; tau_ir = 5000; tau_ei = 27;
    
% initialize
v0 = -70;  
xr(:, :, :, 1) = 1;
voltages(:, :, :, 1) = v0; % 
i1_tensor(:, :, :, 1) = 0.01;
i2_tensor(:, :, :, 1) = 0.001;
theta_tensor(:, :, :, 1) = -50.0;

for iter=1:n_iters
    
    fprintf("\n ------iter numm %d -----\n", iter);

    % std - thalamic
    for i=1:n_thalamic
         thalamic_poisson_spikes_std(iter,i, :) = reshape(poisson_generator(lamda_std_protochol, dt), 1, 1, length(tspan));
    end
    for i=1:n_thalamic
        epsc_thalamic_std(iter,i,:) = reshape(get_g_t_vector(thalamic_poisson_spikes_std(iter,i,:), length(tspan)) * weight_thalamic_to_a1 * xe_thalamic,  1,1,length(tspan));
    end

    % dev - thalamic
    for i=1:n_thalamic
         thalamic_poisson_spikes_dev(iter,i, :) = reshape(poisson_generator(lamda_dev_protochol, dt), 1, 1, length(tspan));
    end
    for i=1:n_thalamic
        epsc_thalamic_dev(iter,i,:) = reshape(get_g_t_vector(thalamic_poisson_spikes_dev(iter,i,:), length(tspan)) * weight_thalamic_to_a1 * xe_thalamic,  1,1,length(tspan));
    end

    % common - thalamic
    for i=1:n_thalamic
         thalamic_poisson_spikes_common(iter,i, :) = reshape(poisson_generator(lamda_common, dt), 1, 1, length(tspan));
    end
    for i=1:n_thalamic
        epsc_thalamic_common(iter,i,:) = reshape(get_g_t_vector(thalamic_poisson_spikes_common(iter,i,:), length(tspan)) * weight_thalamic_to_a1 * xe_thalamic,  1,1,length(tspan));
    end

    % simulation
    for i=2:length(tspan)
	for c=1:n_columns
	     fprintf("\n +++++iter numm %d, column %d +++++++\n", iter, c);
        if c == 2
                epsc_thalamic = epsc_thalamic_std;
        elseif c == 4
                epsc_thalamic = epsc_thalamic_dev;
        else
                epsc_thalamic = epsc_thalamic_common;
        end

		for n=1:n_total_neurons
					
			% voltage sum from excitatory neighbouring columns, will be useful for inhib and exc neurons
			epsc_ex_neuron_back_c2 = 0;
			if c-2 >= 1
				for j=1:n_excitatory
					spike_train_exc = spikes(iter,c-2,j,:);
                    g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
					epsc_ex_neuron_back_c2 = epsc_ex_neuron_back_c2 + g_t*xe(iter, c-2,j,i-1); 
                 end
			end

			epsc_ex_neuron_back_c1 = 0;
			if c-1 >= 1
				for j=1:n_excitatory
					spike_train_exc = spikes(iter,c-1,j,:);
					g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
					epsc_ex_neuron_back_c1 = epsc_ex_neuron_back_c1 + g_t*xe(iter,c-1,j,i-1); 
				end
			end

			epsc_ex_neuron_front_c1 = 0;
			if c+1 <= n_columns
				for j=1:n_excitatory
				    spike_train_exc = spikes(iter,c+1,j,:);
					g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
					epsc_ex_neuron_front_c1 = epsc_ex_neuron_front_c1 + g_t*xe(iter,c+1,j,i-1);
				end
			end	


			epsc_ex_neuron_front_c2 = 0;
			if c+2 <= n_columns
				for j=1:n_excitatory
					spike_train_exc = spikes(iter,c+2,j,:);
					g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
					epsc_ex_neuron_front_c2 = epsc_ex_neuron_front_c2 + g_t*xe(iter,c+2,j,i-1);
				end
			end	


			epsc_ex_own_column = 0;
			for j=1:n_excitatory
				if j == n
					continue;
				end
				spike_train_exc = spikes(iter,c,j,:);
				g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
				epsc_ex_own_column	= epsc_ex_own_column + g_t*xe(iter,c,j,i-1); 
			end

			epsc_inh_own_column = 0;
			for j=n_excitatory+1:n_total_neurons
				if j == n
					continue;
                end
                spike_train_inh = spikes(iter,c,j,:);
				g_t = get_g_t(spike_train_inh, dt, i-1, tspan);
				epsc_inh_own_column = epsc_inh_own_column + g_t*xe(iter,c,j,i-1);
            end
		
          % epsc from thalamic neurons
          epsc_from_thalamic = 0;
          for ttt=1:num_of_input_giving_thalamic
             thalamic_neuron_num = mapping_matrix_thalamic_to_a1(n, ttt);
             epsc_from_thalamic = epsc_from_thalamic + epsc_thalamic(iter,thalamic_neuron_num, i-1);
          end
          
          
		  if n <= n_excitatory
			total_epsc = epsc_ex_neuron_back_c2 * J_ee_2 + ...
						epsc_ex_neuron_back_c1 * J_ee_1 + ...
						epsc_ex_neuron_front_c1 * J_ee_1 + ...
						epsc_ex_neuron_front_c2 * J_ee_2 + ...
						epsc_ex_own_column * J_ee_0 + ...
						epsc_inh_own_column * J_ei;
					
		  else

			total_epsc = epsc_ex_neuron_back_c2 * J_ie_2 + ...
						epsc_ex_neuron_back_c1 * J_ie_1 + ...
						epsc_ex_neuron_front_c1 * J_ie_1 + ...
						epsc_ex_neuron_front_c2 * J_ie_2 + ...
						epsc_ex_own_column * J_ie_0 + ...
						epsc_inh_own_column * J_ii;
          end
          recurrence_epsc_tensor(iter, c, n, i-1) = total_epsc ;
          total_epsc = total_epsc + epsc_from_thalamic; % recurrence + thalamic
          % a temporary statement to check if backcurrent is strong enough
          % to produce current
%           total_epsc = 0;
          
          epsc_tensor(iter, c, n, i-1) = total_epsc;
            
          
% 		  fprintf("c,n - %d %d total espc %f \n", c,n,total_epsc);

			
%             I_background = rand * (1);
                % to see only effect of thalamic in feedforward
                I_background = 0;
%             if i<501 | i >1500
%                 I_background = 0;
%             end
            I_background_tensor(iter, i) = I_background;
            
            feedforward_epsc_tensor(iter,c,n,i-1) = I_background + epsc_from_thalamic;
          epsc_with_background_tensor(iter, c, n, i-1) = total_epsc + I_background;
    
            % calculate voltage using the function
% 			[voltages(iter,c, n, i), u_values(iter,c, n, i)] = calculate_v_u(v_current, u_current, dt, neuron_params_rb_ss, total_epsc, I_background );
		   
            [voltages(iter,c, n, i), i1_tensor(iter,c, n, i), i2_tensor(iter,c, n, i), theta_tensor(iter,c, n, i), spikes(iter,c,n,i)] = calculate_new_state(voltages(iter,c, n, i-1), i1_tensor(iter,c, n, i-1), i2_tensor(iter,c, n, i-1), theta_tensor(iter,c, n, i-1), total_epsc, I_background,dt);
          
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

            xr(iter,c,n,i) = update_xr(M, current_xr, current_xi, tau_re, tau_ir);
            xe(iter,c, n, i) = update_xe(M, current_xr, current_xe, tau_re, tau_ei);
            xi(iter,c, n, i) = update_xi(current_xe, current_xi, tau_ei, tau_ir);
        
        end
    
   
%     fprintf("xr %f, xe %f, xi %f\n", xr(c,n,i),xe(c,n,i), xi(c,n,i));
   % pause(0.4);
    
    end

	
%	break % for testing only one iteration
end

end

save('l4_only1iter.mat');
