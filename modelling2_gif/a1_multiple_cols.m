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
dt = 1;  % 1 dt = 1 ms 
t_simulate = 4000; 
tspan = 0:dt:t_simulate;

% making bins
spike_rate_dt = 1*dt;
spike_rate_length = (length(tspan)-1)/(spike_rate_dt/dt);


% connection strength
weight_reducing_l4 = 0.25; % for now all weights reduced by factor of 0.2
increase_inhibitory_factor = 75;
weight_exc_factor = 15;
exc_to_exc_factor = 6;
inh_to_exc_factor = 1;
inh_to_inh_factor = 1;

J_ee_0 = 6*weight_reducing_l4*weight_exc_factor*exc_to_exc_factor; 
J_ie_0 = 0.5*weight_reducing_l4*weight_exc_factor;
J_ei = -4*weight_reducing_l4*increase_inhibitory_factor*inh_to_exc_factor; 
J_ii = -0.5*weight_reducing_l4*increase_inhibitory_factor*inh_to_inh_factor;
J_ee_1 = 0.045*weight_reducing_l4*weight_exc_factor; 
J_ie_1 = 0.0035*weight_reducing_l4*weight_exc_factor; 
J_ee_2 = 0.015*weight_reducing_l4*weight_exc_factor; 
J_ie_2 = 0.0015*weight_reducing_l4*weight_exc_factor;

% voltages and terms from it are 3d tensors
voltages = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
i1_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
i2_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
theta_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
spikes = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
spike_rates = zeros(n_iters, n_columns, n_total_neurons, spike_rate_length);


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
xr_thalamic_std = zeros(n_iters, n_thalamic, length(tspan));
xe_thalamic_std = zeros(n_iters, n_thalamic, length(tspan));
xi_thalamic_std = zeros(n_iters, n_thalamic, length(tspan));

xr_thalamic_dev = zeros(n_iters, n_thalamic, length(tspan));
xe_thalamic_dev = zeros(n_iters, n_thalamic, length(tspan));
xi_thalamic_dev = zeros(n_iters, n_thalamic, length(tspan));

xr_thalamic_common = zeros(n_iters, n_thalamic, length(tspan));
xe_thalamic_common = zeros(n_iters, n_thalamic, length(tspan));
xi_thalamic_common = zeros(n_iters, n_thalamic, length(tspan));

% mapping from input thalamic neurons to a1 column neurons
all_combinations = nchoosek(1:n_thalamic, num_of_input_giving_thalamic);
mapping_matrix_thalamic_to_a1 = all_combinations(1:n_total_neurons,:);

% protochol
thalamic_poisson_spikes_std = zeros(n_iters, n_thalamic, length(tspan));
thalamic_poisson_spikes_dev = zeros(n_iters, n_thalamic, length(tspan));
thalamic_poisson_spikes_common = zeros(n_iters, n_thalamic, length(tspan));

pre_stimulus_time = 500; post_stimulus_time = 500; 
single_stimulus_duration = 100; gap_duration = 200;

n_tokens = 10;
lamda_i = 0;

prob = 0.8;

freq_stim_in_std = 300;
unfreq_stim_in_std = 1;

freq_stim_in_dev = unfreq_stim_in_std;
unfreq_stim_in_dev = freq_stim_in_std;

lamda_pre_stimulus = lamda_i*ones(1, pre_stimulus_time);
lamda_post_stimulus = lamda_i*ones(1, post_stimulus_time);

% ---  protochols -----
lamda_std = zeros(1, n_tokens*(single_stimulus_duration + gap_duration));
lamda_dev = zeros(1, n_tokens*(single_stimulus_duration + gap_duration));

for tok=1:n_tokens
   
    if randi([1,100]) <= prob*100
        stim_s = freq_stim_in_std;
        stim_d = freq_stim_in_dev;
    else
        stim_s = unfreq_stim_in_std;
        stim_d = unfreq_stim_in_dev;
    end

    ind = (tok-1)*(single_stimulus_duration+gap_duration) + 1;
    
    for i=ind:ind+single_stimulus_duration-1
        lamda_std(1,i) = stim_s;
        lamda_dev(1,i) = stim_d;
    end
    for i=ind+single_stimulus_duration:ind+single_stimulus_duration+gap_duration-1
        lamda_std(1,i) = lamda_i;
        lamda_dev(1,i) = lamda_i;
    end
end

% ----merging pre-stim, stim, post sim
% std
lamda_std_protochol = [];
lamda_std_protochol = [lamda_std_protochol, lamda_pre_stimulus];
lamda_std_protochol = [lamda_std_protochol, lamda_std];
lamda_std_protochol = [lamda_std_protochol, lamda_post_stimulus];
lamda_std_protochol = [lamda_std_protochol, lamda_i]; % to make it n000 to n001

% dev
% merging pre-stim, stim, post sim
lamda_dev_protochol = [];
lamda_dev_protochol = [lamda_dev_protochol, lamda_pre_stimulus];
lamda_dev_protochol = [lamda_dev_protochol, lamda_dev];
lamda_dev_protochol = [lamda_dev_protochol, lamda_post_stimulus];
lamda_dev_protochol = [lamda_dev_protochol, lamda_i]; % to make it n000 to n001

% ---- common protochol ---- thalamic poisson input other than std and dev
lamda_common_protochol = lamda_i*ones(1, length(tspan));


% calculating epsc of each thalamic neuron
epsc_thalamic_std = zeros(n_iters,n_thalamic, length(tspan));
epsc_thalamic_dev = zeros(n_iters,n_thalamic, length(tspan));
epsc_thalamic_common = zeros(n_iters,n_thalamic, length(tspan));

%% generate inhomo poisson spikes for thalamic neurons
thalamic_poisson_spikes = zeros(n_iters, n_thalamic, length(tspan));
lamda = zeros(1, length(tspan));

% %% generate inhomo poisson spikes for thalamic neurons
% thalamic_poisson_spikes = zeros(n_iters, n_thalamic, length(tspan));
% lamda = zeros(1, length(tspan));
% % 100ms - 3-4 spikes, 200ms - 18-20 spikes, 300 - rest - 3-4 spikes
% % WARNING: FOR NOW THIS STIMULS IS HARD CODED, need to adjust acc to
% % t_simulate
% lamda_s = 300; lamda_i = 0;
% for i=1:500
%     lamda(1,i) = lamda_i;
% end
% for i=500:600
%     lamda(1,i) = lamda_s+lamda_i;
% end
% for i=601:length(tspan)
%     lamda(1,i) = lamda_i;
% end

% calculating epsc of each  thalamic neuron
epsc_thalamic = zeros(n_iters,n_thalamic, length(tspan));

weight_thalamic_to_exc_l4 = 550;
weight_thalamic_to_inh_l4 = 900;

%% time constant for synaptic resources
tau_re = 0.6; tau_ir = 700; tau_ei = 15;
tau_re_thalamic = 0.6; tau_ir_thalamic = 2700; tau_ei_thalamic = 35;

% izhikevich neuron params
% for rebound burst and sustained_spike
neuron_params_rb_ss = containers.Map({'a', 'b', 'c', 'd'}, [0.02 0.25 -70 0]); 
% a=0.02; b=0.25; c=-70;  d=0; % this should also cause disinhibition
% a=0.02; b=0.25; c=-55;  d=0.05; % used this on party day in lab
% neuron_params_rb_ss = containers.Map({'a', 'b', 'c', 'd'}, [0.03 0.25 -48 0]); 

% for rebound burst and phasic spike
neuron_params_rb_ps = containers.Map({'a', 'b', 'c', 'd'}, [0.02 0.25 -58 0.5]);
    
% initialize
v0 = -70;  
xr(:, :, :, 1:5) = 1; 
xr_thalamic_std(:,:,1) = 1;
xr_thalamic_dev(:,:,1) = 1;
xr_thalamic_common(:,:,1) = 1;
voltages(:, :, :, 1:5) = v0; % 
i1_tensor(:, :, :, 1:5) = 0.01;
i2_tensor(:, :, :, 1:5) = 0.001;
theta_tensor(:, :, :, 1:5) = -50.0;

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
    
    fprintf("------iter numm %d -----", iter);

    % TODO - seperate them as standard, deviant and common
    % check connections

    % ---- thalamic epsc -- std ---------------
    % thalamic
    for i=1:n_thalamic
         thalamic_poisson_spikes_std(iter,i, :) = reshape(poisson_generator(lamda_std_protochol, dt), 1, 1, length(tspan));
    end
    
    % depressing synapse thalamus -> l4
    for n_thal=1:n_thalamic
        for t=2:length(tspan)
                current_xr_thalamic_std = xr_thalamic_std(iter, n_thal, t-1);
                current_xe_thalamic_std = xe_thalamic_std(iter, n_thal, t-1);
                current_xi_thalamic_std = xi_thalamic_std(iter, n_thal, t-1);
                
                M = 0;
                if thalamic_poisson_spikes_std(iter,n_thal,t) == 1
                    M = 1;
                end

                xr_thalamic_std(iter,n_thal,t) = update_xr(M, current_xr_thalamic_std, current_xi_thalamic_std, tau_re_thalamic, tau_ir_thalamic);
                xe_thalamic_std(iter,n_thal,t) = update_xe(M, current_xr_thalamic_std, current_xe_thalamic_std, tau_re_thalamic, tau_ei_thalamic);
                xi_thalamic_std(iter,n_thal,t) = update_xi(current_xe_thalamic_std, current_xi_thalamic_std, tau_ei_thalamic, tau_ir_thalamic);
            
        end
    end

    

    for i=1:n_thalamic
        epsc_thalamic_std(iter,i,:) = reshape(get_g_t_vector(thalamic_poisson_spikes_std(iter,i,:), length(tspan)) .* reshape(xe_thalamic_std(iter,i,:), 1, length(tspan)),  1,1,length(tspan));
    end


    % ---- thalamic epsc ---- dev ---
    for i=1:n_thalamic
         thalamic_poisson_spikes_dev(iter,i, :) = reshape(poisson_generator(lamda_dev_protochol, dt), 1, 1, length(tspan));
    end
    
    % depressing synapse thalamus -> l4
    for n_thal=1:n_thalamic
        for t=2:length(tspan)
                current_xr_thalamic_dev = xr_thalamic_dev(iter, n_thal, t-1);
                current_xe_thalamic_dev = xe_thalamic_dev(iter, n_thal, t-1);
                current_xi_thalamic_dev = xi_thalamic_dev(iter, n_thal, t-1);
                
                M = 0;
                if thalamic_poisson_spikes_dev(iter,n_thal,t) == 1
                    M = 1;
                end

                xr_thalamic_dev(iter,n_thal,t) = update_xr(M, current_xr_thalamic_dev, current_xi_thalamic_dev, tau_re_thalamic, tau_ir_thalamic);
                xe_thalamic_dev(iter,n_thal,t) = update_xe(M, current_xr_thalamic_dev, current_xe_thalamic_dev, tau_re_thalamic, tau_ei_thalamic);
                xi_thalamic_dev(iter,n_thal,t) = update_xi(current_xe_thalamic_dev, current_xi_thalamic_dev, tau_ei_thalamic, tau_ir_thalamic);
            
        end
    end

    

    for i=1:n_thalamic
        epsc_thalamic_dev(iter,i,:) = reshape(get_g_t_vector(thalamic_poisson_spikes_dev(iter,i,:), length(tspan)) .* reshape(xe_thalamic_dev(iter,i,:), 1, length(tspan)),  1,1,length(tspan));
    end


    % ---- thalamic epsc ---- common ---
    for i=1:n_thalamic
         thalamic_poisson_spikes_common(iter,i, :) = reshape(poisson_generator(lamda_common_protochol, dt), 1, 1, length(tspan));
    end
    
    % depressing synapse thalamus -> l4
    for n_thal=1:n_thalamic
        for t=2:length(tspan)
                current_xr_thalamic_common = xr_thalamic_common(iter, n_thal, t-1);
                current_xe_thalamic_common = xe_thalamic_common(iter, n_thal, t-1);
                current_xi_thalamic_common = xi_thalamic_common(iter, n_thal, t-1);
                
                M = 0;
                if thalamic_poisson_spikes_common(iter,n_thal,t) == 1
                    M = 1;
                end

                xr_thalamic_common(iter,n_thal,t) = update_xr(M, current_xr_thalamic_common, current_xi_thalamic_common, tau_re_thalamic, tau_ir_thalamic);
                xe_thalamic_common(iter,n_thal,t) = update_xe(M, current_xr_thalamic_common, current_xe_thalamic_common, tau_re_thalamic, tau_ei_thalamic);
                xi_thalamic_common(iter,n_thal,t) = update_xi(current_xe_thalamic_common, current_xi_thalamic_common, tau_ei_thalamic, tau_ir_thalamic);
            
        end
    end

    

    for i=1:n_thalamic
        epsc_thalamic_common(iter,i,:) = reshape(get_g_t_vector(thalamic_poisson_spikes_common(iter,i,:), length(tspan)) .* reshape(xe_thalamic_common(iter,i,:), 1, length(tspan)),  1,1,length(tspan));
    end

    % simulation
    for i=6:length(tspan)

	fprintf("i = %d\n", i);
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
                    g_t = get_g_t(spike_train_exc, dt, i-5, tspan);
					epsc_ex_neuron_back_c2 = epsc_ex_neuron_back_c2 + g_t*xe(iter, c-2,j,i-5); 
                 end
			end

			epsc_ex_neuron_back_c1 = 0;
			if c-1 >= 1
				for j=1:n_excitatory
					spike_train_exc = spikes(iter,c-1,j,:);
					g_t = get_g_t(spike_train_exc, dt, i-5, tspan);
					epsc_ex_neuron_back_c1 = epsc_ex_neuron_back_c1 + g_t*xe(iter,c-1,j,i-5); 
				end
			end

			epsc_ex_neuron_front_c1 = 0;
			if c+1 <= n_columns
				for j=1:n_excitatory
				    spike_train_exc = spikes(iter,c+1,j,:);
					g_t = get_g_t(spike_train_exc, dt, i-5, tspan);
					epsc_ex_neuron_front_c1 = epsc_ex_neuron_front_c1 + g_t*xe(iter,c+1,j,i-5);
				end
			end	


			epsc_ex_neuron_front_c2 = 0;
			if c+2 <= n_columns
				for j=1:n_excitatory
					spike_train_exc = spikes(iter,c+2,j,:);
					g_t = get_g_t(spike_train_exc, dt, i-5, tspan);
					epsc_ex_neuron_front_c2 = epsc_ex_neuron_front_c2 + g_t*xe(iter,c+2,j,i-5);
				end
			end	


			epsc_ex_own_column = 0;
			for j=1:n_excitatory
				if j == n
					continue;
				end
				spike_train_exc = spikes(iter,c,j,:);
				g_t = get_g_t(spike_train_exc, dt, i-5, tspan);
				epsc_ex_own_column	= epsc_ex_own_column + g_t*xe(iter,c,j,i-5); 
			end

			epsc_inh_own_column = 0;
			for j=n_excitatory+1:n_total_neurons
				if j == n
					continue;
                end
                spike_train_inh = spikes(iter,c,j,:);
				g_t = get_g_t(spike_train_inh, dt, i-5, tspan);
				epsc_inh_own_column = epsc_inh_own_column + g_t*xe(iter,c,j,i-5);
            end
		
          % epsc from thalamic neurons
          epsc_from_thalamic = 0;
          for ttt=1:num_of_input_giving_thalamic
             thalamic_neuron_num = mapping_matrix_thalamic_to_a1(n, ttt);
             epsc_from_thalamic = epsc_from_thalamic + epsc_thalamic(iter,thalamic_neuron_num, i-1);
          end
          % seperate weights thalamic to exc and inh
          if n <= n_excitatory
            epsc_from_thalamic = epsc_from_thalamic*weight_thalamic_to_exc_l4;
          else
             epsc_from_thalamic = epsc_from_thalamic*weight_thalamic_to_inh_l4;
          end

          thalamic_epsc_tensor(iter,c,n,i-5) = epsc_from_thalamic;
          
                    
		  if n <= n_excitatory
			total_epsc = epsc_ex_neuron_back_c2 * J_ee_2 + ...
						epsc_ex_neuron_back_c1 * J_ee_1 + ...
						epsc_ex_neuron_front_c1 * J_ee_1 + ...
						epsc_ex_neuron_front_c2 * J_ee_2 + ...
						epsc_ex_own_column * J_ee_0 + ...
						epsc_inh_own_column * J_ei;
            recurrence_exc_self_column_epsc_tensor(iter,c,n,i-1) = epsc_ex_own_column * J_ee_0;
            recurrence_exc_neighbour_column_epsc_tensor(iter,c,n,i-1) = epsc_ex_neuron_back_c2 * J_ee_2 + ...
						                                                     epsc_ex_neuron_back_c1 * J_ee_1 + ...
						                                                    epsc_ex_neuron_front_c1 * J_ee_1 + ...
						                                                    epsc_ex_neuron_front_c2 * J_ee_2 ;
            
            recurrence_inh_self_column_epsc_tensor(iter,c,n,i-1) = epsc_inh_own_column * J_ei;
            recurrence_inh_neighbour_column_epsc_tensor(iter,c,n,i-1) = 0; 

		  else

			total_epsc = epsc_ex_neuron_back_c2 * J_ie_2 + ...
						epsc_ex_neuron_back_c1 * J_ie_1 + ...
						epsc_ex_neuron_front_c1 * J_ie_1 + ...
						epsc_ex_neuron_front_c2 * J_ie_2 + ...
						epsc_ex_own_column * J_ie_0 + ...
						epsc_inh_own_column * J_ii;
           
            recurrence_exc_self_column_epsc_tensor(iter,c,n,i-1) = epsc_ex_own_column * J_ie_0;
           recurrence_exc_neighbour_column_epsc_tensor(iter,c,n,i-1) = epsc_ex_neuron_back_c2 * J_ie_2 + ...
						epsc_ex_neuron_back_c1 * J_ie_1 + ...
						epsc_ex_neuron_front_c1 * J_ie_1 + ...
						epsc_ex_neuron_front_c2 * J_ie_2;
            recurrence_inh_self_column_epsc_tensor(iter,c,n,i-1) = epsc_inh_own_column * J_ii;
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

	
%	break % for testing only one iteration
end

end