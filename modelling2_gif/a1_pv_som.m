clear all;
close all;
tic ;
n_iters = 1;

% basic variables;
n_columns = 1;
n_excitatory = 20;
n_pv = 3; n_som = 2;    
n_inhibitory = n_pv + n_som;
n_total_neurons = n_excitatory + n_inhibitory;
n_thalamic = 9; num_of_input_giving_thalamic = 4;

% time step
n_tokens = 10;
pre_stimulus_time = 100; post_stimulus_time = 0;
single_stimulus_duration = 50; gap_duration = 300;

physical_time_in_ms = 1; %dt time step
dt = 1;  % 0.2 dt = 20 ms, so 0 .01 = 1 ms
t_simulate = pre_stimulus_time + post_stimulus_time + n_tokens*(single_stimulus_duration + gap_duration);
tspan = 0:dt:t_simulate;

% making bins of 100ms = 20*dt and calculating spike rate
spike_rate_dt = 1*dt;
spike_rate_length = (length(tspan)-1)/(spike_rate_dt/dt);


% connection strength
% weight_reducing_l4 = 0.25; % for now all weights reduced by factor of 0.2
% increase_inhibitory_factor = 75;
% weight_exc_factor = 15;
% exc_to_exc_factor = 6;
% inh_to_exc_factor = 1;
% inh_to_inh_factor = 1;
% 
% J_ee_0 = 6*weight_reducing_l4*weight_exc_factor*exc_to_exc_factor;
% J_ie_0 = 0.5*weight_reducing_l4*weight_exc_factor;
% J_ei = -4*weight_reducing_l4*increase_inhibitory_factor*inh_to_exc_factor;
% J_ii = -0.5*weight_reducing_l4*increase_inhibitory_factor*inh_to_inh_factor;
% J_ee_1 = 0.045*weight_reducing_l4*weight_exc_factor;
% J_ie_1 = 0.0035*weight_reducing_l4*weight_exc_factor;
% J_ee_2 = 0.015*weight_reducing_l4*weight_exc_factor;
% J_ie_2 = 0.0015*weight_reducing_l4*weight_exc_factor;

% within column
J_ee_0 = 135;
J_pv_e_0 = 1.8750;
J_som_e_0 = 1.8750*3;

J_e_pv = -100;
J_pv_pv = -9.3750;
J_som_pv = 0;

J_e_som = -50;
J_som_som = 0;
J_pv_som = -9.3750;

% other column
J_ee_1 = 0.1687;
J_pv_e_1 = 0.0131;
J_som_e_1 = 0.0131;

J_e_som_1 = -10;

J_ee_2 = 0.0562;
J_pv_e_2 = 0.0056;
J_som_e_2 = 0.0056;

J_e_som_2 = -2;

% synaptic weight matrix - exc to exc - row: presyn, col: postsyn
exc_to_exc_weight_matrix = zeros(n_iters, n_columns, length(tspan),n_excitatory, n_excitatory);
minimum_weight_exc_to_exc = 85;
maximum_weight_exc_to_exc = 150;

% plasticity parameters
Amp_strength = 0.015; Amp_weak = 0.021;
tau_strength = 13; tau_weak = 20;

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
xr = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
xe = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
xi = zeros(n_iters, n_columns, n_total_neurons, length(tspan));

xr_thalamic = zeros(n_iters, n_thalamic, length(tspan));
xe_thalamic = zeros(n_iters, n_thalamic, length(tspan));
xi_thalamic = zeros(n_iters, n_thalamic, length(tspan));

% analysis of weights
num_of_LTPs = zeros(n_iters, n_columns, length(tspan));
num_of_LTDs = zeros(n_iters, n_columns, length(tspan));

% mapping from input thalamic neurons to a1 column neurons
all_combinations = nchoosek(1:n_thalamic, num_of_input_giving_thalamic);
mapping_matrix_thalamic_to_a1 = all_combinations(1:n_total_neurons,:);

%% generate inhomo poisson spikes for thalamic neurons
thalamic_poisson_spikes = zeros(n_iters, n_thalamic, length(tspan));

silence_period_pre_and_post = 100;
pre_protochol_silence = zeros(1,pre_stimulus_time);
post_protochol_silence = zeros(1, post_stimulus_time);

lamda_s = 300;
lamda_i = 0;

lamda = zeros(1, n_tokens*(single_stimulus_duration + gap_duration));
for tok=1:n_tokens
    ind = (tok-1)*(single_stimulus_duration+gap_duration) + 1;

    for i=ind:ind+single_stimulus_duration-1
        lamda(1,i) = lamda_s;
    end

    for i=ind+single_stimulus_duration:ind+single_stimulus_duration+gap_duration-1
        lamda(1,i) = lamda_i;
    end

end

protochol = [pre_protochol_silence, lamda, post_protochol_silence, 0]; % 0 is to make n00 to n01 due to tspan length
% calculating epsc of each  thalamic neuron
epsc_thalamic = zeros(n_iters,n_thalamic, length(tspan));

weight_thalamic_to_exc_l4 = 550;
weight_thalamic_to_pv_l4 = 750;
weight_thalamic_to_som_l4 = 750;

%% time constant for synaptic resources
tau_re = 0.6; tau_ir = 700; tau_ei = 15;
tau_re_thalamic = 0.6; tau_ir_thalamic = 300; tau_ei_thalamic = 35;


% initialize
v0 = -70;
xr(:, :, :, 1:5) = 1;
xr_thalamic(:,:,1) = 1;
voltages(:, :, :, 1:5) = v0; %
i1_tensor(:, :, :, 1:5) = 0.01;
i2_tensor(:, :, :, 1:5) = 0.001;
theta_tensor(:, :, :, 1:5) = -50.0;

J_ee_0_initial = 100;
exc_to_exc_weight_matrix(:, :, 1:5, :,:) = J_ee_0_initial;

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

    % thalamic
    for i=1:n_thalamic
        thalamic_poisson_spikes(iter,i, :) = reshape(poisson_generator(protochol, dt), 1, 1, length(tspan));
    end

    % depressing synapse thalamus -> l4
    for n_thal=1:n_thalamic
        for t=2:length(tspan)
            current_xr_thalamic = xr_thalamic(iter, n_thal, t-1);
            current_xe_thalamic = xe_thalamic(iter, n_thal, t-1);
            current_xi_thalamic = xi_thalamic(iter, n_thal, t-1);

            M = 0;
            if thalamic_poisson_spikes(iter,n_thal,t) == 1
                M = 1;
            end

            xr_thalamic(iter,n_thal,t) = update_xr(M, current_xr_thalamic, current_xi_thalamic, tau_re_thalamic, tau_ir_thalamic);
            xe_thalamic(iter,n_thal,t) = update_xe(M, current_xr_thalamic, current_xe_thalamic, tau_re_thalamic, tau_ei_thalamic);
            xi_thalamic(iter,n_thal,t) = update_xi(current_xe_thalamic, current_xi_thalamic, tau_ei_thalamic, tau_ir_thalamic);

        end
    end



    for i=1:n_thalamic
        epsc_thalamic(iter,i,:) = reshape(get_g_t_vector(thalamic_poisson_spikes(iter,i,:), length(tspan)) .* reshape(xe_thalamic(iter,i,:), 1, length(tspan)),  1,1,length(tspan));
    end

    % simulation
    for i=6:length(tspan)

    	fprintf("i = %d\n", i);
    	for c=1:n_columns

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
    			if n > n_excitatory % if n is inhibitory neuron
                    for j=1:n_excitatory
                        if j == n
                            continue;
                        end
                        spike_train_exc = spikes(iter,c,j,:);
                        g_t = get_g_t(spike_train_exc, dt, i-5, tspan);
                        x_e_presyn_neuron = xe(iter,c,j,i-5);
                        epsc_ex_own_column	= epsc_ex_own_column + g_t*x_e_presyn_neuron;
                    end
                else
                    for j=1:n_excitatory
                        if j == n
                            continue;
                        end
                        spike_train_exc = spikes(iter,c,j,:);
                        g_t = get_g_t(spike_train_exc, dt, i-5, tspan);
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
    				g_t = get_g_t(spike_train_inh, dt, i-5, tspan);
                    epsc_pv_own_column = epsc_pv_own_column + g_t*xe(iter,c,j,i-5);
                end


    			epsc_som_own_column = 0;
    			for j=n_excitatory+n_pv+1:n_total_neurons
    				if j == n
    					continue;
                    end
                    spike_train_inh = spikes(iter,c,j,:);
    				g_t = get_g_t(spike_train_inh, dt, i-5, tspan);
                    epsc_som_own_column = epsc_som_own_column + g_t*xe(iter,c,j,i-5);
                end

                epsc_som_back_c2 = 0;
                if c-2 >= 1
    				for j=n_excitatory+n_pv+1:n_total_neurons
    					spike_train_exc = spikes(iter,c-2,j,:);
                        g_t = get_g_t(spike_train_exc, dt, i-5, tspan);
                        epsc_som_back_c2 = epsc_som_back_c2 + g_t*xe(iter, c-2,j,i-5);
                    end
    			end

                epsc_som_back_c1 = 0;
                if c-1 >= 1
    				for j=n_excitatory+n_pv+1:n_total_neurons
    					spike_train_exc = spikes(iter,c-1,j,:);
                        g_t = get_g_t(spike_train_exc, dt, i-5, tspan);
                        epsc_som_back_c1 = epsc_som_back_c1 + g_t*xe(iter, c-1,j,i-5);
                    end
    			end


                epsc_som_front_c1 = 0;
                if c+1 <= n_columns
    				for j=n_excitatory+n_pv+1:n_total_neurons
    				    spike_train_exc = spikes(iter,c+1,j,:);
    					g_t = get_g_t(spike_train_exc, dt, i-5, tspan);
                        epsc_som_front_c1 = epsc_som_front_c1 + g_t*xe(iter,c+1,j,i-5);
    				end
                end

                epsc_som_front_c2 = 0;
                if c+2 <= n_columns
    				for j=n_excitatory+n_pv+1:n_total_neurons
    					spike_train_exc = spikes(iter,c+2,j,:);
    					g_t = get_g_t(spike_train_exc, dt, i-5, tspan);
                        epsc_som_front_c2 = epsc_som_front_c2 + g_t*xe(iter,c+2,j,i-5);
    				end
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
                elseif n > n_excitatory && n <= n_excitatory + n_pv
                    epsc_from_thalamic = epsc_from_thalamic*weight_thalamic_to_pv_l4;
                else
                    epsc_from_thalamic = epsc_from_thalamic*weight_thalamic_to_som_l4;
                end

                thalamic_epsc_tensor(iter,c,n,i-5) = epsc_from_thalamic;

              
      		  if n <= n_excitatory % excitatory neuron
      			total_epsc = epsc_ex_neuron_back_c2 * J_ee_2 + ...
                    epsc_ex_neuron_back_c1 * J_ee_1 + ...
                    epsc_ex_neuron_front_c1 * J_ee_1 + ...
                    epsc_ex_neuron_front_c2 * J_ee_2 + ...
                    epsc_ex_own_column  + ... % weight already included in weight matrix
                    epsc_som_back_c2 * J_e_som_2 + ...
                    epsc_som_back_c1 * J_e_som_1 + ...
                    epsc_som_front_c1 * J_e_som_1 + ...
                    epsc_som_front_c2 * J_e_som_2 + ...
                    epsc_pv_own_column * J_e_pv + ...
                    epsc_som_own_column * J_e_som;
                recurrence_exc_self_column_epsc_tensor(iter,c,n,i-1) = epsc_ex_own_column;
                recurrence_exc_neighbour_column_epsc_tensor(iter,c,n,i-1) = epsc_ex_neuron_back_c2 * J_ee_2 + ...
                                                                            epsc_ex_neuron_back_c1 * J_ee_1 + ...
                                                                            epsc_ex_neuron_front_c1 * J_ee_1 + ...
                                                                            epsc_ex_neuron_front_c2 * J_ee_2 ;

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

            % ------updating weights using STDP
        
                % if change has occured in LTP or LTD for a synapse, 
                % Don't undo it in the next iteration
            either_LTP_or_LTD_occured = zeros(n_excitatory, n_excitatory);
                
            for N=1:n_excitatory
                
                % N -> postyn : LTD
                for postsyn_neuron=1:n_excitatory
                    if spikes(iter,c,N,i) == 1 % if there is a spike
                        found_spike_in_window_LTD = 0;

                        for postsyn_spike_time=i-1:-1:i-20
                            if postsyn_spike_time >= 1 && spikes(iter,c,postsyn_neuron,postsyn_spike_time) == 1 && spikes(iter,c,postsyn_neuron,i) == 0 
                                exc_to_exc_weight_matrix(iter,c,i,N,postsyn_neuron) = exc_to_exc_weight_matrix(iter,c,i-1,N,postsyn_neuron)*(1 - Amp_weak*exp(-abs(i-postsyn_spike_time)/tau_weak));
                                % clipping weights
%                                 if exc_to_exc_weight_matrix(iter,c,i,N,postsyn_neuron) < minimum_weight_exc_to_exc
%                                     exc_to_exc_weight_matrix(iter,c,i,N,postsyn_neuron) = minimum_weight_exc_to_exc;
%                                 end
%                                 if exc_to_exc_weight_matrix(iter,c,i,N,postsyn_neuron) > maximum_weight_exc_to_exc
%                                     exc_to_exc_weight_matrix(iter,c,i,N,postsyn_neuron) = maximum_weight_exc_to_exc;
%                                 end
                                num_of_LTDs(iter,c,i) = num_of_LTDs(iter,c,i) + 1;
                                
%                                 if N == 5 && postsyn_neuron == 7
%                                     fprintf("\n LTD - !!! - old %f, new %f \n ",exc_to_exc_weight_matrix(iter,c,i-1,N,postsyn_neuron),exc_to_exc_weight_matrix(iter,c,i,N,postsyn_neuron))
%                                     fprintf("\n due to same time ??? %d %d \n", spikes(iter,c,N,i), spikes(iter,c,postsyn_neuron,i))
%                                 end
                                either_LTP_or_LTD_occured(N,postsyn_neuron) = 1;
                                found_spike_in_window_LTD = 1;
                                break;
                            end
                        end

                        if found_spike_in_window_LTD == 0 && either_LTP_or_LTD_occured(N,postsyn_neuron) == 0
                            exc_to_exc_weight_matrix(iter,c,i,N,postsyn_neuron) = exc_to_exc_weight_matrix(iter,c,i-1,N,postsyn_neuron);
%                             if N == 5 && postsyn_neuron == 7
%                                     fprintf("\n LTD - @@@ - old %f, new %f \n ",exc_to_exc_weight_matrix(iter,c,i-1,N,postsyn_neuron),exc_to_exc_weight_matrix(iter,c,i,N,postsyn_neuron))
%                             end
                        end
                    else % if there is no spike

                        if either_LTP_or_LTD_occured(N,postsyn_neuron) == 0
                               exc_to_exc_weight_matrix(iter,c,i,N,postsyn_neuron) = exc_to_exc_weight_matrix(iter,c,i-1,N,postsyn_neuron);
%                                 if N == 5 && postsyn_neuron == 7
%                                     fprintf("\n LTD - $$$ - old %f, new %f \n ",exc_to_exc_weight_matrix(iter,c,i-1,N,postsyn_neuron),exc_to_exc_weight_matrix(iter,c,i,N,postsyn_neuron))
%                                 end
                        end
                        
                    end
                end


                % presyn -> N : LTP
                for presyn_neuron=1:n_excitatory
                    if spikes(iter,c,N,i) == 1 % if there is a spike
                        found_spike_in_window_LTP = 0;

                        for presyn_spike_time=i-1:-1:i-20
                            if presyn_spike_time >= 1 && spikes(iter,c,presyn_neuron,presyn_spike_time) == 1 && spikes(iter,c,presyn_neuron,i) == 0
                                    exc_to_exc_weight_matrix(iter,c,i,presyn_neuron,N) = exc_to_exc_weight_matrix(iter,c,i-1,presyn_neuron,N)*(1 + Amp_strength*exp(-abs(i-presyn_spike_time)/tau_strength));
                                    % clipping weights 
%                                     if exc_to_exc_weight_matrix(iter,c,i,presyn_neuron,N) < minimum_weight_exc_to_exc
%                                         exc_to_exc_weight_matrix(iter,c,i,presyn_neuron,N) = minimum_weight_exc_to_exc;
%                                     end
%                                     if exc_to_exc_weight_matrix(iter,c,i,presyn_neuron,N) > maximum_weight_exc_to_exc
%                                         exc_to_exc_weight_matrix(iter,c,i,presyn_neuron,N) = maximum_weight_exc_to_exc;
%                                     end

                                    num_of_LTPs(iter,c,i) = num_of_LTPs(iter,c,i) + 1; 
%                                     if presyn_neuron == 5 && N == 7
%                                         fprintf("\n LTP - !!! - old %f, new %f \n ",exc_to_exc_weight_matrix(iter,c,i-1,presyn_neuron,N),exc_to_exc_weight_matrix(iter,c,i,presyn_neuron,N))
%                                     end
                                
                                    
                                
                                either_LTP_or_LTD_occured(presyn_neuron,N) = 1;
                                found_spike_in_window_LTP = 1;
                                break
                            end
                        end

                        if found_spike_in_window_LTP == 0 && either_LTP_or_LTD_occured(presyn_neuron,N) == 0
                            exc_to_exc_weight_matrix(iter,c,i,presyn_neuron,N) = exc_to_exc_weight_matrix(iter,c,i-1,presyn_neuron,N);
%                             if presyn_neuron == 5 && N == 7
%                                     fprintf("\n LTP - @@@ - old %f, new %f \n ",exc_to_exc_weight_matrix(iter,c,i-1,presyn_neuron,N),exc_to_exc_weight_matrix(iter,c,i,presyn_neuron,N))
%                             end
                        end
                    else % if there is no spike
                        if either_LTP_or_LTD_occured(presyn_neuron,N) == 0
                            exc_to_exc_weight_matrix(iter,c,i,presyn_neuron,N) = exc_to_exc_weight_matrix(iter,c,i-1,presyn_neuron,N);
%                             if presyn_neuron == 5 && N == 7
%                                     fprintf("\n LTP - $$$ - old %f, new %f \n ",exc_to_exc_weight_matrix(iter,c,i-1,presyn_neuron,N),exc_to_exc_weight_matrix(iter,c,i,presyn_neuron,N))
%                             end
                        end
                     end
                end

                

            end % end of for all excitatory neurons



            %     fprintf("xr %f, xe %f, xi %f\n", xr(c,n,i),xe(c,n,i), xi(c,n,i));
            % pause(0.4);

        end


        %	break % for testing only one iteration
    end

end
toc;
save('batch_1.mat');
return
%% -  analysis single column
% fill the spike rates tensor
for i=1:n_iters
    for n=1:n_total_neurons
        spikes1 = spikes(i, 1, n, :);
        spikes1_reshaped = reshape(spikes1, 1,length(tspan));
        spike_rate1 = spikes_to_spike_rate_neat(spikes1_reshaped, physical_time_in_ms, dt, spike_rate_dt);
        spikes_rate1 = reshape(spike_rate1, 1,1,1,spike_rate_length);
        spike_rates(i,1,n,:) = spikes_rate1;
    end
end

% -- plot total input epsc to l4 (leave I_background)
total_input_epsc = thalamic_epsc_tensor ...
    + recurrence_exc_self_column_epsc_tensor + recurrence_inh_self_column_epsc_tensor ...
    + recurrence_exc_neighbour_column_epsc_tensor + recurrence_inh_neighbour_column_epsc_tensor ;

% if clip is on
% for iter=1:n_iters
%     for col=1:n_columns
%         for n=1:n_total_neurons
%             for t=1:t_simulate
%                 if total_input_epsc(iter,col,n,t) < 0
%                     total_input_epsc(iter,col,n,t) = 0;
%                 end
%             end
%         end
%     end
% end

[mean_input_epsc_exc_for_iters, mean_input_epsc_exc_for_neurons] = get_mean(total_input_epsc(:,:,1:n_excitatory,:), n_iters, n_excitatory, length(tspan)-1, 1);
[mean_input_epsc_inh_for_iters, mean_input_epsc_inh_for_neurons] = get_mean(total_input_epsc(:,:,n_excitatory+1:n_total_neurons,:), n_iters, n_inhibitory, length(tspan)-1, 1);
[mean_input_epsc_all_for_iters, mean_input_epsc_all_for_neurons] = get_mean(total_input_epsc, n_iters, n_total_neurons, length(tspan)-1, 1);
figure
subplot(1,2,1)
hold on

plot(mean_input_epsc_exc_for_neurons)
plot(mean_input_epsc_inh_for_neurons)
plot(mean_input_epsc_all_for_neurons)
legend('total input to exc l4', 'total input to inh l4', 'total input to l4 all')
title('l4 - total input epsc')
hold off

subplot(1,2,2)
plot(mean_input_epsc_all_for_neurons)
title('l4 total input epsc')

grid

% -- plot psth of l4 exc, inh, all
figurea1
subplot(1,2,1)
hold on
[mean_spike_rate_exc_for_iters, mean_spike_rate_exc_for_neurons] = get_mean(spike_rates(:,:,1:n_excitatory,:), n_iters, n_excitatory, spike_rate_length,1);
[mean_spike_rate_inh_for_iters, mean_spike_rate_inh_for_neurons] = get_mean(spike_rates(:,:,n_excitatory+1:n_total_neurons,:), n_iters, n_inhibitory, spike_rate_length,1);
[mean_spike_rate_for_iters, mean_spike_rate_for_neurons] = get_mean(spike_rates, n_iters, n_total_neurons, spike_rate_length,1);

n_bins = spike_rate_dt/dt;

plot(mean_spike_rate_exc_for_neurons*(n_bins*physical_time_in_ms*0.001))
plot(mean_spike_rate_inh_for_neurons*(n_bins*physical_time_in_ms*0.001))
plot(mean_spike_rate_for_neurons*(n_bins*physical_time_in_ms*0.001));
legend('psth l4 exc', 'psth l4 inh','psth l4 all')
hold off

subplot(1,2,2)
plot(mean_spike_rate_for_neurons*(n_bins*physical_time_in_ms*0.001))
title('psth l4 all')
grid

recurrence_exc_epsc = recurrence_exc_self_column_epsc_tensor + recurrence_exc_neighbour_column_epsc_tensor;
recurrence_inh_epsc = recurrence_inh_self_column_epsc_tensor + recurrence_inh_neighbour_column_epsc_tensor;
recurrence_epsc = recurrence_exc_epsc + recurrence_inh_epsc;

figure
hold on
[mean_epsc_exc_recurrence_for_iters, mean_epsc_exc_reccurence_for_neurons] = get_mean(recurrence_exc_epsc, n_iters, n_total_neurons, length(tspan)-1, 1);
[mean_epsc_inh_recurrence_for_iters, mean_epsc_inh_reccurence_for_neurons] = get_mean(recurrence_inh_epsc, n_iters, n_total_neurons, length(tspan)-1, 1);

plot(mean_epsc_exc_reccurence_for_neurons)
plot(mean_epsc_inh_reccurence_for_neurons)
legend('exc epsc to all l4', 'inh epsc to all l4')
title('epscs-exc and inh')
hold off
grid
figure
subplot(1,2,1)
hold on
[mean_recurrence_epsc_exc_for_iters, mean_recurrence_epsc_exc_for_neurons] = get_mean(recurrence_epsc(:,:,1:n_excitatory,:), n_iters, n_excitatory, length(tspan)-1,1);
[mean_recurrence_epsc_inh_for_iters, mean_recurrence_epsc_inh_for_neurons] = get_mean(recurrence_epsc(:,:,n_excitatory+1:n_total_neurons,:), n_iters, n_inhibitory, length(tspan)-1,1);
[mean_recurrence_epsc_all_for_iters, mean_recurrence_epsc_all_for_neurons] = get_mean(recurrence_epsc, n_iters, n_total_neurons, length(tspan)-1, 1);

plot(mean_recurrence_epsc_exc_for_neurons)
plot(mean_recurrence_epsc_inh_for_neurons)
plot(mean_recurrence_epsc_all_for_neurons)
legend('recurrence epsc input to exc in l4', 'recurrence epsc to inh in l4', 'recurrence epsc to all in l4')
title('recurrrence epscs')
hold off

subplot(1,2,2)
plot(mean_recurrence_epsc_all_for_neurons)
title('mean recurrence epsc to all l4 neurons')


grid

% --- thalamic epsc to l4 --
[mean_thalamic_epsc_for_iters, mean_thalamic_epsc_for_neurons] = get_mean(thalamic_epsc_tensor, n_iters, n_total_neurons, length(tspan)-1, 1);
figure
plot(mean_thalamic_epsc_for_neurons)
title('thalamic epsc to l4 all')
grid

figure
[mean_spike_rate_for_iters, mean_spike_rate_for_neurons] = get_mean(spike_rates, n_iters, n_total_neurons, spike_rate_length,1);
n_bins = spike_rate_dt/dt;
hold on
mean_input_epsc_extended = zeros(1, length(tspan));
mean_input_epsc_extended(1, 2:length(tspan)) = mean_input_epsc_all_for_neurons;
mean_input_epsc_binned = spikes_to_spike_rate_neat(mean_input_epsc_extended, 1, dt, spike_rate_dt);
plot(mean_input_epsc_binned)
plot(mean_spike_rate_for_neurons/0.001)

plot(protochol);
hold off
title('total input and psth')
legend('total input epsc', 'psth l4', 'protochol')
grid


% 5 bins - psth and epsc
spike_rate_dt_5 = 5*dt;
spike_rate_length_5 = (length(tspan)-1)/(spike_rate_dt_5/dt);
spike_rates_5 = zeros(n_iters, n_columns, n_total_neurons, spike_rate_length_5);

for i=1:n_iters
    for n=1:n_total_neurons
        spikes1 = spikes(i, 1, n, :);
        spikes1_reshaped = reshape(spikes1, 1,length(tspan));
        spike_rate1 = spikes_to_spike_rate_neat(spikes1_reshaped, physical_time_in_ms, dt, spike_rate_dt_5);
        spikes_rate1 = reshape(spike_rate1, 1,1,1,spike_rate_length_5);
        spike_rates_5(i,1,n,:) = spikes_rate1;
    end
end

figure
[mean_spike_rate_for_iters, mean_spike_rate_for_neurons] = get_mean(spike_rates_5, n_iters, n_total_neurons, spike_rate_length_5,1);
n_bins = spike_rate_dt_5/dt;
hold on
mean_input_epsc_extended = zeros(1, length(tspan));
mean_input_epsc_extended(1, 2:length(tspan)) = mean_input_epsc_all_for_neurons;
mean_input_epsc_binned = spikes_to_spike_rate_neat(mean_input_epsc_extended, 1, dt, spike_rate_dt_5);
plot(mean_input_epsc_binned)
plot(mean_spike_rate_for_neurons/(n_bins*0.001))

protochol_binned = spikes_to_spike_rate_neat(protochol, 1, dt, spike_rate_dt_5);
plot(protochol_binned);
hold off
title('5 bins total input and psth')
legend('total input epsc', 'psth l4','protochol')
grid

% for a random neuron
random_neuron = 23;
col = 1;
figure
hold on
mean_spike_rate_of_random_neuron = zeros(1, spike_rate_length);
for t=1:spike_rate_length
    mean_spike_rate_of_random_neuron(1,t) = sum(spike_rates(:, col, random_neuron,t))/n_iters;
end
plot(mean_spike_rate_of_random_neuron);

% adjusting from length(tspan)-1 to length(tspan)
epsc_all_iters_for_random_neuron = squeeze(total_input_epsc(:, col, random_neuron, :));
actual_epsc_input_to_random_neuron_size_adjusted = zeros(n_iters, length(tspan));
for iter=1:n_iters
    actual_epsc_input_to_random_neuron_size_adjusted(iter, 2:length(tspan)) = epsc_all_iters_for_random_neuron(iter, :);
end

% taking avg iters
random_epsc_iters_avg = zeros(1, length(tspan));
for t=1:length(tspan)
    random_epsc_iters_avg(1,t) = sum(actual_epsc_input_to_random_neuron_size_adjusted(:,t))/n_iters;
end

epsc_to_random_neuron_binned = spikes_to_spike_rate_neat(random_epsc_iters_avg, physical_time_in_ms, dt, spike_rate_dt);
plot(epsc_to_random_neuron_binned*(n_bins*physical_time_in_ms*0.001))

hold off
title('random neuron epsc and psth')
legend('psth', 'epsc')
grid

% raster of l4
% figure
%     c = 1;
%     spike_reshaped = reshape(spikes(:,c,:,:),  n_iters*n_total_neurons, length(tspan));
%     imagesc(spike_reshaped);
%     title('all iters raster l4')
% grid

for iter=1:n_iters
    figure(iter*100 + 77)
    hold on
    c = 1;
    spike_reshaped = reshape(spikes(iter,c,:,:),  n_total_neurons, length(tspan));
    imagesc(spike_reshaped);
    title('single iter raster l4')
    plot(protochol/10)
    hold off
    grid
end

%% - end of analysis single colum -
% var_tensor, n_iters, n_neurons, time_length,column_index
return
% -------------------------------
% -------- psth -------------
% fill the spike rates tensor - re run
% spike_rate_dt = 5*dt;
% spike_rate_length = (length(tspan)-1)/(spike_rate_dt/dt);
% spike_rates = zeros(n_iters, n_columns, n_total_neurons, spike_rate_length);

for i=1:n_iters
    for n=1:n_total_neurons
        spikes1 = spikes(i, 1, n, :);
        spikes1_reshaped = reshape(spikes1, 1,length(tspan));
        spike_rate1 = spikes_to_spike_rate_neat(spikes1_reshaped, physical_time_in_ms, dt, spike_rate_dt);
        spikes_rate1 = reshape(spike_rate1, 1,1,1,spike_rate_length);
        spike_rates(i,1,n,:) = spikes_rate1;
    end
end


% get a mean of all spikes
spike_rate_l4 = zeros(n_total_neurons, spike_rate_length);
for i=1:spike_rate_length
    for n=1:n_total_neurons
        spike_rate_l4(n, i) = sum(spike_rates(:,1,n,i))/n_iters;
    end
end
% mean psth of all neurons
spike_rate_l4_all = zeros(1, spike_rate_length);
for i=1:spike_rate_length
    spike_rate_l4_all(1, i) = sum(spike_rate_l4(:, i))/n_total_neurons;
end


figure
plot(spike_rate_l4_all);
title('psth of l4  all neurons')
grid

figure
n_bins = spike_rate_dt/dt;
plot(spike_rate_l4_all*(n_bins*physical_time_in_ms*0.001));
title('num of spikes of l4  all neurons')
grid
% -------- psth end -------------

% testing spikes of all 25 neurons
x = squeeze(spikes);
x1 = reshape(x, n_iters*n_total_neurons, length(tspan));
figure
imagesc(x1)
grid

return

%% ---- random neuron voltage plot
figure
x = reshape(thalamic_poisson_spikes, n_iters*n_thalamic, length(tspan));
imagesc(x);
title('thalamic poision spikes')
grid

x2 = reshape(feedforward_epsc_tensor, n_iters*n_total_neurons, length(tspan)-1);
figure
imagesc(x2);
title('feedforward epsc tensor')
grid

x1 = reshape(recurrence_epsc_tensor, n_iters*n_total_neurons, length(tspan)-1);
figure
imagesc(x1);
title('epsc recurrence of 25 neurons')
grid

nnn = 7;
x = voltages(1,1,nnn,:);
t = theta_tensor(1,1,nnn,:);
x = squeeze(x);
t = squeeze(t);
s = spikes(1,1,nnn,:);
s = squeeze(s)*-10;

figure
hold on
plot(x)
plot(t)
plot(s);
title('voltage and threshold of random neuron')
legend('voltage', 'threshold','spikes')
hold off
grid
%% ---- graph end

% --------- raster plot ----
% testing spikes of all 25 neurons
x = squeeze(spikes);
x1 = reshape(x, n_iters*n_total_neurons, length(tspan));
figure
imagesc(x1)
grid


xx = squeeze(recurrence_epsc_tensor);
yy = zeros(1, length(tspan)-1);
for i=1:length(tspan)-1
    yy(1,i) = sum(xx(:, i))/n_total_neurons;
end
figure
plot(yy)
title('recurrence epsc avg')
grid

xx = squeeze(xe);
yy = zeros(1, length(tspan)-1);
for i=1:length(tspan)-1
    yy(1,i) = sum(xx(:, i))/n_total_neurons;
end
figure
plot(yy)
title('xe')
grid

% raster psth
figure
imagesc(spike_rate_l4);
title('raster psth')
grid

% -------- recurence ------------
% get a mean of all spikes
recurrence_avg = zeros(n_total_neurons, length(tspan)-1);
for i=1:length(tspan)-1
    for n=1:n_total_neurons
        recurrence_avg(n, i) = sum(recurrence_epsc_tensor(:,1,n,i))/n_iters;
    end
end

% mean psth of all neurons
recurrence_avg_all = zeros(1, length(tspan)-1);
for i=1:length(tspan)-1
    recurrence_avg_all(1, i) = sum(recurrence_avg(:, i))/n_total_neurons;
end

figure
plot(recurrence_avg_all);
title('recurrence of l4  all neurons')
grid

% -------------------------------

% ------ epsc thalamic-----
epsc_thalamic_avg = zeros(n_thalamic, length(tspan)-1);
for n=1:n_thalamic
    for t=1:length(tspan)-1
        epsc_thalamic_avg(n,t) = sum(epsc_thalamic(:, n, t))/n_iters;
    end
end

epsc_thalamic_avg_all = zeros(1, length(tspan)-1);
for t=1:length(tspan)-1
    epsc_thalamic_avg_all(1,t) = sum(epsc_thalamic_avg(:,t))/n_thalamic;
end

figure
plot(epsc_thalamic_avg_all);
title('epsc thalamic avg all')
grid

% ---- epsc thalamic end ----


% ------ epsc thalamic-----
e1 = squeeze(epsc_with_background_tensor);
e_avg = zeros(1, length(tspan)-1);
for t=1:length(tspan)-1
    e_avg(1, t) = sum(e1(:, t))/n_total_neurons;
end
figure
plot(e_avg)
title('epsc with background')
grid

% raster of epsc with background
e_squezed = squeeze(epsc_with_background_tensor);
e_squeezed_reshaped = reshape(e_squezed, n_iters*n_total_neurons, length(tspan)-1);
figure
imagesc(e_squeezed_reshaped)
title('raster of epsc with background')
grid

% ---- epsc end ----

% ------ epsc excitatory -------
% for 1 column
c=1
epsc_exc_avged_iters = zeros(n_total_neurons, length(tspan)-1);
for t=1:length(tspan)-1
    for n=1:n_total_neurons
        epsc_exc_avged_iters(n,t) = sum(epsc_exc_tensor(:,c,n,t))/n_iters;
    end
end

% avg for all neurons
epsc_exc_aveg_neurons = zeros(1, length(tspan)-1);
for t=1:length(tspan)-1
    epsc_exc_aveg_neurons(1,t) = sum(epsc_exc_avged_iters(:,t))/n_total_neurons;
end

figure
plot(epsc_exc_aveg_neurons);
title('recurrence exc epsc avg l4 neurons')
grid
% --- - epsc inhi ----
epsc_inh_avged_iters = zeros(n_total_neurons, length(tspan)-1);
for t=1:length(tspan)-1
    for n=1:n_total_neurons
        epsc_inh_avged_iters(n,t) = sum(epsc_inh_tensor(:,c,n,t))/n_iters;
    end
end

% avg for all neurons
epsc_inh_aveg_neurons = zeros(1, length(tspan)-1);
for t=1:length(tspan)-1
    epsc_inh_aveg_neurons(1,t) = sum(epsc_inh_avged_iters(:,t))/n_total_neurons;
end

figure
plot(epsc_inh_aveg_neurons);
title('recurrence inh epsc avg l4 neurons')
grid


return



% hold on - epsc and voltage of random neuron
figure
hold on
voltage1 = voltages(1,1,10,:);
v1 = reshape(voltage1, 1,length(tspan));
plot(v1);

epsc_from_thalamic = 0;
n=10;
for ttt=1:num_of_input_giving_thalamic
    thalamic_neuron_num = mapping_matrix_thalamic_to_a1(n, ttt);
    epsc_from_thalamic = epsc_from_thalamic + epsc_thalamic(1,thalamic_neuron_num, :);
end
ee = reshape(epsc_from_thalamic,1,length(tspan));
ef = ee + I_background_tensor(1, :);
plot(ee);
hold off

title('voltage and epsc')
grid




%------- excitatory and inhibitory------
% spike_rate_l4 - 25 x spike_rate_length
psth_exc_avg = zeros(1, spike_rate_length);
for i=1:spike_rate_length
    psth_exc_avg(1,i) = sum(spike_rate_l4(1:n_excitatory,i))/n_excitatory;
end


psth_inh_avg = zeros(1, spike_rate_length);
for i=1:spike_rate_length
    psth_inh_avg(1,i) = sum(spike_rate_l4(n_excitatory+1:n_total_neurons,i))/n_inhibitory;
end

figure
hold on
plot(spike_rate_l4_all);
plot(psth_exc_avg);
plot(psth_inh_avg);
legend('all', 'exc', 'inh')
hold off
grid
% --------------------------------------

% ------ an excitatory input and psth -----
random_exc = 10;
random_neuron_psth = spike_rate_l4(random_exc, :);
epsc_tensor_for_random = zeros(1,spike_rate_length);
epsc_squeeze = squeeze(epsc_tensor);
x = zeros(1, length(tspan));
y = zeros(1, length(tspan)-1);
for i=2:length(tspan)-1
    x(1,i) = sum(epsc_squeeze(:,random_exc, i))/n_iters;
end

for k=2:length(tspan)
    y(k) = sum(I_background_tensor(:, k))/n_iters;
end

x1 = spikes_to_spike_rate_neat(x,physical_time_in_ms, dt, spike_rate_dt);
y1 = spikes_to_spike_rate_neat(y,physical_time_in_ms, dt, spike_rate_dt);
z1 = x1 + y1;
figure
hold on
plot(z1);
plot(random_neuron_psth);
legend('input curent', 'psth of random neuron')

hold off
grid
% ---------------------

% ------ epsc thalamic-----
epsc_thalamic_squeezed = squeeze(epsc_thalamic);
epsc_thalamic_avg = zeros(1, length(tspan));
for t=1:length(tspan)
    epsc_thalamic_avg(1, t) = sum(epsc_thalamic_squeezed(:, t))/n_thalamic;
end
figure
plot(epsc_thalamic_avg)
title('epsc thalamic avg')
grid
% ---- epsc thalamic end ----
return











%------ l4 neurons------
% fill the spikes tensor
for i=1:n_iters
    for n=1:n_total_neurons
        voltage1 = voltages(i, 1, n, :);
        voltage1_reshaped = reshape(voltage1, 1, length(tspan));
        spikes1 = voltage_to_spikes(voltage1_reshaped);
        spikes1 = reshape(spikes1, 1, 1, 1,length(tspan));
        spikes(i,1,n,:) = spikes1;
    end
end



% epsc thalamic
epsc_thalamic_squeezed = squeeze(epsc_thalamic);
epsc_thalamic_avg = zeros(1, length(tspan));
for t=1:length(tspan)
    epsc_thalamic_avg(1, t) = sum(epsc_thalamic_squeezed(:, t))/n_thalamic;
end
figure
plot(epsc_thalamic_avg)
title('epsc thalamic avg')
grid

% epsc weight
epsc_tensor_squeeze = squeeze(epsc_tensor);
epsc_all_avg = zeros(1, length(tspan)-1);
for i=1:length(tspan)-1
    epsc_all_avg(1, i) = sum(epsc_tensor_squeeze(:, i))/n_total_neurons;
end
figure
plot(epsc_all_avg)
title('epsc all neurons avg')
grid
return

% epsc input average
epsc_squeeze = squeeze(epsc_tensor);
epsc_average = zeros(n_total_neurons, length(tspan)-1);
for n=1:n_total_neurons
    for t=2:length(tspan)-1
        epsc_average(n, t) = sum(epsc_squeeze(n,t))/n_iters + I_background_tensor(1,t);
    end
end

% random neuron epsc and avg voltage
voltage_avg = zeros(n_total_neurons, length(tspan)-1);
voltages_squeeze = squeeze(voltages);
for n=1:n_total_neurons
    for t=2:length(tspan)
        voltage_avg(n,t) = sum(voltages_squeeze(:,n,t))/n_iters;
    end
end


figure
hold on
plot(epsc_average(10,:))
title('random neuron epsc')
v_avg = voltage_avg(10,:);
plot(v_avg);
hold off
grid

epsc_avg_l4 = zeros(1,length(tspan)-1);
for t=1:length(tspan)-1
    epsc_avg_l4(1,t) = sum(epsc_average(:,t))/n_total_neurons;
end
figure
plot(epsc_avg_l4);
title('epsc avg l4')
grid

% testing if really thalamic inputs vary
figure
hold on
for n=1:n_total_neurons
    epsc_from_thalamic = 0;
    for ttt=1:num_of_input_giving_thalamic
        thalamic_neuron_num = mapping_matrix_thalamic_to_a1(n, ttt);
        epsc_from_thalamic = epsc_from_thalamic + epsc_thalamic(1,thalamic_neuron_num, :);
    end
    plot(reshape(epsc_from_thalamic, 1, length(tspan)));

end
hold off
grid

figure
hold on
voltage1 = voltages(1,1,10,:);
v1 = reshape(voltage1, 1,length(tspan));
plot(v1);

epsc_from_thalamic = 0;
n=10;
for ttt=1:num_of_input_giving_thalamic
    thalamic_neuron_num = mapping_matrix_thalamic_to_a1(n, ttt);
    epsc_from_thalamic = epsc_from_thalamic + epsc_thalamic(1,thalamic_neuron_num, :);
end
ee = reshape(epsc_from_thalamic,1,length(tspan));
ef = ee + I_background_tensor(1, :);
plot(ee);
hold off
grid


%





x = squeeze(spikes);
x1 = reshape(x, 25, 2001);
figure
imagesc(x1)
grid

% x = reshape(spikes, 25, 2001);
% figure
%     imagesc(x)
%     title('l4')
% grid

return


% fill the spike rates tensor - re run
for i=1:n_iters
    for n=1:n_total_neurons
        spikes1 = spikes(i, 1, n, :);
        spikes1_reshaped = reshape(spikes1, 1,length(tspan));
        spike_rate1 = spikes_to_spike_rate_neat(spikes1_reshaped, physical_time_in_ms, dt, spike_rate_dt);
        spikes_rate1 = reshape(spike_rate1, 1,1,1,spike_rate_length);
        spike_rates(i,1,n,:) = spikes_rate1;
    end
end


% get a mean of all spikes
spike_rate_l4 = zeros(n_total_neurons, spike_rate_length);
for i=1:spike_rate_length
    for n=1:n_total_neurons
        spike_rate_l4(n, i) = sum(spike_rates(:,1,n,i))/n_iters;
    end
end

% check spike rate of some random l4 neuron
figure
test_neuron = 10;
spike_rate_test_l4 = spike_rate_l4(test_neuron, :);
plot(spike_rate_test_l4);
grid


%  raster plot
raster_all = [];
figure
spikes_test = reshape(spikes(:, 1, n, :), n_iters,length(spikes(:, 1, n, :)));
imagesc(spikes_test)
title('raster of l4 all neuron')
grid

% input to test neuron
figure
current = reshape(epsc_tensor(1,1,test_neuron, :), 1,length(epsc_tensor(10,1,test_neuron, :)));
plot(current);
title('epsc to test')
grid

% mean psth of all neurons
spike_rate_l4_all = zeros(1, spike_rate_length);
for i=1:spike_rate_length
    spike_rate_l4_all(1, i) = sum(spike_rate_l4(:, i))/n_total_neurons;
end

figure
plot(spike_rate_l4_all);
title('psth of l4  all neurons')
grid

% ---- thalamic neurons ------
thalamic_spike_rates = zeros(n_iters, n_thalamic, spike_rate_length);
for i=1:n_iters
    for n=1:n_thalamic
        spikes1 = thalamic_poisson_spikes(i,n,:);
        spike_rate1 = spikes_to_spike_rate_neat(spikes1, physical_time_in_ms, dt, spike_rate_dt);
        spike_rate1 = reshape(spike_rate1, 1,1,spike_rate_length);
        thalamic_spike_rates(i,n,:) = spike_rate1;
    end
end

thalamic_spike_rate_avg = zeros(n_thalamic, spike_rate_length);
for i=1:spike_rate_length
    for n=1:n_thalamic
        thalamic_spike_rate_avg(n,i) = sum(thalamic_spike_rates(:,n,i))/n_iters;
    end
end

thalamic_spike_rate_all_avg = zeros(1, spike_rate_length);
for i=1:spike_rate_length
    thalamic_spike_rate_all_avg(1, i) = sum(thalamic_spike_rate_avg(:,i))/n_thalamic;
end

figure
plot(thalamic_spike_rate_all_avg);
title('psth of all thalamic neurons')
grid
% ========== remove later after checking =============

% ff=[];
% for ii=1:25
%     qw=mean_spikes(ii,1:2000);
%     ff=[ff;mean(reshape(qw,5,400))];
% end
% figure(457)
%     plot(mean(ff))
%     title('psth of l4')
% grid

% ff=[];
% for ii=1:9
%     qw=thalamic_mean_spikes(ii,1:2000);
%     ff=[ff;mean(reshape(qw,5,400))];
% end
% figure(754)
%     plot(mean(ff))
%     title('psth of thalmic')
% grid


%{
% see later

% average of iterations
mean_spikes = zeros(n_total_neurons, length(tspan));
for t=1:length(tspan)
    sum_spikes = ;
end

testing_column = 1;
testing_neuron = 15;
thalamic_testing_neuron = 6;

figure(1)
    vvv = reshape(voltages(1, testing_neuron, :),1, length(tspan));
    plot(vvv);
    title('voltage of neuron')
grid

figure(2)
    eee = reshape(epsc_tensor(1, testing_neuron, :),1, length(tspan)-1);
    eee = eee + I_background;
    plot(eee);
    title('total curret')
grid




% population behaviour - mean spike rate of column with time
n_bins = spike_rate_dt/dt;
original_length_spikes = length(tspan);

allneurons_spike_rates = zeros(n_total_neurons, (original_length_spikes-1)/n_bins);
for i=1:n_total_neurons
     spikes1 = voltage_to_spikes(voltages(1, i, :));
     spike_rates1 = spikes_to_spike_rate_neat(spikes1, physical_time_in_ms, dt, spike_rate_dt);
     allneurons_spike_rates(i,:) = spike_rates1;
end
population_psth = zeros(1, (original_length_spikes-1)/n_bins);
for i=1:(original_length_spikes-1)/n_bins
    population_psth(1,i) = sum(allneurons_spike_rates(:,i))/n_total_neurons;
end


% psth of l4 neurons
ff=[];
for ii=1:25
    qw=spikes_2d_matrix_l4(ii,1:2500);
    ff=[ff;mean(reshape(qw,5,500))];
end
figure(457)
    plot(mean(ff))
    title('psth of l4')
grid


% spikes of l4 neurons
spikes_2d_matrix_l4 = zeros(n_total_neurons, length(tspan));
for i=1:n_total_neurons
    spikes_2d_matrix_l4(i,:) = voltage_to_spikes(voltages(1,i,:));
end
figure(770)
    imagesc(spikes_2d_matrix_l4);
grid

figure(9)
    plot(population_psth);
    title('psth of all neurons in column')
grid

% thalamic neurons - mean spike rate
thalamic_neurons_spikes_rates = zeros(n_thalamic, (original_length_spikes-1)/n_bins);
for i=1:n_thalamic
   spikes11 = poisson_generator(lamda, dt, length(tspan));
    spike_rates11 = spikes_to_spike_rate_neat(spikes11, physical_time_in_ms, dt, spike_rate_dt);
   thalamic_neurons_spikes_rates(i,:) = spike_rates11;
end
thalamic_population_psth = zeros(1, (original_length_spikes-1)/n_bins);
for i=1:(original_length_spikes-1)/n_bins
    thalamic_population_psth(1,i) = sum(thalamic_neurons_spikes_rates(:,i))/n_thalamic;
end



figure(3)
    spikes1 = voltage_to_spikes(voltages(testing_column, testing_neuron, :));
    stem(tspan, spikes1);
    title('spikes')
grid

figure(7)
    spikes1 = voltage_to_spikes(voltages(testing_column, testing_neuron, :));
	spike_rates1 = spikes_to_spike_rate_neat(spikes1, physical_time_in_ms, dt, spike_rate_dt);
    plot(spike_rates1);
    title('spike rate')
grid

figure(57)
    % for a single column
   reshaped_spikes =  zeros(n_total_neurons, length(tspan));
    for i=1:n_total_neurons
        reshaped_voltage = reshape(voltages(1,i,:), 1, length(tspan));
        reshaped_spikes(i,:) = voltage_to_spikes(reshaped_voltage);
    end
    imagesc(reshaped_spikes);
    title('raster of l4')
grid


figure(75)
    reshaped_spikes = reshape(spikes, n_total_neurons, length(tspan));
    imagesc(thalamic_poisson_spikes);
    title('raster of thalamic')
grid
%}



%--- most probably useless in multiple iterations -----
% --- thalamic mean spikes
