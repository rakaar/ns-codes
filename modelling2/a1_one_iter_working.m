close all;

% basic variables
n_columns = 1;
n_excitatory = 20; 
n_inhibitory = 5; 
n_total_neurons = n_excitatory + n_inhibitory;
n_thalamic = 9; num_of_input_giving_thalamic = 4;
    
% time step
physical_time_in_ms = 1; %dt time step 
dt = 0.01;  % 0.2 dt = 20 ms, so 0.01 = 1 ms 
t_simulate = 20; % x100 ms = x0.1s 
tspan = 0:dt:t_simulate;

% making bins of 100ms = 20*dt and calculating spike rate
spike_rate_dt = 5*dt;
tspan_spike_rates = 0:spike_rate_dt:t_simulate;

% currents
I_background = 8.0;

% connection strength
weight_reducing_l4 = 5; % for now all weights reduced by factor of 0.2
J_ee_0 = 6*weight_reducing_l4; 
J_ie_0 = 0.5*weight_reducing_l4;
J_ei = -4*weight_reducing_l4; 
J_ii = -0.5*weight_reducing_l4;
J_ee_1 = 0.045*weight_reducing_l4; 
J_ie_1 = 0.0035*weight_reducing_l4; 
J_ee_2 = 0.015*weight_reducing_l4; 
J_ie_2 = 0.0015*weight_reducing_l4;

% voltages and terms from it are 3d tensors
voltages = zeros(n_columns, n_total_neurons, length(tspan));
u_values = zeros(n_columns, n_total_neurons, length(tspan));
spikes = zeros(n_columns, n_total_neurons, length(tspan));
spike_rates = zeros(n_columns, n_total_neurons, length(tspan_spike_rates));
epsc_tensor = zeros(n_columns, n_total_neurons, length(tspan)-1);

% synaptic resources
xr = zeros(n_columns, n_total_neurons, length(tspan));
xe = zeros(n_columns, n_total_neurons, length(tspan));
xi = zeros(n_columns, n_total_neurons, length(tspan));

% mapping from input thalamic neurons to a1 column neurons
all_combinations = nchoosek(1:n_thalamic, num_of_input_giving_thalamic);
mapping_matrix_thalamic_to_a1 = all_combinations(1:n_total_neurons,:);

%% generate inhomo poisson spikes for thalamic neurons
thalamic_poisson_spikes = zeros(n_thalamic, length(tspan));
lamda = zeros(1, length(tspan));
% 100ms - 3-4 spikes, 200ms - 18-20 spikes, 300 - rest - 3-4 spikes
% WARNING: FOR NOW THIS STIMULS IS HARD CODED, need to adjust acc to
% t_simulate
lamda_s = 200; lamda_i = 4;
for i=1:500
    lamda(1,i) = lamda_i;
end
for i=501:1500
    lamda(1,i) = lamda_s+lamda_i;
end
for i=1501:length(tspan)
    lamda(1,i) = lamda_i;
end
for i=1:n_thalamic
    thalamic_poisson_spikes(i, :) = poisson_generator(lamda, dt, length(tspan));
end
% calculating epsc of each thalamic neuron
weight_thalamic_to_a1 = 250; xe_thalamic = 1;
epsc_thalamic = zeros(n_thalamic, length(tspan));
for i=1:n_thalamic
    epsc_thalamic(i,:) = get_g_t_vector(thalamic_poisson_spikes(i,:), length(tspan)) * weight_thalamic_to_a1 * xe_thalamic;
end

%% time constant for synaptic resources
tau_re = 0.9; tau_ir = 5000; tau_ei = 27;
% izhikevich neuron params
% for rebound burst and sustained_spike
neuron_params_rb_ss = containers.Map({'a', 'b', 'c', 'd'}, [0.03 0.25 -52 0]); 
% neuron_params_rb_ss = containers.Map({'a', 'b', 'c', 'd'}, [0.03 0.25 -48 0]); 

% for rebound burst and phasic spike
neuron_params_rb_ps = containers.Map({'a', 'b', 'c', 'd'}, [0.02 0.25 -58 0.5]);

% initialize
v0 = -64; u0 = neuron_params_rb_ss('b')*v0; 
xr(:, :, 1) = 1;
voltages(:, :, 1) = v0; % 
u_values(:, :, 1) = u0; % 


for i=2:length(tspan)
	fprintf("i = %d\n", i);
	for c=1:n_columns
	
		for n=1:n_total_neurons
					
			% voltage sum from excitatory neighbouring columns, will be useful for inhib and exc neurons
			epsc_ex_neuron_back_c2 = 0;
			if c-2 >= 1
				for j=1:n_excitatory
					spike_train_exc = voltage_to_spikes(voltages(c-2,j,:));
					g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
					epsc_ex_neuron_back_c2 = epsc_ex_neuron_back_c2 + g_t*xe(c-2,j,i-1); 
                 end
			end

			epsc_ex_neuron_back_c1 = 0;
			if c-1 >= 1
				for j=1:n_excitatory
					spike_train_exc = voltage_to_spikes(voltages(c-1,j,:));
					g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
					epsc_ex_neuron_back_c1 = epsc_ex_neuron_back_c1 + g_t*xe(c-1,j,i-1); 
				end
			end

			epsc_ex_neuron_front_c1 = 0;
			if c+1 <= n_columns
				for j=1:n_excitatory
					spike_train_exc = voltage_to_spikes(voltages(c+1,j,:));
					g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
					epsc_ex_neuron_front_c1 = epsc_ex_neuron_front_c1 + g_t*xe(c+1,j,i-1);
				end
			end	


			epsc_ex_neuron_front_c2 = 0;
			if c+2 <= n_columns
				for j=1:n_excitatory
					spike_train_exc = voltage_to_spikes(voltages(c+2,j,:));
					g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
					epsc_ex_neuron_front_c2 = epsc_ex_neuron_front_c2 + g_t*xe(c+2,j,i-1);
				end
			end	


			epsc_ex_own_column = 0;
			for j=1:n_excitatory
				if j == n
					continue;
				end
				spike_train_exc = voltage_to_spikes(voltages(c,j,:));
				g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
				epsc_ex_own_column	= epsc_ex_own_column + g_t*xe(c,j,i-1); 
			end

			epsc_inh_own_column = 0;
			for j=n_excitatory+1:n_total_neurons
				if j == n
					continue;
				end
				spike_train_inh = voltage_to_spikes(voltages(c,j,:));
				g_t = get_g_t(spike_train_inh, dt, i-1, tspan);
				epsc_inh_own_column = epsc_inh_own_column + g_t*xe(c,j,i-1);
            end
		
          % epsc from thalamic neurons
          epsc_from_thalamic = 0;
          for ttt=1:num_of_input_giving_thalamic
             thalamic_neuron_num = mapping_matrix_thalamic_to_a1(n, ttt);
             epsc_from_thalamic = epsc_from_thalamic + epsc_thalamic(thalamic_neuron_num, i-1);
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
          total_epsc = total_epsc + epsc_from_thalamic;
          epsc_tensor(c, n, i-1) = total_epsc;
          
% 		  fprintf("c,n - %d %d total espc %f \n", c,n,total_epsc);

			v_current = voltages(c,n,i-1);
			u_current = u_values(c, n, i-1);
			if v_current == 30
                v_current = neuron_params_rb_ss('c');
				u_current = u_current + neuron_params_rb_ss('d');
			end
			% calculate voltage using the function
			[voltages(c, n, i), u_values(c, n, i)] = calculate_v_u(v_current, u_current, dt, neuron_params_rb_ss, total_epsc, I_background );
						
			M = 0;
			if voltages(c, n, i) == 30
				M = 1;
			end
		%	fprintf("voltage returned from function is %f \n", voltages(c,n,i));
    
            % update synaptic resources
            current_xr = xr(c, n, i-1);
            current_xe = xe(c, n, i-1);
            current_xi = xi(c, n, i-1);

            xr(c,n,i) = update_xr(M, current_xr, current_xi, tau_re, tau_ir);
            xe(c, n, i) = update_xe(M, current_xr, current_xe, tau_re, tau_ei);
            xi(c, n, i) = update_xi(current_xe, current_xi, tau_ei, tau_ir);
        
        end
    
   
%     fprintf("xr %f, xe %f, xi %f\n", xr(c,n,i),xe(c,n,i), xi(c,n,i));
   % pause(0.4);
    
    end

	
%	break % for testing only one iteration
end

testing_column = 1;testing_neuron = 15;thalamic_testing_neuron = 6;

figure
    vvv = reshape(voltages(1, testing_neuron, :),1, length(tspan));
    plot(vvv);
    title('voltage of neuron')
grid

figure
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

% thalamic spikes psth
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

figure
    stem(thalamic_population_psth)
    title('psth of thalmic')
grid

figure
    plot(population_psth);
    title('psth of l4')
grid


% spikes of l4 neurons
spikes_2d_matrix_l4 = zeros(n_total_neurons, length(tspan));
for i=1:n_total_neurons
    spikes_2d_matrix_l4(i,:) = voltage_to_spikes(voltages(1,i,:));
end

figure
    imagesc(spikes_2d_matrix_l4);
    title('raster l4')
grid

figure
    imagesc(thalamic_poisson_spikes);
    title('raster of thalamic')
grid


% figure
%     spikes1 = voltage_to_spikes(voltages(testing_column, testing_neuron, :));
%     stem(tspan, spikes1);
%     title('spikes')
% grid
%  
% figure
%     spikes1 = voltage_to_spikes(voltages(testing_column, testing_neuron, :));
% 	spike_rates1 = spikes_to_spike_rate_neat(spikes1, physical_time_in_ms, dt, spike_rate_dt);
%     plot(spike_rates1);
%     title('spike rate')
% grid

% psth of l4 neurons
% ff=[];
% for ii=1:25
%     qw=spikes_2d_matrix_l4(ii,1:2500);
%     ff=[ff;mean(reshape(qw,5,500))];
% end
% figure
%     plot(mean(ff))
%     title('psth of l4')
% grid




