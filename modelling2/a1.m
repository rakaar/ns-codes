% basic variables
n_columns = 5;
n_excitatory = 20; 
n_inhibitory = 5; 
n_total_neurons = n_excitatory + n_inhibitory;
    
% time step
physical_time_in_ms = 1; %dt time step 
dt = 0.01;  % 0.2 dt = 20 ms, so 0.01 = 1 ms 
t_simulate = 1000; % x100 ms = x0.1s 
tspan = 0:dt:t_simulate;

% making bins of 100ms = 20*dt and calculating spike rate
spike_rate_dt = 5*dt;
tspan_spike_rates = 0:spike_rate_dt:t_simulate;

% currents
I = 0; I_background = 1;

% connection strength
J_ee_0 = 6; J_ie_0 = 0.5;
J_ei = -4; J_ii = -0.5;
J_ee_1 = 0.045; J_ie_1 = 0.0035; J_ee_2 = 0.015; J_ie_2 = 0.0015;

% voltages and terms from it are 3d tensors
voltages = zeros(n_columns, n_total_neurons, length(tspan));
u_values = zeros(n_columns, n_total_neurons, length(tspan));
spikes = zeros(n_columns, n_total_neurons, length(tspan));
spike_rates = zeros(n_columns, n_total_neurons, length(tspan_spike_rates));

% synaptic resources
xr = zeros(n_columns, n_total_neurons, length(tspan));
xe = zeros(n_columns, n_total_neurons, length(tspan));
xi = zeros(n_columns, n_total_neurons, length(tspan));

% izhikevich neuron params
% for rebound burst and sustained_spike
neuron_params_rb_ss = containers.Map({'a', 'b', 'c', 'd'}, [0.03 0.25 -52 0]); 
% for rebound burst and phasic spike
neuron_params_rb_ps = containers.Map({'a', 'b', 'c', 'd'}, [0.02 0.25 -58 0.5]);

% initialize 
v0 = -64; u0 = neuron_params_rb_ss('b')*v0; 
xr(:, :, 1) = 1;
voltages(:, :, 1) = v0; % 
u_values(:, :, 1) = u0; % 

%{
disp('firing all neurons')
% example for firing a neuron, we will fire all with a bacground noise
[voltage_val, ~] = neuron_fire(dt, t_simulate, 0,0,0,0,0,0, neuron_params_rb_ss);
% reshaping and assining the voltage values 
repeated_voltage = repmat(voltage_val, [n_columns  1 n_total_neurons]);
repeated_voltage = reshape(repeated_voltage, n_columns, n_total_neurons, length(tspan));

voltages(:, :, :) = repeated_voltage;
spikes(3, 10, :) = reshape(voltage_to_spikes(voltages(3, 10, :)) ,1, 1, length(tspan));
spikes = repmat(spikes(3,10,:), [n_columns  1 n_total_neurons]);
spikes = reshape(spikes, n_columns, n_total_neurons, length(tspan));
% since all are same, i can take any neuron and repeat copies of it
disp('spike rates')

spike_rates_single_neuron = spikes_to_spike_rate(dt, spike_rate_dt, t_simulate,physical_time_in_ms, spikes(3, 10, :));
repeated_spike_rate = repmat(spike_rates_single_neuron, [n_columns  1 n_total_neurons]);
repeated_spike_rate = reshape(repeated_spike_rate, n_columns, n_total_neurons, length(tspan_spike_rates));
spike_rates(:, :, :) = repeated_spike_rate;



figure(1)
     plot(tspan, reshape(voltages(3, 10, :), 1, length(tspan)));
    
    title('c 3 n 10 voltage')
 grid

figure(2)
    stem(tspan, reshape(spikes(3, 10, :), 1, length(tspan)));
    
    title('c 3 n 10 spikes')
grid

figure(3)
    plot(tspan_spike_rates, reshape(spike_rates(3, 10, :),  1,length(tspan_spike_rates)))
    
    title('c 3 n 10 spike rate/s')
grid
%}

for i=2:floor(t_simulate/dt)
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
		
		  total_epsc = 0;
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


			% calculate voltage using the function
			[voltages(c, n, i) u_values(c, n, i)] = calculate_v_u(voltages(c, n, i-1) ,u_values(c,n,i-1), dt, neuron_params_rb_ss, total_epsc, I_background );
			
			if voltages(c, n, i) == 30
				voltages(c, n, i+1) = neuron_params_rb_ss('c');
				u_values(c, n, i+1) = u_values(c, n, i) + neuron_params_rb_ss('d');
			end
	
		end
	end
	break % for testing only one iteration
end

