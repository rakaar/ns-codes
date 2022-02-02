% basic variables
n_columns = 5;
n_excitatory = 100; 
n_inhibitory = 25; 
n_total_neurons = n_excitatory + n_inhibitory;
    
% time step
physical_time_in_ms = 1; %dt time step 
dt = 0.01;  % 0.2 dt = 20 ms, so 0.01 = 1 ms 
t_simulate = 1000; % x100 ms = x0.1s 
tspan = 0:dt:t_simulate;

% making bins of 100ms = 20*dt and calculating spike rate
spike_rate_dt = 100*dt;
tspan_spike_rates = 0:spike_rate_dt:t_simulate;

% voltages and terms from it are 3d tensors
voltages = zeros(n_columns, n_total_neurons, length(tspan));
spikes = zeros(n_columns, n_total_neurons, length(tspan));
spike_rates = zeros(n_columns, n_total_neurons, length(tspan_spike_rates));

% izhikevich neurons
% neuron params: 
% for rebound burst and sustained_spike
neuron_params_rb_ss = containers.Map({'a', 'b', 'c', 'd'}, [0.03 0.25 -52 0]); 
% for rebound burst and phasic spike
neuron_params_rb_ps = containers.Map({'a', 'b', 'c', 'd'}, [0.02 0.25 -58 0.5]);

disp('firing all neurons')
% example for firing a neuron, we will fire all with a bacground noise
[voltage_val, ~] = neuron_fire(dt, t_simulate, 0,0,0,0,0,0, neuron_params_rb_ss);
% reshaping and assining the voltage values 
repeated_voltage = repmat(voltage_val, [n_columns  1 n_total_neurons]);
repeated_voltage = reshape(repeated_voltage, n_columns, n_total_neurons, length(tspan));

voltages(:, :, :) = repeated_voltage;
spikes(3, 10, :) = reshape(voltage_to_spikes(voltages(3, 10, :)) ,1, 1, length(tspan));

% since all are same, i can take any neuron and repeat copies of it
disp('spike rates')

spike_rates_single_neuron = spikes_to_spike_rate(dt, spike_rate_dt, t_simulate,physical_time_in_ms, spikes(3, 10, :));
repeated_spike_rate = repmat(spike_rates_single_neuron, [n_columns  1 n_total_neurons]);
repeated_spike_rate = reshape(repeated_spike_rate, n_columns, n_total_neurons, length(tspan_spike_rates));
disp('size of spike_rates')
disp(size(spike_rates))
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

