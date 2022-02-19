% basic variables
n_columns = 5;
n_excitatory = 100; 
n_inhibitory = 25; 
n_total_neurons = n_excitatory + n_inhibitory;
    
% time step
physical_time_in_ms = 20; %dt time step 
dt = 0.2;  % 20 ms as per paper https://www.izhikevich.org/publications/whichmod.pdf and code http://www.izhikevich.org/publications/figure1.m
t_simulate = 1000; % x100 ms = x0.1s 
tspan = 0:dt:t_simulate;

% making bins of 100ms = 20*dt and calculating spike rate
spike_rate_dt = 20*dt;
tspan_spike_rates = 0:spike_rate_dt:t_simulate;


% Background input drawn from a uniform distribution (neurons are indexed by input strength):
bg_low_E  = -9.9; % lowest background input (in Hz) to excitatory population
bg_high_E = 9.9; % highest background input (in Hz) to excitatory population
bg_low_I  = bg_low_E; % lowest background input (in Hz) to inhibitory population
bg_high_I = bg_high_E; % highest background input (in Hz) to inhibitory population

% time constants 
tau_E     = 0.001; % excitatory neurons' time constant (in seconds)
tau_I     = 0.001; % inhibitory neurons' time constant (in seconds) 
tau_ref_E = 0.003; % tau refractory of excitatory neurons (in seconds)
tau_ref_I = 0.003; % tau refractory of inhibitory neurons (in seconds)
tau_rec   = 0.800; % recovery time constant of intracortical synapses (in seconds)
tau_rec_s = 0.300; % recovery time constant of sensory input synapses (in seconds)

% connection strength
J_ee_0 = 6; J_ie_0 = 0.5;
J_ei = -4; J_ii = -0.5;
J_ee_1 = 0.045; J_ie_1 = 0.0035; J_ee_2 = 0.015; J_ie_2 = 0.0015;

% utilization fraction
U = 0.5; Us = 0.7;


% voltages and terms from it are 3d tensors
voltages = zeros(n_columns, n_total_neurons, length(tspan));
spikes = zeros(n_columns, n_total_neurons, length(tspan));
spike_rates = zeros(n_columns, n_total_neurons, length(tspan_spike_rates));

% resources as per paper
resources_x_or_y = zeros(n_columns, n_total_neurons, length(tspan_spike_rates));

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

% start simulating using spike rate
for i=2:floor(t_simulate/spike_rate_dt)
    % for each column
    for c=1:n_columns
        % a)for excitatory neurons
        excitatory_input_from_neighbours = 0;
        % get the neighbours
        if c-2 >= 1
            excitatory_input_from_neighbours =  excitatory_input_from_neighbours + (J_ee_2/n_excitatory) * U * sum(resources_x_or_y(c-2, 1:n_excitatory, i).*spike_rates(c-2, 1:n_excitatory, i), 'all');
        end
        
        if c+2 <= n_columns
            excitatory_input_from_neighbours =  excitatory_input_from_neighbours + (J_ee_2/n_excitatory) * U * sum(resources_x_or_y(c+2, 1:n_excitatory, i).*spike_rates(c+2, 1:n_excitatory, i), 'all');
        end
        
        if c-1 >= 1
            excitatory_input_from_neighbours =  excitatory_input_from_neighbours + (J_ee_1/n_excitatory) * U *sum(resources_x_or_y(c-1, 1:n_excitatory, i).*spike_rates(c-1, 1:n_excitatory, i), 'all');
        end
        
        if c+1 <= n_columns
            excitatory_input_from_neighbours =  excitatory_input_from_neighbours + (J_ee_1/n_excitatory) * U *sum(resources_x_or_y(c+1, 1:n_excitatory, i).*spike_rates(c+1, 1:n_excitatory, i), 'all');
        end
        
        excitatory_input_from_its_column = (J_ee_0/n_excitatory)*U*sum(resources_x_or_y(c, 1:n_excitatory, i).*spike_rates(c,1:n_excitatory, i), 'all');
        
        inhibitory_input_to_excitatory_own_column = (J_ei/n_inhibitory)*U*sum(resources_x_or_y(c, n_excitatory+1:n_total_neurons, i).*spike_rates(c, n_excitatory+1:n_total_neurons, i), 'all');
        
        % external inputs - e in the equation
        e_step = (bg_high_E - bg_low_E)/(n_excitatory - 1);
        external_input_to_excitatory = bg_low_E:e_step:bg_high_E;
        
        
        excitatory_input_total = excitatory_input_from_its_column + excitatory_input_from_neighbours + inhibitory_input_to_excitatory_own_column;
        excitatory_input_vector = repmat([excitatory_input_total], 1, n_excitatory) + external_input_to_excitatory;
        non_linear_gain_excitatory_input_vector = non_linear_gain_vector(excitatory_input_vector);
        
        tau_ref_e_term_vector = 1 - tau_ref_E*spike_rates(c, 1:n_excitatory, i);
        
        spike_rates(c, 1:n_excitatory, i) = spike_rates(c, 1:n_excitatory, i-1) ...
                                            + spike_rate_dt*((1/tau_E)*(-spikes(c, 1:n_excitatory, i) + tau_ref_e_term_vector.*non_linear_gain_excitatory_input_vector));
        
      
        % b) inhibitory neurons
        excitatory_input_from_neighbours = 0;
        if c-2 >= 1
            excitatory_input_from_neighbours =  excitatory_input_from_neighbours + (J_ie_2/n_excitatory) *  sum(spike_rates(c-2, 1:n_excitatory, i), 'all');
        end
        
        if c+2 <= n_columns
            excitatory_input_from_neighbours =  excitatory_input_from_neighbours + (J_ie_2/n_excitatory) *  sum(spike_rates(c+2, 1:n_excitatory, i), 'all');
        end
        
        if c-1 >= 1
            excitatory_input_from_neighbours =  excitatory_input_from_neighbours + (J_ie_1/n_excitatory) * sum(spike_rates(c-1, 1:n_excitatory, i), 'all');
        end
        
        if c+1 <= n_columns
            excitatory_input_from_neighbours =  excitatory_input_from_neighbours + (J_ie_1/n_excitatory) * sum(spike_rates(c+1, 1:n_excitatory, i), 'all');
        end
        
        % external input to inhibitory neurons
        i_step = (bg_high_I - bg_low_I)/(n_inhibitory - 1);
        external_input_to_inhibitory = bg_low_I:i_step:bg_high_I;
        
        excitatory_input_from_its_column = (J_ie_0/n_excitatory)*sum(spike_rates(c,1:n_excitatory, i), 'all');
        
        excitatory_input_total = excitatory_input_from_its_column + excitatory_input_from_neighbours;
        
        inhibitory_input_from_its_column = (J_ii/n_inhibitory)*sum(spike_rates(c, n_excitatory+1:n_total_neurons, i), 'all');
        input_to_inhibitory = inhibitory_input_from_its_column + excitatory_input_total;
        input_to_inhibitory_vector = repmat([input_to_inhibitory], 1, n_inhibitory) + external_input_to_inhibitory;
        non_linear_gain_input_to_inhibitory_vector = non_linear_gain_vector(input_to_inhibitory_vector);
        
        tau_ref_i_term_vector = 1 - tau_ref_I*spike_rates(c, n_excitatory+1:n_total_neurons, i);
        
        spike_rates(c, n_excitatory+1:n_total_neurons, i) = spike_rates(c, n_excitatory+1:n_total_neurons, i-1) ...
                                                            + spike_rate_dt* ((1/tau_I)*(-spike_rates(c, n_excitatory+1:n_total_neurons, i) + tau_ref_i_term_vector.*non_linear_gain_input_to_inhibitory_vector));
    
        
        % resources update
        resources_x_or_y(c, 1:n_excitatory, i) = resources_x_or_y(c, 1:n_excitatory, i-1) ...
                                                   + spike_rate_dt*((1-resources_x_or_y(c,1:n_excitatory,i))/tau_ref_E  - U*resources_x_or_y(c, 1:n_excitatory, i).*spike_rates(c, 1:n_excitatory, i));
                                               
       resources_x_or_y(c, n_excitatory+1:n_total_neurons, i) = resources_x_or_y(c, n_excitatory+1:n_total_neurons, i-1) ...
                                                   + spike_rate_dt*((1-resources_x_or_y(c,n_excitatory+1:n_total_neurons,i))/tau_ref_I  - U*resources_x_or_y(c, n_excitatory+1:n_total_neurons, i).*spike_rates(c, n_excitatory+1:n_total_neurons, i));
                                               
        
    end

        
    % break % only for testing one time step
end

figure(11)
    plot(reshape(spike_rates(3, 99, :)  , 1,length(tspan_spike_rates)))
    title('3, 99 spike rate')
grid

figure(12)
    plot(reshape(spike_rates(2, 124, :)  , 1,length(tspan_spike_rates)))
    title('2, 124 spike rate')
grid


figure(13)
    plot(reshape(spike_rates(1, 30, :)  , 1,length(tspan_spike_rates)))
    title('1, 30 spike rate')
grid
