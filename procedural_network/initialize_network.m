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

% making bins of 100ms = 5*dt and calculating spike rate
spike_rate_dt = 20*dt;
tspan_spike_rates = 0:spike_rate_dt:t_simulate;

% voltages is a 3d tensor
voltages = zeros(n_columns, n_total_neurons, length(tspan));
spikes = zeros(n_columns, n_total_neurons, length(tspan));

% spike rate
spike_rates = zeros(n_columns, n_total_neurons, length(tspan_spike_rates));

% resources as per paper
resources_x_or_y = zeros(n_columns, n_total_neurons, length(tspan_spike_rates));

% neuron params: 
% for rebound burst and sustained_spike
neuron_params_rb_ss = containers.Map({'a', 'b', 'c', 'd'}, [0.03 0.25 -52 0]); 
% for rebound burst and phasic spike
neuron_params_rb_ps = containers.Map({'a', 'b', 'c', 'd'}, [0.02 0.25 -58 0.5]);

% weights for all neurons within a column
weight_matrix = zeros(n_total_neurons, n_total_neurons);

Jee = 6;
Jie = 0.5;
Jei = -4;
Jii = -0.5;

for i=1:1:n_excitatory
    for j=1:1:n_excitatory
        weight_matrix(i,j) = Jee;
    end
end

for i=1:1:n_inhibitory
    for j=1:1:n_inhibitory
        weight_matrix(i,j) = Jii;
    end
end

for i=1:1:n_excitatory
    for j=1:1:n_inhibitory
        weight_matrix(i,j) = Jei;
    end
end


for i=1:1:n_inhibitory
    for j=1:1:n_excitatory
        weight_matrix(i,j) = Jie;
    end
end

% syaptic resources for all neurons - xr, xe, xi
tau_re = 0.9; tau_ei = 27; tau_ir = 5000;
synaptic_resources = zeros(n_total_neurons, 3, t_simulate/dt + 1);

disp('firing all neurons')
% example for firing a neuron, we will fire all with a bacground noise
[voltage_val, ~] = neuron_fire(dt, t_simulate, 0,0,0,0,0,0, neuron_params_rb_ss);

% testing repmat 
repeated_voltage = repmat(voltage_val, [n_columns  1 n_total_neurons]);
repeated_voltage = reshape(repeated_voltage, n_columns, n_total_neurons, length(tspan));

voltages(:, :, :) = repeated_voltage;



figure(1)
    plot(tspan, reshape(voltages(3, 10, :), 1, length(tspan)));
    
    title('c 3 n 10 voltage')
grid

spikes(3, 10, :) = reshape(voltage_to_spikes(voltages(3, 10, :)) ,1, 1, length(tspan));

figure(2)
    stem(tspan, reshape(spikes(3, 10, :), 1, length(tspan)));
    
    title('c 3 n 10 spikes')
grid

% since all are same, i can take any neuron and repeat copies of it
disp('spike rates')

spike_rates_single_neuron = spikes_to_spike_rate(dt, spike_rate_dt, t_simulate,physical_time_in_ms, spikes(3, 10, :));
repeated_spike_rate = repmat(spike_rates_single_neuron, [n_columns  1 n_total_neurons]);
repeated_spike_rate = reshape(repeated_spike_rate, n_columns, n_total_neurons, length(tspan_spike_rates));
spike_rates(:, :, :) = repeated_spike_rate;



figure(3)
    plot(tspan_spike_rates, reshape(spike_rates(3, 10, :),  1,length(tspan_spike_rates)))
    
    title('c 3 n 10 spike rate/s')
grid

% start simulating using spike rate
t_simulate_test = t_simulate * 5;
for i=2:floor(t_simulate_test/spike_rate_dt)
    % for each column
    for c=1:n_columns
        % a)for excitatory neurons
        excitatory_input_from_neighbours = 0;
        % get the neighbours
        if c-2 >= 1
            excitatory_input_from_neighbours =  excitatory_input_from_neighbours + (J_ee_2/n_excitatory) * U * sum(spike_rates(c-2, 1:n_excitatory, i), 'all');
        end
        
        if c+2 <= n_columns
            excitatory_input_from_neighbours =  excitatory_input_from_neighbours + (J_ee_2/n_excitatory) * U * sum(spike_rates(c+2, 1:n_excitatory, i), 'all');
        end
        
        if c-1 >= 1
            excitatory_input_from_neighbours =  excitatory_input_from_neighbours + (J_ee_1/n_excitatory) * U *sum(spike_rates(c-1, 1:n_excitatory, i), 'all');
        end
        
        if c+1 <= n_columns
            excitatory_input_from_neighbours =  excitatory_input_from_neighbours + (J_ee_1/n_excitatory) * U *sum(spike_rates(c+1, 1:n_excitatory, i), 'all');
        end
        
        excitatory_input_from_its_column = (J_ee_0/n_excitatory)*sum(spike_rates(c,1:n_excitatory, i), 'all');
        
        % external inputs - e in the equation
        e_step = (bg_high_E - bg_low_E)/(n_excitatory - 1);
        external_input_to_excitatory = bg_low_E:e_step:bg_high_E;
        
        
        excitatory_input_total = excitatory_input_from_its_column + excitatory_input_from_neighbours;
        excitatory_input_vector = repmat([excitatory_input_total], 1, n_excitatory) + external_input_to_excitatory;
        non_linear_gain_excitatory_input_vector = non_linear_gain_vector(excitatory_input_vector);
        spike_rates(c, 1:n_excitatory, i) = spike_rates(c, 1:n_excitatory, i) + non_linear_gain_excitatory_input_vector;
        
      
        % inhibitory neurons
        excitatory_input_from_neighbours = 0;
        if c-2 >= 1
            excitatory_input_from_neighbours =  excitatory_input_from_neighbours + (J_ie_2/n_excitatory) * U * sum(spike_rates(c-2, 1:n_excitatory, i), 'all');
        end
        
        if c+2 <= n_columns
            excitatory_input_from_neighbours =  excitatory_input_from_neighbours + (J_ie_2/n_excitatory) * U * sum(spike_rates(c+2, 1:n_excitatory, i), 'all');
        end
        
        if c-1 >= 1
            excitatory_input_from_neighbours =  excitatory_input_from_neighbours + (J_ie_1/n_excitatory) * U *sum(spike_rates(c-1, 1:n_excitatory, i), 'all');
        end
        
        if c+1 <= n_columns
            excitatory_input_from_neighbours =  excitatory_input_from_neighbours + (J_ie_1/n_excitatory) * U *sum(spike_rates(c+1, 1:n_excitatory, i), 'all');
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
        spike_rates(c, n_excitatory+1:n_total_neurons, i) = spike_rates(c, n_excitatory+1:n_total_neurons, i) + non_linear_gain_input_to_inhibitory_vector;
    end

        
    break % only for testing one time step
end

%{

%================== testing ======================
% fire a neuron 10 in column 3
[voltage_val, ~] = neuron_fire(dt, t_simulate, 20,200,-20, 350,30,2, neuron_params_rb_ss);
voltages(3, 10, :) = reshape(voltage_val, 1, 1, length(tspan));
figure(1)
    plot(tspan, reshape(voltages(3, 10, :), 1, length(tspan)));
grid

% fire 3 more random neurons for testing purpose
random_neurons = [94 28 39];
for i=1:3
    [voltage_val, ~] = neuron_fire(dt, t_simulate, 20,200,-20, 350,30,2, neuron_params_rb_ss);
    voltages(3, random_neurons(i), :) = reshape(voltage_val, 1, 1, length(tspan));
end

% checking for 10
check_10_xe = [];
% resources depletion for 94 28 39 10
for i=1:n_total_neurons
    for j=1:length(tspan)
        M = 0;
        if voltages(3, i, j) == 30
            M = 1;
        end 
       
        
        if j==1
            x_r = 1;
            x_e = 0;
            x_i = 0;
        else
             x_r = synaptic_resources(i, 1, j-1); 
             x_e = synaptic_resources(i, 2, j-1);
             x_i = synaptic_resources(i, 3, j-1);
        end
        
        synaptic_resources(i, 1, j) = x_r + (-M*(x_r/tau_re)) + (x_i/tau_ir);
        synaptic_resources(i, 2, j) = x_e + (M*(x_r/tau_re)) - (x_e/tau_ei);
        synaptic_resources(i ,3, j) = x_i + (x_e/tau_ei) - (x_i/tau_ir);
    
        if i==94
            check_10_xe = [check_10_xe, synaptic_resources(i,2,j)];
        end
    
    end
end


figure(12)
    plot(check_10_xe)
    title('xe of 94th neuron')
grid

figure(2)
    plot(tspan, reshape(voltages(3, 28, :), 1, length(tspan)));
    title('voltage of 28thn neuron')
grid


% synapse test
% all neurons add to column 3 neuron 1

for i=2:n_total_neurons
    x = voltage_to_spikes(reshape(voltages(3, i, :), 1, length(tspan)));
    g = get_g_t(x);
    g = g(1,1:length(x));
    w = weight_matrix(1,i);
    xe = synaptic_resources(i,2, :);
    xe = reshape(xe, 1, length(x));
      
    voltages(3, 1, :) =  voltages(3, 1, :) + reshape(w*xe.*shift_1(g).*shift_1(x).*shift_1(x), 1, 1, length(tspan));
   
end

voltages(3, 1, :) = reshape(voltages(3, 1, :), 1,1,length(tspan));
% -- what about the rules of decreasing voltage ???

figure(3)
    plot(tspan, reshape(voltages(3, 1, :), 1, length(tspan)));
    title('voltage of neuron 1 ')
grid
%}