% basic variables
n_columns = 21;
n_excitatory = 100; 
n_inhibitory = 25; 
n_total_neurons = n_excitatory + n_inhibitory;
    
% time step
dt = 0.2; 
t_simulate = 700;
tspan = 0:dt:t_simulate;

% Voltage: key is column number, value is voltage values
voltages = containers.Map('KeyType','int32', 'ValueType','any');
for k=1:1:n_columns
    voltages(k) = zeros(n_total_neurons, t_simulate/dt + 1);
end

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


% fire a neuron 10 in column 11
column11 = voltages(11);
[ column11(10, :), ~] = neuron_fire(dt, t_simulate, 20,200,-20, 350,30,2, neuron_params_rb_ss);

figure(1)
    plot(tspan, column11(10, :))
grid


