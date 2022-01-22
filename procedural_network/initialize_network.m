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

% syaptic resources for all neurons - xr, xe, xi
tau_re = 0.9; tau_ei = 27; tau_ir = 5000;
synaptic_resources = zeros(125,3);
synaptic_resources(:,1) = 1;

%================== testing ======================
% fire a neuron 10 in column 11
column11 = voltages(11);
[ column11(10, :), ~] = neuron_fire(dt, t_simulate, 20,200,-20, 350,30,2, neuron_params_rb_ss);

figure(1)
    plot(tspan, column11(10, :))
grid

% fire 3 more random neurons for testing purpose
random_neurons = [94 28 39];
for i=1:3
    [column11(random_neurons(i),:), ~] = neuron_fire(dt, t_simulate, 20,200,-20, 350,30,2, neuron_params_rb_ss);
end

% checking for 10
check_10_xe = [];
% resources depletion for 94 28 39 10
for i=1:n_total_neurons
    for j=1:length(column11)
        M = 0;
        if column11(i,j) == 30
            M = 1;
        end 
        x_r = synaptic_resources(i, 1); 
        x_e = synaptic_resources(i, 2);
        x_i = synaptic_resources(i, 3);
        synaptic_resources(i, 1) = x_r + (-M*(x_r/tau_re)) + (x_i/tau_ir);
        synaptic_resources(i, 2) = x_e + (M*(x_r/tau_re)) - (x_e/tau_ei);
        synaptic_resources(i ,3) = x_i + (x_e/tau_ei) - (x_i/tau_ir);
    
    if i==94
        check_10_xe = [check_10_xe, synaptic_resources(i,2)];
    end
    
    end
end

figure(12)
    plot(check_10_xe)
    title('xe of 10th neuron')
grid

figure(2)
    plot(tspan, column11(28, :))
grid

figure(5)
    plot(synaptic_resources(:, 2)); % xe of all neurons
grid

% synapse test
% all neurons add to column 11 neuron 1
volt_1 = 0;
for i=2:n_total_neurons
    x = voltage_to_spikes(column11(i,:));
    g = get_g_t(x);
    g = g(1,1:length(x));
    w = weight_matrix(1,i);
    xe = synaptic_resources(i,2);
      
    column11(1, :) =  column11(1, :) + w*shift_1(g).*shift_1(x).*shift_1(x);
   
end
% -- what about the rules of decreasing voltage ???

figure(3)
    plot(tspan, column11(1,:))
    title('voltage of neuron 1 ')
grid