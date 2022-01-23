% basic variables
n_columns = 21;
n_excitatory = 100; 
n_inhibitory = 25; 
n_total_neurons = n_excitatory + n_inhibitory;
    
% time step
dt = 0.2;  % 20 ms as per paper https://www.izhikevich.org/publications/whichmod.pdf and code http://www.izhikevich.org/publications/figure1.m
t_simulate = 1000; % x100 ms = x0.1s 
tspan = 0:dt:t_simulate;

% voltages is a 3d tensor
voltages = zeros(n_columns, n_total_neurons, length(tspan));

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


[voltage_val, ~] = neuron_fire(dt, t_simulate, 0,0,0,0,0,0, neuron_params_rb_ss);
voltages(11, 10, :) = reshape(voltage_val, 1, 1, length(tspan));
figure(1)
    plot(tspan, reshape(voltages(11, 10, :), 1, length(tspan)));
grid
%{

%================== testing ======================
% fire a neuron 10 in column 11
[voltage_val, ~] = neuron_fire(dt, t_simulate, 20,200,-20, 350,30,2, neuron_params_rb_ss);
voltages(11, 10, :) = reshape(voltage_val, 1, 1, length(tspan));
figure(1)
    plot(tspan, reshape(voltages(11, 10, :), 1, length(tspan)));
grid

% fire 3 more random neurons for testing purpose
random_neurons = [94 28 39];
for i=1:3
    [voltage_val, ~] = neuron_fire(dt, t_simulate, 20,200,-20, 350,30,2, neuron_params_rb_ss);
    voltages(11, random_neurons(i), :) = reshape(voltage_val, 1, 1, length(tspan));
end

% checking for 10
check_10_xe = [];
% resources depletion for 94 28 39 10
for i=1:n_total_neurons
    for j=1:length(tspan)
        M = 0;
        if voltages(11, i, j) == 30
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
    plot(tspan, reshape(voltages(11, 28, :), 1, length(tspan)));
    title('voltage of 28thn neuron')
grid


% synapse test
% all neurons add to column 11 neuron 1

for i=2:n_total_neurons
    x = voltage_to_spikes(reshape(voltages(11, i, :), 1, length(tspan)));
    g = get_g_t(x);
    g = g(1,1:length(x));
    w = weight_matrix(1,i);
    xe = synaptic_resources(i,2, :);
    xe = reshape(xe, 1, length(x));
      
    voltages(11, 1, :) =  voltages(11, 1, :) + reshape(w*xe.*shift_1(g).*shift_1(x).*shift_1(x), 1, 1, length(tspan));
   
end

voltages(11, 1, :) = reshape(voltages(11, 1, :), 1,1,length(tspan));
% -- what about the rules of decreasing voltage ???

figure(3)
    plot(tspan, reshape(voltages(11, 1, :), 1, length(tspan)));
    title('voltage of neuron 1 ')
grid
%}