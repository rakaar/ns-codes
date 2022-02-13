% basic variables
n_columns = 1;
n_excitatory = 20; 
n_inhibitory = 5; 
n_total_neurons = n_excitatory + n_inhibitory;
n_thalamic = 9; num_of_input_giving_thalamic = 4;
    
% time step
physical_time_in_ms = 1; %dt time step 
dt = 0.01;  % 0.2 dt = 20 ms, so 0.01 = 1 ms 
t_simulate = 10; % x100 ms = x0.1s 
tspan = 0:dt:t_simulate;

% making bins of 100ms = 20*dt and calculating spike rate
spike_rate_dt = 5*dt;
tspan_spike_rates = 0:spike_rate_dt:t_simulate;

% currents
I = 0; I_background = 2;

% connection strength
J_ee_0 = 6; J_ie_0 = 0.5;
J_ei = -4; J_ii = -0.5;
J_ee_1 = 0.045; J_ie_1 = 0.0035; J_ee_2 = 0.015; J_ie_2 = 0.0015;

% voltages and terms from it are 3d tensors
voltages = zeros(n_columns, n_total_neurons, length(tspan));
u_values = zeros(n_columns, n_total_neurons, length(tspan));
spikes = zeros(n_columns, n_total_neurons, length(tspan));
spike_rates = zeros(n_columns, n_total_neurons, length(tspan_spike_rates));
epsc_tensor = zeros(n_columns, n_total_neurons, length(tspan)-1);
epsc_thalamic = zeros(n_columns, n_total_neurons, length(tspan)-1);

% synaptic resources
xr = zeros(n_columns, n_total_neurons, length(tspan));
xe = zeros(n_columns, n_total_neurons, length(tspan));
xi = zeros(n_columns, n_total_neurons, length(tspan));

% mapping from input thalamic neurons to a1 column neurons
all_combinations = nchoosek(1:n_thalamic, num_of_input_giving_thalamic);
mapping_matrix_thalamic_to_a1 = all_combinations(1:n_total_neurons,:);

% generate inhomo poisson spikes for thalamic neurons
thalamic_poisson_spikes = zeros(n_thalamic, length(tspan));
lamda = zeros(1, length(tspan));
% 100ms - 3-4 spikes, 200ms - 18-20 spikes, 300 - rest - 3-4 spikes
% WARNING: FOR NOW THIS STIMULS IS HARD CODED
lamda_s = 3; lamda_i = 20;
for i=1:100
    lamda(1,i) = lamda_s;
end
for i=101:200
    lamda(1,i) = lamda_s+lamda_i;
end
for i=201:length(tspan)
    lamda(1,i) = lamda_s;
end
% lamda_s = 3+25 ; % tweaked so that a1 neurons have 3-4 spikes/s

for i=1:n_thalamic
    thalamic_poisson_spikes(i, :) = poisson_generator(lamda, dt, length(tspan));
end

% time constant for synaptic resources
tau_re = 0.9; tau_ir = 5000; tau_ei = 27;
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

% disp('firing all neurons')
% % example for firing a neuron, we will fire all with a bacground noise
% [voltage_val, ~] = neuron_fire(dt, t_simulate, 0,0,0,0,0,0, neuron_params_rb_ss);
% % reshaping and assining the voltage values 
% repeated_voltage = repmat(voltage_val, [n_columns  1 n_total_neurons]);
% repeated_voltage = reshape(repeated_voltage, n_columns, n_total_neurons, length(tspan));
% 
% voltages(:, :, :) = repeated_voltage;
% spikes(3, 10, :) = reshape(voltage_to_spikes(voltages(3, 10, :)) ,1, 1, length(tspan));
% spikes = repmat(spikes(3,10,:), [n_columns  1 n_total_neurons]);
% spikes = reshape(spikes, n_columns, n_total_neurons, length(tspan));
% % since all are same, i can take any neuron and repeat copies of it
% disp('spike rates')
% 
% spike_rates_single_neuron = spikes_to_spike_rate(dt, spike_rate_dt, t_simulate,physical_time_in_ms, spikes(3, 10, :));
% repeated_spike_rate = repmat(spike_rates_single_neuron, [n_columns  1 n_total_neurons]);
% repeated_spike_rate = reshape(repeated_spike_rate, n_columns, n_total_neurons, length(tspan_spike_rates));
% spike_rates(:, :, :) = repeated_spike_rate;



% figure(1)
%      plot(tspan, reshape(voltages(3, 10, :), 1, length(tspan)));
%     
%     title('c 3 n 10 voltage')
%  grid
% 
% figure(2)
%     stem(tspan, reshape(spikes(3, 10, :), 1, length(tspan)));
%     
%     title('c 3 n 10 spikes')
% grid
% 
% figure(3)
%     plot(tspan_spike_rates, reshape(spike_rates(3, 10, :),  1,length(tspan_spike_rates)))
%     
%     title('c 3 n 10 spike rate/s')
% grid

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
                    
%                     if (g_t*xe(c-2,j,i-1) ~= 0)
%                         fprintf("non0 epscontri %f, c %d, j %d", g_t*xe(c-2,j,i-1),c,j);
%                     end
                    
                    
				end
			end

			epsc_ex_neuron_back_c1 = 0;
			if c-1 >= 1
				for j=1:n_excitatory
					spike_train_exc = voltage_to_spikes(voltages(c-1,j,:));
					g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
					epsc_ex_neuron_back_c1 = epsc_ex_neuron_back_c1 + g_t*xe(c-1,j,i-1); 
                    
%                      if (g_t*xe(c-1,j,i-1) ~= 0)
%                         fprintf("non0 epscontri %f, c %d, j %d", g_t*xe(c-1,j,i-1),c,j);
%                     end
				end
			end

			epsc_ex_neuron_front_c1 = 0;
			if c+1 <= n_columns
				for j=1:n_excitatory
					spike_train_exc = voltage_to_spikes(voltages(c+1,j,:));
					g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
					epsc_ex_neuron_front_c1 = epsc_ex_neuron_front_c1 + g_t*xe(c+1,j,i-1);
                    
%                     if (g_t*xe(c+1,j,i-1) ~= 0)
%                         fprintf("non0 epscontri %f, c %d, j %d ", g_t*xe(c+1,j,i-1),c,j);
%                     end
                    
				end
			end	


			epsc_ex_neuron_front_c2 = 0;
			if c+2 <= n_columns
				for j=1:n_excitatory
					spike_train_exc = voltage_to_spikes(voltages(c+2,j,:));
					g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
					epsc_ex_neuron_front_c2 = epsc_ex_neuron_front_c2 + g_t*xe(c+2,j,i-1);
                    
%                     if (g_t*xe(c+2,j,i-1) ~= 0)
%                         fprintf("non0 epscontri %f, c %d, j %d ", g_t*xe(c+2,j,i-1),c,j);
%                     end
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
                
%                 if (g_t*xe(c,j,i-1) ~= 0)
%                         fprintf("non0 epscontri %f, c %d, j %d ", g_t*xe(c,j,i-1),c,j);
%                 end
			end

			epsc_inh_own_column = 0;
			for j=n_excitatory+1:n_total_neurons
				if j == n
					continue;
				end
				spike_train_inh = voltage_to_spikes(voltages(c,j,:));
				g_t = get_g_t(spike_train_inh, dt, i-1, tspan);
				epsc_inh_own_column = epsc_inh_own_column + g_t*xe(c,j,i-1);
                
%                 if (g_t*xe(c,j,i-1) ~= 0)
%                         fprintf("non0 epscontri %f, c %d, j %d\n", g_t*xe(c,j,i-1),c,j);
%                 end
                
            end
		
          % epsc from thalamic neurons
          weight_thalamic_to_a1 = 0.2; xe_thalamic = 1;
          epsc_from_thalamic = 0;
          
          for ttt=1:num_of_input_giving_thalamic
             thalamic_neuron_num = mapping_matrix_thalamic_to_a1(n, ttt);
             g_t = get_g_t(thalamic_poisson_spikes(thalamic_neuron_num, :), dt, i-1, tspan);
             epsc_from_thalamic = epsc_from_thalamic + g_t*weight_thalamic_to_a1*xe_thalamic;
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
          epsc_thalamic(c, n, i-1) = epsc_from_thalamic;
          
% 		  fprintf("c,n - %d %d total espc %f \n", c,n,total_epsc);
% 			
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

testing_column = 1;
testing_neuron = 6;
thalamic_testing_neuron = 6;

figure(36894)
    stem(thalamic_poisson_spikes(thalamic_testing_neuron,:));
    title(' 5 thalamic poisson spikes')
grid

% epsc due to 5th thalamic neuron
epsc_thalamic_5 = zeros(1, length(tspan));

 for i=2:floor(t_simulate/dt)
    g_t = get_g_t(thalamic_poisson_spikes(thalamic_testing_neuron, :), dt, i-1, tspan);
    epsc_thalamic_5(1,i) = epsc_thalamic_5(1,i) + g_t*weight_thalamic_to_a1*xe_thalamic;
end
figure(53643)
   plot(epsc_thalamic_5)
   title('5 epsc thalamic')
grid

figure(355343)
    hold on
        stem(thalamic_poisson_spikes(thalamic_testing_neuron,:));
        plot(epsc_thalamic_5);
    hold off
grid


% population behaviour - mean spike rate of column with time
allneurons_spike_rates = zeros(n_total_neurons, length(tspan_spike_rates));
for i=1:n_total_neurons
     spikes1 = voltage_to_spikes(voltages(1, i, :));
     spike_rates1 = spikes_to_spike_rate(dt, spike_rate_dt, t_simulate, physical_time_in_ms, spikes1);
     allneurons_spike_rates(i,:) = spike_rates1;
end
population_psth = zeros(1, length(tspan_spike_rates));
for i=1:length(tspan_spike_rates)
    population_psth(1,i) = sum(allneurons_spike_rates(:,i))/n_total_neurons;
end
figure(4546)
    plot( population_psth);
    title('population psth');
grid

% thalamic neurons - mean spike rate
thalamic_neurons_spikes_rates = zeros(n_thalamic, length(tspan_spike_rates));
for i=1:n_thalamic
   spikes11 = poisson_generator(lamda, dt, length(tspan));
   spike_rates11 = spikes_to_spike_rate(dt, spike_rate_dt, t_simulate, physical_time_in_ms, spikes11);
   thalamic_neurons_spikes_rates(i,:) = spike_rates11;
end
thalamic_population_psth = zeros(1, length(tspan_spike_rates));
for i=1:length(tspan_spike_rates)
    thalamic_population_psth(1,i) = sum(thalamic_neurons_spikes_rates(:,i))/n_thalamic;
end
figure(46578)
    plot(thalamic_population_psth);
    title('thalamic populaton psth')
grid


figure(1991)
    reshaped_epsc = reshape(epsc_tensor(testing_column, testing_neuron, :), 1, length(tspan)-1);
    plot(reshaped_epsc);
    title('epsc');
grid

figure(199374)
    reshaped_epsc_thalamic = reshape(epsc_thalamic(testing_column, testing_neuron, :), 1, length(tspan)-1);
    plot(reshaped_epsc_thalamic);
    title('epsc thalamic');
grid

%checking if epsc of all neurons is same or not
epsc_neurons_sum = zeros(1,n_total_neurons);
for i=1:n_total_neurons
    epsc_neurons_sum(1,i) = sum(epsc_thalamic(1, i, :));
end

figure(46753)
    stem(1:n_total_neurons, epsc_neurons_sum);
    title('epsc thalamic of all neurons 1-ntotal')
grid

figure(23111)
    reshaped_xe = reshape(xe(testing_column,testing_neuron,:), 1, length(tspan));
    plot(reshaped_xe);
    title('xe');
grid

figure(23112)
    reshaped_xr = reshape(xr(testing_column,testing_neuron,:), 1, length(tspan));
     plot(reshaped_xr);
    title('xr');
grid

figure(23113)
    reshaped_xi = reshape(xi(testing_column,testing_neuron,:), 1, length(tspan));
     plot(reshaped_xi);
    title('xi');
grid

figure(100)
	reshaped_voltage = reshape(voltages(testing_column, testing_neuron, :), 1, length(tspan));
	plot(tspan, reshaped_voltage);
	title('voltage')
grid

figure(345)
    spikes1 = voltage_to_spikes(voltages(testing_column, testing_neuron, :));
    stem(tspan, spikes1);
    title('spikes')
grid

figure(1234)
    spikes1 = voltage_to_spikes(voltages(testing_column, testing_neuron, :));
	spike_rates1 = spikes_to_spike_rate(dt, spike_rate_dt, t_simulate, physical_time_in_ms, spikes1);
    plot(tspan_spike_rates, spike_rates1);
    title('spike rate')
grid

figure(87556)
    spike_rates1 = spikes_to_spike_rate(dt, spike_rate_dt, t_simulate, physical_time_in_ms, thalamic_poisson_spikes(thalamic_testing_neuron,:));
    plot(spike_rates1);
    title('spike rate thal')
grid

% for i=1:n_total_neurons
% 	spikes1 = voltage_to_spikes(voltages(1, i, :));
% 	spike_rates1 = spikes_to_spike_rate(dt, spike_rate_dt, t_simulate, physical_time_in_ms, spikes1);
% 
% 	fprintf("mean sprate %d is %f \n", i, sum(spike_rates1)/length(spike_rates1));
% end


% calculate mean spike rate of all columns


x = zeros(1,25);	
for i=1:25
     v = voltages(1,i,:);
    s = voltage_to_spikes(v);
    sss = sum(s);
    x(1,i) = sss;
end
figure(34)
    stem(1:25, x)
grid
   
   
    
% 	for n=1:n_total_neurons
%         voltage_neuron = voltages(1,n,:);
% 		spikes_neuron = voltage_to_spikes(voltage_neuron);
% 		spike_rate_neuron = spikes_to_spike_rate(dt, spike_rate_dt, t_simulate,physical_time_in_ms, spikes_neuron);
%         len_of_spike_rate = length(spike_rate_neuron);
%         sum_of_spike_rates = sum(spike_rate_neuron, 'all');
%         mean_spike_rates(1,n) = sum_of_spike_rates/len_of_spike_rate;
%     end
		
    
%     
%     figure(12)
%         plot(1:n_total_neurons, ,mean_spike_rates);
%         title("neurons-x, mean spikerate - y")
%     grid



%{
spikes_random_neuron = spikes(3, 5, 7);
spike_rate_random = spikes_to_spike_rate(dt, spike_rate_dt, t_simulate, physical_time_in_ms, spikes_random_neuron);

figure(100)
	plot(spike_rate_random);
	title('spike rate of random neuron')
grid

figure(101)
	plot(spikes_random_neuron);
	title("spikes of random neuron");
grid
%}