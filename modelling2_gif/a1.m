clear all;
close all;


n_iters = 1;

% basic variables;
n_columns = 1
n_excitatory = 20; 
n_inhibitory = 5; 
n_total_neurons = n_excitatory + n_inhibitory;
n_thalamic = 9; num_of_input_giving_thalamic = 4;
    
% time step
physical_time_in_ms = 1; %dt time step 
dt = 1;  % 0.2 dt = 20 ms, so 0.01 = 1 ms 
t_simulate = 1000; 
tspan = 0:dt:t_simulate;

% making bins of 100ms = 20*dt and calculating spike rate
spike_rate_dt = 5*dt;
spike_rate_length = (length(tspan)-1)/(spike_rate_dt/dt);


% connection strength
weight_reducing_l4 = 10; % for now all weights reduced by factor of 0.2
J_ee_0 = 6*weight_reducing_l4; 
J_ie_0 = 0.5*weight_reducing_l4;
J_ei = -4*weight_reducing_l4; 
J_ii = -0.5*weight_reducing_l4;
J_ee_1 = 0.045*weight_reducing_l4; 
J_ie_1 = 0.0035*weight_reducing_l4; 
J_ee_2 = 0.015*weight_reducing_l4; 
J_ie_2 = 0.0015*weight_reducing_l4;

% voltages and terms from it are 3d tensors
voltages = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
i1_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
i2_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
theta_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
spikes = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
spike_rates = zeros(n_iters, n_columns, n_total_neurons, spike_rate_length);
recurrence_epsc_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan)-1);
epsc_tensor = zeros(n_iters, n_columns, n_total_neurons, length(tspan)-1);
I_background_tensor = zeros(n_iters, length(tspan));

% synaptic resources
xr = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
xe = zeros(n_iters, n_columns, n_total_neurons, length(tspan));
xi = zeros(n_iters, n_columns, n_total_neurons, length(tspan));

% mapping from input thalamic neurons to a1 column neurons
all_combinations = nchoosek(1:n_thalamic, num_of_input_giving_thalamic);
mapping_matrix_thalamic_to_a1 = all_combinations(1:n_total_neurons,:);

%% generate inhomo poisson spikes for thalamic neurons
thalamic_poisson_spikes = zeros(n_iters, n_thalamic, length(tspan));
lamda = zeros(1, length(tspan));
% 100ms - 3-4 spikes, 200ms - 18-20 spikes, 300 - rest - 3-4 spikes
% WARNING: FOR NOW THIS STIMULS IS HARD CODED, need to adjust acc to
% t_simulate
lamda_s = 400; lamda_i = 4;
for i=1:500
    lamda(1,i) = lamda_i;
end
for i=501:600
    lamda(1,i) = lamda_s+lamda_i;
end
for i=601:length(tspan)
    lamda(1,i) = lamda_i;
end

% calculating epsc of each thalamic neuron
weight_thalamic_to_a1 = 10; xe_thalamic = 1;
epsc_thalamic = zeros(n_iters,n_thalamic, length(tspan));

%% time constant for synaptic resources
tau_re = 0.9; tau_ir = 5000; tau_ei = 27;
% izhikevich neuron params
% for rebound burst and sustained_spike
neuron_params_rb_ss = containers.Map({'a', 'b', 'c', 'd'}, [0.02 0.25 -70 0]); 
% a=0.02; b=0.25; c=-70;  d=0; % this should also cause disinhibition
% a=0.02; b=0.25; c=-55;  d=0.05; % used this on party day in lab
% neuron_params_rb_ss = containers.Map({'a', 'b', 'c', 'd'}, [0.03 0.25 -48 0]); 

% for rebound burst and phasic spike
neuron_params_rb_ps = containers.Map({'a', 'b', 'c', 'd'}, [0.02 0.25 -58 0.5]);
    
% initialize
v0 = -70;  
xr(:, :, :, 1) = 1;
voltages(:, :, :, 1) = v0; % 
i1_tensor(:, :, :, 1) = 0.01;
i2_tensor(:, :, :, 1) = 0.001;
theta_tensor(:, :, :, 1) = -50.0;

for iter=1:n_iters
    
    fprintf("------iter numm %d -----", iter);

    % thalamic
    for i=1:n_thalamic
         thalamic_poisson_spikes(iter,i, :) = reshape(poisson_generator(lamda, dt, length(tspan)), 1, 1, length(tspan));
    end
    for i=1:n_thalamic
        epsc_thalamic(iter,i,:) = reshape(get_g_t_vector(thalamic_poisson_spikes(iter,i,:), length(tspan)) * weight_thalamic_to_a1 * xe_thalamic,  1,1,length(tspan));
    end

    % simulation
    for i=2:length(tspan)
	fprintf("i = %d\n", i);
	for c=1:n_columns
	
		for n=1:n_total_neurons
					
			% voltage sum from excitatory neighbouring columns, will be useful for inhib and exc neurons
			epsc_ex_neuron_back_c2 = 0;
			if c-2 >= 1
				for j=1:n_excitatory
					spike_train_exc = spikes(iter,c-2,j,:);
                    g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
					epsc_ex_neuron_back_c2 = epsc_ex_neuron_back_c2 + g_t*xe(iter, c-2,j,i-1); 
                 end
			end

			epsc_ex_neuron_back_c1 = 0;
			if c-1 >= 1
				for j=1:n_excitatory
					spike_train_exc = spikes(iter,c-1,j,:);
					g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
					epsc_ex_neuron_back_c1 = epsc_ex_neuron_back_c1 + g_t*xe(iter,c-1,j,i-1); 
				end
			end

			epsc_ex_neuron_front_c1 = 0;
			if c+1 <= n_columns
				for j=1:n_excitatory
				    spike_train_exc = spikes(iter,c+1,j,:);
					g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
					epsc_ex_neuron_front_c1 = epsc_ex_neuron_front_c1 + g_t*xe(iter,c+1,j,i-1);
				end
			end	


			epsc_ex_neuron_front_c2 = 0;
			if c+2 <= n_columns
				for j=1:n_excitatory
					spike_train_exc = spikes(iter,c+2,j,:);
					g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
					epsc_ex_neuron_front_c2 = epsc_ex_neuron_front_c2 + g_t*xe(iter,c+2,j,i-1);
				end
			end	


			epsc_ex_own_column = 0;
			for j=1:n_excitatory
				if j == n
					continue;
				end
				spike_train_exc = spikes(iter,c,j,:);
				g_t = get_g_t(spike_train_exc, dt, i-1, tspan);
				epsc_ex_own_column	= epsc_ex_own_column + g_t*xe(iter,c,j,i-1); 
			end

			epsc_inh_own_column = 0;
			for j=n_excitatory+1:n_total_neurons
				if j == n
					continue;
                end
                spike_train_inh = spikes(iter,c,j,:);
				g_t = get_g_t(spike_train_inh, dt, i-1, tspan);
				epsc_inh_own_column = epsc_inh_own_column + g_t*xe(iter,c,j,i-1);
            end
		
          % epsc from thalamic neurons
          epsc_from_thalamic = 0;
          for ttt=1:num_of_input_giving_thalamic
             thalamic_neuron_num = mapping_matrix_thalamic_to_a1(n, ttt);
             epsc_from_thalamic = epsc_from_thalamic + epsc_thalamic(iter,thalamic_neuron_num, i-1);
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
          recurrence_epsc_tensor(iter, c, n, i-1) = total_epsc ;
          total_epsc = total_epsc + epsc_from_thalamic;
          % a temporary statement to check if backcurrent is strong enough
          % to produce current
%           total_epsc = 0;
          
          epsc_tensor(iter, c, n, i-1) = total_epsc;
          
% 		  fprintf("c,n - %d %d total espc %f \n", c,n,total_epsc);

			
            I_background = rand * (1);
%             if i<501 | i >1500
%                 I_background = 0;
%             end
            I_background_tensor(iter, i) = I_background;
			
            
            % calculate voltage using the function
% 			[voltages(iter,c, n, i), u_values(iter,c, n, i)] = calculate_v_u(v_current, u_current, dt, neuron_params_rb_ss, total_epsc, I_background );
			[voltages(iter,c, n, i), i1_tensor(iter,c, n, i), i2_tensor(iter,c, n, i), theta_tensor(iter,c, n, i), spikes(iter,c,n,i)] = calculate_new_state(voltages(iter,c, n, i-1), i1_tensor(iter,c, n, i-1), i2_tensor(iter,c, n, i-1), theta_tensor(iter,c, n, i-1), total_epsc, I_background,dt);

			M = 0;
			if spikes(iter,c,n,i) == 1
				M = 1;
			end
		
            
            
            %	fprintf("voltage returned from function is %f \n", voltages(c,n,i));
    
            % update synaptic resources
            current_xr = xr(iter,c, n, i-1);
            current_xe = xe(iter,c, n, i-1);
            current_xi = xi(iter,c, n, i-1);

            xr(iter,c,n,i) = update_xr(M, current_xr, current_xi, tau_re, tau_ir);
            xe(iter,c, n, i) = update_xe(M, current_xr, current_xe, tau_re, tau_ei);
            xi(iter,c, n, i) = update_xi(current_xe, current_xi, tau_ei, tau_ir);
        
        end
    
   
%     fprintf("xr %f, xe %f, xi %f\n", xr(c,n,i),xe(c,n,i), xi(c,n,i));
   % pause(0.4);
    
    end

	
%	break % for testing only one iteration
end

end

%% ---- random neuron voltage plot
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
% -------------------------------
% -------- psth -------------
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

% mean psth of all neurons
spike_rate_l4_all = zeros(1, spike_rate_length);
for i=1:spike_rate_length
    spike_rate_l4_all(1, i) = sum(spike_rate_l4(:, i))/n_total_neurons;
end

figure
    plot(spike_rate_l4_all);
    title('psth of l4  all neurons')
grid
% -------- psth end -------------

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

% testing spikes of all 25 neurons
x = squeeze(spikes);
x1 = reshape(x, n_iters*n_total_neurons, length(tspan));
figure
    imagesc(x1)
grid

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
