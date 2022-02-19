function [v_new u_new] = calculate_v_u(v, u, dt, neuron_params, epsc, I_background) 
    
    % params
    a = neuron_params('a'); 
    b = neuron_params('b');
    c = neuron_params('c');
    d = neuron_params('d');

    v_new  = v + dt*(0.04*v^2 + 5*v + 140 - u + epsc + I_background);
	u_new = u + dt*a*(b*v-u);
    if v >= 30
         v_new = 30;
    else
end
