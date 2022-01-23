function [VV, uu] = neuron_fire(dt, t_simulate, tp, tp_span, positive_current, tn, tn_span, negative_current,  neuron_params)
    tspan = 0:dt:t_simulate;
    
    % params
    a = neuron_params('a'); 
    b = neuron_params('b');
    c = neuron_params('c');
    d = neuron_params('d');
    
    % initlialize 
    V=-64;  u=b*V;
    
    % noise
    I_background = 3; % generates 2.5 spikes/s
    
    VV = zeros(1, length(tspan));
    uu = zeros(1, length(tspan));
    
        for t=tspan
            if (t>tp) & (t < tp + tp_span)
                I = positive_current;
            elseif (t > tn) & (t < tn + tn_span)
                I = negative_current;
            else
                I = 0;
            end
    
             V = V + dt*(0.04*V^2+5*V+140-u+I+I_background);
             u = u + dt*a*(b*V-u);
             if V > 30
                 VV(1, round(t/dt + 1))=30;
                 V = c;
                 u = u + d;
             else
                VV(1, round(t/dt + 1))=V;
             end
            uu(1, round(t/dt + 1))=u;
    
        end
end
