% simulate single izhikevich neurons
function [voltage_arr feedback_var_arr] = neuron(a, b, c, d, input_negative_current, input_positive_current)
    % basic variables
    T_end = 700;
    dt = 0.2; 
    tspan = 0:dt:T_end; 

    % time when negative current starts and its duration
    T1=20; negative_current_span = 200;
    % time when positive current starts and its duration
    T2 = 350; positive_current_span = 30;

    % initializing V and u
    V=-64;  u=b*V;
    VV=zeros(1, length(tspan));  
    uu=zeros(1, length(tspan));

    for t=tspan
        if (t>T1) & (t < T1+negative_current_span) 
            I=input_negative_current;
        elseif (t>T2) & (t < T2+positive_current_span)
            I=input_positive_current;
        else
            I=0;
        end;
        
        % updating voltage and feedback variable
        V = V + dt*(0.04*V^2+5*V+140-u+I);
        u = u + dt*a*(b*V-u);
        if V > 30
            VV(1, round(t/dt + 1))=30;
            V = c;
            u = u + d;
        else
            VV(1, round(t/dt + 1))=V;
        end;
        uu(1, round(t/dt + 1))=u;
    end;

    voltage_arr = VV; feedback_var_arr = uu;
end