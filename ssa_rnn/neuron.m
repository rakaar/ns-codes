classdef neuron
    properties
        % params
        a;b;c;d;
        
        % basic variables
        T_end = 700;
        dt = 0.2; 
       
        % v and u
        voltage;
        feedback_var;
    end
    
    methods
        function obj = neuron(a,b,c,d)
            obj.a = a;
            obj.b = b;
            obj.c = c;
            obj.d = d;
        end
        
        function obj = fire(obj, negative_current, T1, negative_current_span, positive_current, T2, positive_current_span)
            
            a=obj.a; b=obj.b; c=obj.c; d=obj.d;
            
            % initializing V and u
            V=-64;  u=b*V;
            
            dt = obj.dt;
            tspan = 0:dt:obj.T_end;
           
            
            
            voltage = zeros(1, length(tspan));
            feedback_var = zeros(1, length(tspan));
     
            disp(obj.T_end)
            disp(max(tspan))
            for t=tspan
                if (t>T1) & (t < T1+negative_current_span) 
                    I = negative_current;
                elseif (t>T2) & (t < T2+positive_current_span)
                    I = positive_current;
                else
                    I=0;
                end;

                % updating voltage and feedback variable
                V = V + dt*(0.04*V^2+5*V+140-u+I);
                u = u + dt*a*(b*V-u);
                if V > 30
                    voltage(1, round(t/dt + 1))=30;
                    V = c;
                    u = u + d;
                else
                    voltage(1, round(t/dt + 1))=V;
                end
                feedback_var(1, round(t/dt + 1))=u;
            end
            
            disp("voltage size")
            disp(size(voltage))
            
            obj.voltage = voltage; obj.feedback_var = feedback_var;

            

        end
    end
    
end