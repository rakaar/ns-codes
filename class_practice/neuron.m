classdef Neuron
    properties
        a = 0;
        voltage_arr = zeros(1,2);
    end
    
    methods
        function obj = Neuron(a_val)
            obj.a = a_val;
        end
        
        function voltage_values = fire_neuron()
            voltage_arr = zeros(1,2);
        end
    end

end