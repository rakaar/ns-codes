classdef column
    properties
        N_excitatory = 100;
        N_inhibitory = 25;
        
        weights;
    end
    
    methods
        function obj = column(weights)
            obj.weights = weights;
        end
    end
end