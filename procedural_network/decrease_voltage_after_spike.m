function new_voltage = decrease_voltage_after_spike(original_voltage)
    new_voltage = reshape(original_voltage, 1, length(original_voltage));
    threshold = 0.05; t_spike = 0; tau_ref = 2;beta = 5;
    
    for i=1:length(original_voltage)
        if new_voltage(1, i) >= threshold
            t_spike = i;
            for j=i:i+19
                new_voltage(1,i) = new_voltage(1,i)*(1-(beta*exp(-(j-t_spike)/tau_ref)));
            end
            
            i = j;
        end
    end
    
end