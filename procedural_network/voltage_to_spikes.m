function spikes = voltage_to_spikes(voltage)
    spikes = zeros(1, length(voltage));
    voltage = reshape(voltage, 1, length(voltage));
    
    for i=1:length(voltage)
        if voltage(1,i) == 30
            spikes(1, i) = 1;
        end
    end

end