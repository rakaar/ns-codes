function spikes = voltage_to_spikes(voltage)
    spikes = zeros(1, length(voltage));
    voltage = reshape(voltage, 1, length(voltage));
    % 1 - exp(-(2*0.001)) =     0.0020
    for i=1:length(voltage)
        if voltage(1,i) == 30 || rand <= 0.0020
            spikes(1, i) = 1;
        end
    end

end
