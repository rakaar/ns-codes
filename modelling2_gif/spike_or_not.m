function is_spike = spike_or_not(voltage_val)
    if voltage_val == 30 || rand <= 0.0020
            is_spike = 1;
    end
end
