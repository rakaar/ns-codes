function poisson_spike_train = poisson_generator(lambda, dt)
    % lambda is spike rate
    period = length(lambda);
    poisson_spike_train = zeros(1, period);
    

    rng('shuffle');
    for i=1:period
        x = rand;
        if x <= 1 - exp(-(lambda(i)*0.001))
            poisson_spike_train(1,i) = 1;
        end
    end
end