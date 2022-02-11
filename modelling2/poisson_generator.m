function poisson_spike_train = poisson_generator(lambda_s, dt, period)
    % lambda_s is spontaneous spike rate
    poisson_spike_train = zeros(1, period);
    rng('shuffle');
    for i=1:period
        x = rand;
        if x <= 1 - exp(-(lambda_s*dt))
            poisson_spike_train(1,i) = 1;
        end
    end
end