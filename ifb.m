function project
    v_0 = -58; h_0 = 0;   
    [t r] = ode15s(@integrate_fire_burst, [0 200], [v_0 h_0]);

    figure(1)
        plot(t,r(:,1));
    grid

end

function result = integrate_fire_burst(t,r)
    % vars
    c = 2;
    i0 = 2;
    i1 = 1;
    f=1;

    g_leak = 0.035;
    v_theta = -35;
    v_reset = -50;
    v_leak = -65;
    v_h = -60;
    v_t = 120;

    g_t = 0.07;

    tau_h_plus = 100;
    tau_h_minus = 20;

    v = r(1);
    h = r(2);

    if v >= v_theta
        disp("here ????")
        v = v_reset;
    end
    result = zeros(2,1); 
    result(1) = (1/c) * ( ...
         (i0 + i1*cos(2000*pi*f*t)) ...
         - (g_leak*(v-v_leak))  ...
         -(g_t * heaviside(v-v_h) * h * (v-v_t)) ...
         );
    if v >= v_h
        result(2) = -h/tau_h_minus;
    else
        result(2) = (1-h)/tau_h_plus;
    end
end
