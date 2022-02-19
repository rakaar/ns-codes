function project
    v_0 = -58; h_0 = 1;   
    [t r] = ode45(@integrate_fire_burst, [0 200], [v_0 h_0]);

    figure(15)
        T1=20;
        plot(t, r(:,1), [0 T1 T1 (T1+5) (T1+5) 200],-55+[0 0 -5 -5 0 0]);
        axis([0 200 -80 -30]);
    grid

    figure(25)
        plot(t, r(:,2));
    grid

end

function result = integrate_fire_burst(t,r)
    % vars
    c = 2;
    i1 = 0;
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

    
    T1=20;
    if (t>T1) & (t < T1+5) 
        i0=-2;
    else
        i0=0;
    end;

    v = r(1);
    h = r(2);

    % if v >= v_theta
    %     disp("here ????")
    %     v = v_reset;
    % end
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

    if v >= v_theta
        disp("here ????")
        v = v_reset;
    end

end
