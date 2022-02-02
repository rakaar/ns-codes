function project
    close all;
    b = 0.25;
    v_0 = -64; % source - http://www.columbia.edu/cu/appliedneuroshp/Spring2018/Spring18SHPAppliedNeuroLec5.pdf
    u_0 = b * v_0;
    [t r] = ode15s(@izhikevich_model, [0 200], [v_0 u_0]);
    figure(10)
        plot(t, r(:,1));
        title("izhikevich v vs t");
    grid

    figure(12)
        plot(t, r(:,2));
        title("izhikevich u vs t");
    grid
end

function result = izhikevich_model(t,r)
    % vars
    c1 = 0.04; c2 = 5; c3=140; c4=1; c5=1;
    a=0.03;b=0.25;c=-60;d=4;
    if t < 25 & t > 20
        i = -15;
    else
        i = 0;
    end

    % r(1) is v, r(2) is u
    if r(1) > 30
        fprintf("\n crossed 30 %f \n", r(1));
        r(1) = c;
        r(2) = r(2) + d;
    end
    fprintf("\n now r1 is  %f \n", r(1));

    result = zeros(2,1); % dv/dt, du/dt
    % result(1) = c1*(r(1)^2) + (c2*r(1)) + c3 - c4*r(2) + c5*i;
    result = zeros(2,1); % dv/dt, du/dt
    result(1) = 0.04*(r(1)^2) + 5*r(1) + 140 - r(2) + i;
    result(2) = a * ((b * r(1)) - r(2));

    

end