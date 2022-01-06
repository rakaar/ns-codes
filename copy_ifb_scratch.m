function project

% vars
c = 2;
% i0 = 2;
i1 = 0; % we are not bothered about sinosidal curents
f=1;

g_leak = 0.035;
v_theta = -35; % try for -46, -47, -48
v_reset = -50 ;
v_leak = -40;
v_h = -45;
v_t = 120;

g_t = 0.07;

tau_h_plus = 40;
tau_h_minus = 80;

v_arr = [];
h_arr = [];

dt = 0.001; t_inital=0; t_final=200;
T1=20;

v_0 = v_leak; h_0 = 0.05;
v = v_0; h  = h_0;

for t=t_inital:dt:t_final
    if (t>T1) & (t < T1+20) 
        % i0=-0.3;
        i0 = -0.25;
    elseif (t>=45) & (t<=55)
        i0=0;
    else
        i0=0;
    end;

    v = v + dt*(...
    (1/c) * ( ...
         (i0 + i1*cos(2000*pi*f*t)) ...
         - (g_leak*(v-v_leak))  ...
         -(g_t * my_heaviside(v,v_h) * h * (v-v_t)) ...
         )...
    );
    if v >= v_h
        h = h + dt*(-h/tau_h_minus);
    else
        h = h + dt*((1-h)/tau_h_plus);
    end

    if v  >= v_theta
        v = v_reset;
        v_arr = [v_arr, 30];
    else
        v_arr = [v_arr, v];
    end

    h_arr = [h_arr, h];
end

figure(12)
    plot(t_inital:dt:t_final, v_arr, [0 T1 T1 (T1+20) (T1+20) t_final],-55+[0 0 -5 -5 0 0]);
    axis([0 t_final -80 -30]);
grid

figure(13)
    hold on
    scatter(t_inital:dt:t_final,h_arr);
    plot([0 T1 T1 (T1+5) (T1+5) t_final], [0 0 -0.02 -0.02 0 0]);
    hold off
grid



end % project end

function heaviside_value = my_heaviside(v,v_h)
    if v >= v_h
        heaviside_value = 1;
    else
        heaviside_value = 0;
    end
end