function project

% vars
c = 2;
% i0 = 2;
i1 = 0; % we are not bothered about sinosidal curents
f=1;

g_leak = 0.035;
v_theta = -48; % try for -46, -47, -48
v_reset = -50 ;
v_leak = -25;
v_h = -60;
v_t = 250;

g_t = 0.01;

tau_h_plus = 90;
tau_h_minus = 100;

v_arr = [];
h_arr = [];

dt = 0.001; t_inital=0; t_final=200;
T1=20;

v_0 = -70; h_0 = 0.1;
v = v_0; h  = h_0;

for t=t_inital:dt:t_final
    if (t>T1) & (t < T1+5) 
        i0=-2;
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
    plot(t_inital:dt:t_final, v_arr, [0 T1 T1 (T1+5) (T1+5) t_final],-55+[0 0 -5 -5 0 0]);
    axis([0 t_final -80 -30]);
grid

figure(13)
    hold on
    scatter(t_inital:dt:t_final,h_arr);
    plot([0 T1 T1 (T1+5) (T1+5) t_final], [0 0 -0.2 -0.2 0 0]);
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