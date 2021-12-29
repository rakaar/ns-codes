
% vars
c = 2;
i0 = 2;
i1 = 0; % we are not bothered about sinosidal curents
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

v_arr = [];
h_arr = [];

dt = 0.001; t_inital=0; t_final=200;
T1=20;

v_0 = -58; h_0 = 1;
v = v_0; h  = h_0;

for t=t_inital:dt:t_final
    % if (t>T1) & (t < T1+5) 
    %     I=-15;
    % else
    %     I=0;
    % end;

    v = v + dt*(...
    (1/c) * ( ...
         (i0 + i1*cos(2000*pi*f*t)) ...
         - (g_leak*(v-v_leak))  ...
         -(g_t * heaviside(v-v_h) * h * (v-v_t)) ...
         )...
    );
    if v >= v_h
        h = h + (-h/tau_h_minus);
    else
        h = h + ((1-h)/tau_h_plus);
    end

    if v  >= v_theta
        v = v_reset;
        v_arr = [v_arr, v_theta];
    else
        v_arr = [v_arr, v];
    end

    h_arr = [h_arr, h];
end

figure(12)
    plot(t_inital:dt:t_final, v_arr);
grid