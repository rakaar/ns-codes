a=0.03; b=0.25; c=-52;  d=0;% original in paper for rebound burst

T = 500;
dt = 0.2;  tspan = 0:dt:T;
T1=20; % time when current introducded
current_span = 100; % time for which current is introduced

positive_current_t_inital = 300;
positive_current_t_final = 350;

VV = zeros(1,length(tspan));
uu = zeros(1,length(tspan));

% initial conditions
V=-64;  u=b*V;

for t=tspan
    if (t>T1) & (t < T1+current_span) 
        I=-5;
    elseif (t >= positive_current_t_inital) & (t <= positive_current_t_final)
        I=0;
    else
        I=0;
    end;
    V = V + dt*(0.04*V^2+5*V+140-u+I);
    u = u + dt*a*(b*V-u);
    if V > 30
        VV(1, fix((t+dt)/dt))=30;
        V = c;
        u = u + d;
    else
        VV(1,fix((t+dt)/dt))=V;
    end;
    uu(1, fix((t+dt)/dt))=u;
end;

figure(1)
    plot(tspan,VV,[0 T1 T1 (T1+current_span) (T1+current_span) T],-85+[0 0 -5 -5 0 0]);
    axis([0 max(tspan) -90 30])
    title(['current is -',num2str(I)])
grid

    

