% 0.03 0.25 -52 0
%  a=0.03; b=0.25; c=-52;  d=0;
a=0.02; b=0.25; c=-55;  d=0.05;

V=-64;  u=b*V;
T_end = 25;
tau = 0.01;  
tspan = 0:tau:T_end;

T1=20;
current_span = 25;

VV=zeros(1, length(tspan));  
uu=zeros(1, length(tspan));

T2 = 0.5;
current_span2 = 2.5;
for t=tspan
    if (t>T1) & (t < T1+current_span) 
        I=0;
    elseif (t>T2) & (t < current_span2)
        I=0.6;
    else
        I=0;
    end;
    V = V + tau*(0.04*V^2+5*V+140-u+I);
    u = u + tau*a*(b*V-u);
    if V > 30
        VV(1, round(t/tau + 1))=30;
        V = c;
        u = u + d;
    else
        VV(1, round(t/tau + 1))=V;
    end;
    uu(1, round(t/tau + 1))=u;
end;


plot(tspan,VV,[0 T1 T1 (T1+current_span) (T1+current_span) T2 T2 (T2+current_span2) (T2+current_span2) max(tspan)],-85+[0 0 -5 -5 0 0 5 5 0 0]);
axis([0 max(tspan) -90 30])
title('(N) rebound burst and rebound burst');
