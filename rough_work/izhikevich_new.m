tic
a=0.03; b=0.25; c=-52;  d=0;
V=-64;  u=b*V;
T_end = 200;
tau = 0.2;  tspan = 0:tau:T_end;
T1=20;
current_span = 30;

VV=zeros(1, length(tspan));  
uu=zeros(1, length(tspan));

T2 = 150;
current_span2 = 30;

for t=tspan
    if (t>T1) & (t < T1+current_span) 
        I=-15;
    elseif (t>T2) & (t < T2+current_span2)
        I=30;
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
toc

plot(tspan,VV,[0 T1 T1 (T1+current_span) (T1+current_span) max(tspan)],-85+[0 0 -5 -5 0 0]);
axis([0 max(tspan) -90 30])
title('(N) rebound burst');
