a=0.02; b=0.25; c=-55;  d=0.05;
V=-64;  u=b*V;
VV=[];  uu=[];
tau = 0.01;  tspan = 0:tau:200;
T1=20;
for t=tspan
    if (t>T1) 
        I=0.6;
    else
        I=0;
    end;
    V = V + tau*(0.04*V^2+5*V+140-u+I);
    u = u + tau*a*(b*V-u);
    if V > 30
        VV(end+1)=30;
        V = c;
        u = u + d;
    else
        VV(end+1)=V;
    end;
    uu(end+1)=u;
end;
plot(tspan,VV,[0 T1 T1 max(tspan)],-90+[0 0 10 10]);
axis([0 max(tspan) -90 30])
axis off;
title('(D) phasic bursting');

