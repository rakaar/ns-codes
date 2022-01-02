% for i=-5:1:5
% a=0.5; b=0.5; c=-52;  d=0; % for early spike, but premature spike also hapen
 a=0.03; b=0.25; c=-60;  d=0; % original in paper for rebound burst
% a=0.03; b=0.25; c=-52;  d=0;
% a=0.03; b=0.25; c=-52; d=0;
 V=-64;  u=b*V;
    VV=[];  uu=[];
    tau = 0.001;  tspan = 0:tau:100;
    T1=20;
    for t=tspan
        if (t>T1) & (t < T1+5) 
            I=-1200;
        elseif (t >= 70) & (t <=80)
            I=10;
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
    figure(floor(100+1))
            plot(tspan,VV,[0 T1 T1 (T1+5) (T1+5) max(tspan)],-85+[0 0 -5 -5 0 0]);
            axis([0 max(tspan) -90 30])
            title('voltage vs t,(blue - voltage, red - current)')
    grid
% end
