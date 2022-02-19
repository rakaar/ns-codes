% for i=-5:1:5
a=0.02; b=0.05; c=-60;  d=30;


    V=-64;  u=b*V;
    VV=[];  uu=[];
    tau = 0.2;  tspan = 0:tau:100;
    T1=20;
    for t=tspan
        if (t>T1) & (t < T1+5) 
            I=-15;

        else
            I=0;
        end;

        if t >= T1
            a=0.5; b=0.5; c=-52;  d=0;
        else
            a=0.02; b=0.01; c=-60; d=30;
        end
        

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
        % title([' value ', num2str(i)]);
        title('voltage vs t,(blue - voltage, red - current)')
    grid
% end
