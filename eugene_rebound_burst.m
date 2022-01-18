% for i=-5:1:5
% a=0.5; b=0.5; c=-52;  d=0; % for early spike, but premature spike also hapen
% a=0.03; b=0.25; c=-52;  d=0;
% a=0.03; b=0.25; c=-52; d=0;
tic
a=0.03; b=0.25; c=-55;  d=0; % original in paper for rebound burst

% currents = [20 50 100 500 800 1200];
%for cur=1:length(currents)
    V=-64;  u=b*V;
    VV = zeros(1,100*1000+1);
    uu = zeros(1, 100*1000 + 1);
    VV=[];  uu=[];
    tau = 0.001;  tspan = 0:tau:500;
    T1=20;
    for t=tspan
        if (t>T1) & (t < T1+25) 
            I=-10;
         elseif (t >= 80) & (t <=90)
            I=0;
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

    toc
    figure(1)
            plot(tspan,VV,[0 T1 T1 (T1+100) (T1+100) max(tspan)],-85+[0 0 -5 -5 0 0]);
            axis([0 max(tspan) -90 30])
            title(['current is -',num2str(currents(1,cur))])
    grid

    
%end

