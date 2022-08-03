% % 0.03 0.25 -52 0
% %  a=0.03; b=0.25; c=-52;  d=0;
% % a=0.03; b=0.25; c=-70;  d=0;
% % a=0.02; b=0.25; c=-55; d=0.05;
% % a = 0.02; b = 0.25; c = -70; d = 0.05
% % a=0.03
% % b=0.25
% % c=-48 
% % d=0
% % a=0.02 
% % b=0.25 
% % c=-58 
% % d=0
% % a=0.02
% % b=0.25 
% % c=-58
% % d=0.5
% % a=0.03; b=0.25; c=-52;  d=0;
% a=0.03; b=0.25; c=-52;  d=0;
% V=-64;  u=b*V;
% T_end = 700;
% tau = 0.2;  
% tspan = 0:tau:T_end;
% 
% % T1=20;
% current_span = 100;
% 
% VV=zeros(1, length(tspan));  
% uu=zeros(1, length(tspan));
% 
% T2 = 300;
% current_span2 = 100;
% for t=tspan
%     if (t>T1) & (t < T1+current_span) 
%         I=0;
%     elseif (t>T2) & (t < current_span2)
% %         a=0.03; b=0.25; c=-52;  d=0;
%         I=-15;
%     else
%         I=0;
%     end;
%     V = V + tau*(0.04*V^2+5*V+140-u+I);
%     u = u + tau*a*(b*V-u);
%     if V > 30
%         VV(1, round(t/tau + 1))=30;
%         V = c;
%         u = u + d;
%     else
%         VV(1, round(t/tau + 1))=V;
%     end;
%     uu(1, round(t/tau + 1))=u;
% end;
% 
% 
% plot(tspan,VV,[0 T1 T1 (T1+current_span) (T1+current_span) T2 T2 (T2+current_span2) (T2+current_span2) max(tspan)],-85+[0 0 -5 -5 0 0 5 5 0 0]);
% axis([0 max(tspan) -90 30])
% title('(N) rebound burst and rebound burst');
% a=0.03; b=0.25; c=-52;  d=0;
a=0.03; b=0.25; c=-52;  d=0;
V=-64;  u=b*V;
VV=[];  uu=[];
tau = 0.2;  tspan = 0:tau:500;
T1=20;
for t=tspan
    if (t>T1) & (t < T1+5) 
        I=-15;
    elseif (t>200) & (t<250) 
        I=1;
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
% plot(tspan,VV,[0 T1 T1 (T1+5) (T1+5) max(tspan)],-85+[0 0 -5 -5 0 0]);
plot(tspan,VV,[0 T1 T1 (T1+5) (T1+5) T2 T2 (200+50) (200+50) max(tspan)],-85+[0 0 -5 -5 0 0 5 5 0 0]);
axis([0 max(tspan) -90 30])

title('(N) rebound burst');
