i1 = zeros(1, length(tspan));
i2 = zeros(1, length(tspan));
v = zeros(1, length(tspan));
theta = zeros(1, length(tspan));
spikes = zeros(1, length(tspan));
iext = zeros(1, length(tspan));
iext(950:1000) = positive_current; % 2,3,4

% params for phasic burst and rebound burst
% a = 0.009; A1 = 15; A2 = -0.6; 

% params fo tonic spike through spike for +ve current
% a = 0; A1 = 0; A2 = 0; 

% 
% a = 0.005; A1 = 10; A2 = -0.6;
% 
% % adjusting for realistic voltage values
% C = 3;
% b = 0.001;
% a = 0.01;

% IC=(0.01, 0.001, -70.0, -50.0)
% initial conditions
i1(1) = 0.01;
i2(1) = 0.001;
v(1) =-70.0 ;
theta(1) = -50.0;
latest_spike_time = 1;
% temporrary
x = [];
for t=2:length(tspan)
    i1(t) = i1(t-1) + dt*(-k1*i1(t-1));
%     if t > 940 && t < 960
%         disp("-----------------")
%         disp(dt*(-k1*i1(t-1)))
%     end
    i2(t) = i2(t-1) + dt*(-k2*i2(t-1));
    v(t) = v(t-1) + dt*(1/C)*(iext(t-1)+ i1(t-1) + i2(t-1) - G*(v(t-1) - El));
%     if t > 1000 && t < 1050
%         disp(dt*(1/C)*(iext(t-1)+ i1(t-1) + i2(t-1) - G*(v(t-1) - El)))
%         disp(i1(t-1))
%         disp(i2(t-1))
%         disp(i2(t-1))
%         disp(-G*(v(t-1) - El))
%         disp('------------------------')
%     end
    theta(t) = theta(t-1) + dt*( a*(v(t-1) - El) - b*(theta(t-1) - ThetaInf)  );
    time_constant = 5;
    if t - latest_spike_time <= 5 && latest_spike_time ~= -1
        theta(t) = theta(t-1) + (Thetar-theta(t-1))*exp(-(t+1-latest_spike_time)/time_constant);
    end
    v_limit = -100;
    if v(t) < v_limit
        v(t)= v_limit;
    end
    
    if v(t) > theta(t)
        i1(t) = R1*i1(t) + A1;
        i2(t) = R2*i2(t) + A2;
        v(t) = Vr;
%         theta(t) = max(Thetar, theta(t));
        theta(t) = theta(t-1) + (Thetar-theta(t-1))*exp(-(t+1-latest_spike_time)/time_constant);
        spikes(t) = 1;
        latest_spike_time = t;
       x = [x, t];
    end

    if spikes(t-1) == 1 && spikes(t) == 1
        spikes(t) = 0;
    end
end

figure
    hold on
        plot(v)
        plot(theta)
        plot(spikes*100)
        plot(i1)
        plot(i2)
        plot(iext)
        legend('v', 'theta', 'spikes','i1','i2', 'iext')
        title('phasic burst - voltage and threshold and spikes')
    hold off
grid

disp("phasic burst")
a = x(1);b = length(x);c=x(end);
disp(a)
disp(c)
disp(b)
