% parameters
% 	# Define simulation and model parameters dictionary
% params = {'k1': 0.2,             # I₁ decay term
%           'k2': 0.02,            # I₂ decay term
%           'b': 0.01,             # Θ decay term
%           'R1': 0.0,             # I₁ update constant
%           'R2': 1.0,             # I₂ update constant
%           'El': -70.0,           # Reverse potential
%           'Vr': -70.0,           # Resting potential
%           'Thetar': -60.0,       # Resting threshold
%           'a': 0.000,            # Specifier of spike response
%           'A1': 0.0,             # I₁ additive update constant
%           'A2': 0.0,             # I₂ additive update constant
%           'G': 0.05,             # Membrane potential decay term
%           'C': 1.0,              # Membrane capacitance
%           'ThetaInf': -50.0}       # Reverse threshold
close all;clear all;

k1=0.2;k2=0.02;b=0.01;R1=0.0;R2=1.0;
El=-70.0;Vr=-70.0;Thetar=-60.0;G=0.05;C=1.0;ThetaInf=-50.0;

iext_val = -25;
% iext_val = -5;
% a = 0.005; A1 = 10; A2 = -0.6;
% a = 0.020; A1 = 5; A2 = -0.1; % with exp inc threshold, causing large num of spikes

a = 0.020; A1 = 5; A2 = -0.15; % with inc A2, rebound spikes decrease

% a = 0.020; A1 = 0; A2 = -0.1; % phasic sort

t_simulate = 2000;
dt = 1; % may be 1ms
tspan = 0:dt:t_simulate;

i1 = zeros(1, length(tspan));
i2 = zeros(1, length(tspan));
v = zeros(1, length(tspan));
theta = zeros(1, length(tspan));
spikes = zeros(1, length(tspan));
iext = zeros(1, length(tspan));
iext(950:1000) = iext_val;
% iext(950:1000) = 5; % 2,3,4

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
        title('voltage and threshold and spikes')
    hold off
grid

a = x(1);b = length(x);c=x(end);
disp(a)
disp(c)
disp(b)


