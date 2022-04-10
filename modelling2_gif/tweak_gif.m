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
close all;
clear all;

k1=0.2;k2=0.02;b=0.01;R1=0.0;R2=1.0;
El=-70.0;Vr=-70.0;Thetar=-60.0;G=0.05;C=1.0;ThetaInf=-50.0;

t_simulate = 2000;
dt = 1; % may be 1ms
tspan = 0:dt:t_simulate;

i1 = zeros(1, length(tspan));
i2 = zeros(1, length(tspan));
v = zeros(1, length(tspan));
theta = zeros(1, length(tspan));
spikes = zeros(1, length(tspan));
iext = zeros(1, length(tspan));
iext(1000:1010) = -50; % -10 is so low for the current params
% iext(500:600) = 10;


% params for phasic burst and rebound burst
% a = 0.01; A1 = 25; A2 = -0.6; 
% a = 0.01; A1 = 20; A2 = -0.6; 
%  a = 0.01; A1 = 12; A2 = -0.6;
%   a = 0.019; A1 = 5; A2 = -0.6;
%   a = 0.019; A1 = 5; A2 = -0.6;
% params fo tonic spike through spike for +ve current
% a = 0; A1 = 0; A2 = 0; 
a = 0.005; A1 = 10; A2 = -0.6;
% testing if these show rebound burst
%  a = 0.009; A1 = 5; A2 = -0.6;

% IC=(0.01, 0.001, -70.0, -50.0)
% initial conditions
i1(1) = 0.01;
i2(1) = 0.001;
v(1) =-70.0 ;
theta(1) = -50.0;

% temporrary
x = [];
for t=2:length(tspan)
    i1(t) = i1(t-1) + dt*(-k1*i1(t-1));
    i2(t) = i2(t-1) + dt*(-k2*i2(t-1));
    v(t) = v(t-1) + dt*(1/C)*(iext(t-1)+ i1(t-1) + i2(t-1) - G*(v(t-1) - El));
    theta(t) = theta(t-1) + dt*( a*(v(t-1) - El) - b*(theta(t-1) - ThetaInf)  );

    if v(t) > theta(t)
        i1(t) = R1*i1(t) + A1;
        i2(t) = R2*i2(t) + A2;
        v(t) = Vr;
        theta(t) = max(Thetar, theta(t));

        spikes(t) = 1;
       x = [x, t];
    end
end

figure
    hold on 
        plot(v)
        plot(theta)
        stem(spikes*50)
        title('voltage and theta')
        legend('v','theta','spikes')
        
    hold off
grid

a = x(1);b = length(x);c=x(end);
disp(a)
disp(c)
disp(b)
