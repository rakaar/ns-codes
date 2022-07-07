clear all;close all;
k1=0.2;k2=0.02;b=0.01;R1=0.0;R2=1.0;
El=-70.0;Vr=-70.0;Thetar=-60.0;G=0.05;C=1.0;ThetaInf=-50.0;

% a = 0.005; A1 = 10; A2 = -0.6;
% a = 0.020; A1 = 5; A2 = -0.1; % with exp inc threshold, causing large num of spikes

% positive_current = 4;
% negative_current = -25;

current = 5; 
% +ve - 2,3,4
% -ve - -15 to -24

a = 0.020; A1 = abs(current)/5; A2 = min(-0.1, -0.4 + abs(-current/10)); % with inc A2, rebound spikes decrease

disp(A2)
% a = 0.020; A1 = 0; A2 = -0.1; % phasic sort

t_simulate = 2000;
dt = 1; % may be 1ms
tspan = 0:dt:t_simulate;
