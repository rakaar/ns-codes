% time 
T_end = 700;
dt = 0.2; 
tspan = 0:dt:T_end;

n1 = neuron(0.03, 0.25, -52, 0);
 
negative_current=-20;
T1=20; 
negative_curent_span=200;

positive_current=2;
T2=350;
positive_current_span = 30;

n1 = n1.fire(negative_current, T1, negative_curent_span, positive_current, T2, positive_current_span);

figure(1)
    plot(tspan, n1.voltage)
grid