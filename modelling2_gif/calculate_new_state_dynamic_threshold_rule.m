function [v, i1, i2, theta, is_spike] = calculate_new_state_dynamic_threshold_rule(v_old, i1_old, i2_old, theta_old,total_epsc, I_background,dt,t,t_spike)
     is_spike = 0;
    
%      background element causing spike
% some how adding this background element results prevents disinhibtion
%      if rand <= 0.0025
%         is_spike = 1;
%      end

    
    k1=0.2;k2=0.02;b=0.01;R1=0.0;R2=1.0;
    El=-70.0;Vr=-70.0;Thetar=-60.0;G=0.05;C=1.0;ThetaInf=-50.0;

    a = 0.020;




iext = total_epsc + I_background;

     i1 = i1_old + dt*(-k1*i1_old);
    i2 = i2_old + dt*(-k2*i2_old);
    v = v_old + dt*(1/C)*(iext+ i1_old + i2_old - G*(v_old - El));
    

    % reaching threshold slowly
    % current time - last spike <= 5
    
    theta = theta_old + dt*( a*(v_old - El) - b*(theta_old - ThetaInf)  );
     


    time_constant = 5;
    if t - t_spike <= 5 && t_spike ~= -1
        theta = theta_old + (Thetar-theta_old)*exp(-(t+1-t_spike)/time_constant);
    end
    
     % abs/rel refractory period
    if t - t_spike >= 1 && t - t_spike <= 5
        v = v_old + 100*exp(-(t-t_spike)/1);
    end

    if v >= theta
        A1 = abs(iext)/5;
        A2 = min(-0.1, -0.4 + abs(-iext/10)); 
        i1 = R1*i1 + A1;
        i2 = R2*i2 + A2;
        %         v = Vr; % w/o refractory period
        v = Vr - 50; % abs/rel refractory period
        theta = theta_old + (Thetar-theta_old)*exp(-(t+1-t_spike)/time_constant);
        is_spike = 1;
    end

    

    % biological limits on threshold and voltage values
%     if v < -100
%         v = -100;
%     end
% 
%     if theta > 60
%         theta = 60;
%     end
end