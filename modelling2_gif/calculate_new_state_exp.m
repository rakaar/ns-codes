function [v, i1, i2, theta, is_spike] = calculate_new_state_exp(v_old, i1_old, i2_old, theta_old,total_epsc, I_background,dt)
     is_spike = 0;
    
%      background element causing spike
% some how adding this background element results prevents disinhibtion
%      if rand <= 0.001
%         is_spike = 1;
%      end

    k1=0.2;k2=0.02;b=0.01;R1=0.0;R2=1.0;
    El=-70.0;Vr=-70.0;Thetar=-40.0;G=0.05;C=1.0;ThetaInf=-50.0;

%     a = 0.009; A1 = 10; A2 = -0.6; 
%     a = 0.020; A1 = 10; A2 = -0.6;
%     a = 0.090; A1 = 15; A2 = -0.6; % to get rebound burst very soon
%       a = 0.2; A1 = 25; A2 = -0.6; % to get rebound burst for a small period of negative current 
%     a = 0.01; A1 = 25; A2 = -0.6; % with the above params, we get only in the beginning spike
%      a = 0.01; A1 = 12; A2 = -0.6;
% a = 0.01; A1 = 25; A2 = -0.6;
%   a = 0.005; A1 = 10; A2 = -0.6;
  a = 0.005; A1 = 12; A2 = -0.6;
iext = total_epsc + I_background;

     i1 = i1_old + dt*(-k1*i1_old);
    i2 = i2_old + dt*(-k2*i2_old);
    v = v_old + dt*(1/C)*(iext+ i1_old + i2_old - G*(v_old - El));
    theta = theta_old + dt*( a*(v_old - El) - b*(theta_old - ThetaInf)  );

    if v > theta
        i1 = R1*i1 + A1;
        i2 = R2*i2 + A2;
        v = Vr;
        theta = max(Thetar, theta);
        is_spike = 1;
    end

    % if its a spike in prev bin, v should be vr
%     if v_old == Vr
%         is_spike = 0;
%     end
end