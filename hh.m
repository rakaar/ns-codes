function project

    alpha_m = (-0.1 * (-60 + 35))/(exp(-(-60 + 35)/10) - 1);
    beta_m = 4;

    h_fixed = 0.596121;
    alpha_n = (0.1)/(exp(1) - 1);
    beta_n = 0.125; 

    alpha_h = 0.07 ;
    beta_h = 1/(1 + exp(3));

    m_eq = alpha_m/(alpha_m + beta_m);
    h_eq = alpha_h/(alpha_h + beta_h);
    n_eq = alpha_n/(alpha_n + beta_n);

    fprintf("m h n are %f %f %f \n", m_eq, h_eq, n_eq);

     % vars
     g_k_bar = 36;
     e_k = -72;
 
     g_na_bar = 120;
     e_na = 55;
 
     g_l = 0.3;
     e_l = -49.401079;
     
 
     c = 1;
     
     v_eq = -60;

     syms v m h n
     X = vpasolve([
        (g_k_bar * (n^4) * (v - e_k))  + (g_na_bar * (m^3) * h * (v - e_na)) + (g_l * (v - e_l)) == 0,
        ((-0.1 * (v+35))/(exp(-(v+35)/10) -1))*(1-m) - (4 * exp(-(v+60)/18))*(m) == 0,
        (0.07 * exp(-(v+60)/20))*(1-h) - (1/(exp(-(v+30)/10) + 1))*(h) == 0,
        ((-0.01 * (v+50))/(exp(-(v+50)/10) - 1))*(1-n) - (0.125 * (exp(-(v+60)/80)))*(n) == 0
    ], [v,m,h,n]);
    disp("equilibrium points for hh ");
    fprintf("v m h n %f %f %f %f \n", X.v, X.m, X.h, X.n);


    syms v1 m1 n1
     X1 = vpasolve([
        (g_k_bar * (n1^4) * (v1 - e_k))  + (g_na_bar * (m1^3) * 1 * (v1 - e_na)) + (g_l * (v1 - e_l)) == 0,
        ((-0.1 * (v1+35))/(exp(-(v1+35)/10) -1))*(1-m1) - (4 * exp(-(v1+60)/18))*(m1) == 0,
        ((-0.01 * (v1+50))/(exp(-(v1+50)/10) - 1))*(1-n1) - (0.125 * (exp(-(v1+60)/80)))*(n1) == 0
    ], [v1,m1,n1]);
    disp("equilibrium points for hh ");
    fprintf("v m h n %f %f %f %f \n", X1.v1, X1.m1,  X1.n1);

    syms v2 m2 n2
     X2 = vpasolve([
        (g_k_bar * (n2^4) * (v2 - e_k))  + (g_na_bar * (m2^3) * h_fixed * (v2 - e_na)) + (g_l * (v2 - e_l)) == 0,
        ((-0.1 * (v2+35))/(exp(-(v2+35)/10) -1))*(1-m2) - (4 * exp(-(v2+60)/18))*(m2) == 0,
        ((-0.01 * (v2+50))/(exp(-(v2+50)/10) - 1))*(1-n2) - (0.125 * (exp(-(v2+60)/80)))*(n2) == 0
    ], [v2,m2,n2]);
    disp("equilibrium points for hh ");
    fprintf("v m h n %f %f %f %f \n", X2.v2, X2.m2,  X2.n2);


    figure(13)
        subplot(5,1,1)
        [t_h,r_h] = ode15s(@hh4,[0 200],[double(X.v), double(X.m), double(X.h), double(X.n)]);
        plot(t_h, r_h(:,1));
        title('v vs t')

        subplot(5,1,2)
        plot(t_h, r_h(:,3));
        title('h vs t')

        subplot(5,1,3)
        [t_h,r_h] = ode15s(@hh4_without_h,[0 200],[double(X1.v1), double(X1.m1), double(X1.n1)]);
        plot(t_h, r_h(:,1));
        title('without h(h=1), v vs t')

        subplot(5,1,4)
        [t_h,r_h] = ode15s(@hh4_with_fixed_h,[0 200],[double(X2.v2), double(X2.m2), double(X2.n2)]);
        plot(t_h, r_h(:,1));
        title('without h(h=0.596121), v vs t')

        subplot(5,1,5)
        plot([0 40 40 70 70 200], [0 0 1 1 0 0]);
        title('current')
    grid

    figure(15)
        [t_h,r_h] = ode15s(@hh4_without_h,[0 200],[double(X1.v1), double(X1.m1), double(X1.n1)]);
        subplot(3,1,1)
        plot(t_h, r_h(:,1));
        title('h=1 - v vs t')

        subplot(3,1,2)
        plot(t_h, r_h(:,2));
        title('h=1 -m vs t')
        

        subplot(3,1,3)
        plot(t_h, r_h(:,3));
        title('h=1 -n vs t')
    grid

    figure(14)
        [t_h,r_h] = ode15s(@hh4,[0 200],[double(X.v), double(X.m), double(X.h), double(X.n)]);
        subplot(4,1,1)
        plot(t_h, r_h(:,1));
        title('regular-v vs t')

        subplot(4,1,2)
        plot(t_h, r_h(:,2));
        title('regular-m vs t')

        subplot(4,1,3)
        plot(t_h, r_h(:,3));
        title('regular-h vs t')

        subplot(4,1,4)
        plot(t_h, r_h(:,4));
        title('regular-n vs t')
    grid

    %% only fixed h
    
    figure(21)
        [t_h,r_h] = ode15s(@hh4_with_fixed_h,[0 200],[double(X2.v2), double(X2.m2), double(X2.n2)]);
        subplot(3,1,1)
        plot(t_h, r_h(:,1));
        title('fixed h-v vs t')

        subplot(3,1,2)
        plot(t_h, r_h(:,2));
        title('fixed h-m vs t')

        subplot(3,1,3)
        plot(t_h, r_h(:,3));
        title('fixed h-n vs t')

        
    grid

    %% only fixed h
end % end of project

function result = hh4(t,r)

    % vars
    g_k_bar = 36;
    e_k = -72;

    g_na_bar = 120;
    e_na = 55;

    g_l = 0.3;
    e_l = -49.401079;

    c = 1;

    if t >=40 & t<=70
        iext=50;
    else
        iext = 0;
    end

   

    if r(1) == -35
        alpha_m = 1;
    else
        alpha_m = (-0.1 * (r(1) + 35))/(exp(-(r(1) + 35)/10) - 1);
    end
    beta_m = 4 * exp(-(r(1) + 60)/18);


    if r(1) == -50 
        alpha_n = 0.1;
    else
        alpha_n = (-0.01 * (r(1) + 50))/(exp(-(r(1) + 50)/10) - 1);
    end
    beta_n = 0.125 * exp(-(r(1) + 60)/80);

    alpha_h = 0.07 * exp(-(r(1) + 60)/20);
    beta_h = 1/(1 + exp(-(r(1)+30)/10));
    
    result = zeros(4,1); % v,m,h,n
    result(1) = (1/c) * ( iext - (g_k_bar * r(4)^4 * (r(1) - e_k)) - (g_na_bar * r(2)^3 * r(3) * (r(1) - e_na)) - (g_l * (r(1) - e_l)) );
    result(2) = (alpha_m * (1 - r(2))) - (beta_m * r(2));   
    result(3) = (alpha_h * (1 - r(3))) - (beta_h * r(3));
    result(4) = (alpha_n * (1 - r(4))) - (beta_n * r(4));
end

function result = hh4_without_h(t,r)

    % vars
    g_k_bar = 36;
    e_k = -72;

    g_na_bar = 120;
    e_na = 55;

    g_l = 0.3;
    e_l = -49.401079;

    c = 1;

    if t >=40 & t<=70
        iext=50;
    else
        iext = 0;
    end

   

    if r(1) == -35
        alpha_m = 1;
    else
        alpha_m = (-0.1 * (r(1) + 35))/(exp(-(r(1) + 35)/10) - 1);
    end
    beta_m = 4 * exp(-(r(1) + 60)/18);


    if r(1) == -50 
        alpha_n = 0.1;
    else
        alpha_n = (-0.01 * (r(1) + 50))/(exp(-(r(1) + 50)/10) - 1);
    end
    beta_n = 0.125 * exp(-(r(1) + 60)/80);

    
    
    result = zeros(3,1); % v,m,n
    result(1) = (1/c) * ( iext - (g_k_bar * r(3)^4 * (r(1) - e_k)) - (g_na_bar * r(2)^3 * 1 * (r(1) - e_na)) - (g_l * (r(1) - e_l)) );
    result(2) = (alpha_m * (1 - r(2))) - (beta_m * r(2));   
    result(3) = (alpha_n * (1 - r(3))) - (beta_n * r(3));
end

function result = hh4_with_fixed_h(t,r)
    h_fixed = 0.15;
    % vars
    g_k_bar = 36;
    e_k = -72;

    g_na_bar = 120;
    e_na = 55;

    g_l = 0.3;
    e_l = -49.401079;

    c = 1;

    if t >=40 & t<=70
        iext=50;
    else
        iext = 0;
    end

   

    if r(1) == -35
        alpha_m = 1;
    else
        alpha_m = (-0.1 * (r(1) + 35))/(exp(-(r(1) + 35)/10) - 1);
    end
    beta_m = 4 * exp(-(r(1) + 60)/18);


    if r(1) == -50 
        alpha_n = 0.1;
    else
        alpha_n = (-0.01 * (r(1) + 50))/(exp(-(r(1) + 50)/10) - 1);
    end
    beta_n = 0.125 * exp(-(r(1) + 60)/80);

    
    
    result = zeros(3,1); % v,m,n
    result(1) = (1/c) * ( iext - (g_k_bar * r(3)^4 * (r(1) - e_k)) - (g_na_bar * r(2)^3 * h_fixed * (r(1) - e_na)) - (g_l * (r(1) - e_l)) );
    result(2) = (alpha_m * (1 - r(2))) - (beta_m * r(2));   
    result(3) = (alpha_n * (1 - r(3))) - (beta_n * r(3));
end
