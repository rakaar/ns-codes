function param_tune
    a=0.03; b=0.25; c=-52;  d=0;
    V=-64;  u=b*V;
    VV=[];  uu=[];
    tau = 0.001;  tspan = 0:tau:100;
    T1=20;
    for t=tspan
        if (t>T1) & (t < T1+5) 
            I=-15;
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
    
    fprintf("a=%f b=%f c=%f d=%f \n", a, b, c, d);
    num_of_spikes = get_num_of_spikes(VV);
    first_spike = get_first_spike_time(VV,1);
    first_spike_after_negative_i_cuts = get_first_spike_time(VV, 25000);
    fprintf("%d spikes \n", num_of_spikes);
    fprintf("%first spike-> from 0: %f, from 25ms: %f \n", first_spike, first_spike_after_negative_i_cuts);


    plot(tspan,VV,[0 T1 T1 (T1+5) (T1+5) max(tspan)],-85+[0 0 -5 -5 0 0]);
    axis([0 max(tspan) -90 30])
    title([' value ', num2str(i)]);
end

function num_of_spikes = get_num_of_spikes(voltage_arr)
    len = length(voltage_arr);
    num_of_spikes = 0;

    for i=1:len
        if voltage_arr(i) >= 30
            num_of_spikes = num_of_spikes + 1;
        end
    end
end

function spike_time = get_first_spike_time(voltage_arr, from_time)
    len = length(voltage_arr);
    for i=from_time:len
        if voltage_arr(i) >= 30
            spike_time = i/1000;
            break;
        end
    end
end


