function loop_and_see
    standard_a=0.03; standard_b=0.25; standard_c=-52;  standard_d=0;
    
    disp('aaaaaaaaaaaaaaaaaaaaaaaaaaaaa')
    first_spike_arr = [];
    num_of_spikes_arr = [];
    for i=-0.10:0.01:0.1
        voltage_array = run_izhikevich(standard_a+i, standard_b, standard_c, standard_d);
        first_spike_arr = [first_spike_arr, get_first_spike_time(voltage_array, 1)];
        num_of_spikes_arr = [num_of_spikes_arr, get_num_of_spikes(voltage_array)];
    end
    figure(10)
        scatter(-0.10:0.01:0.1, first_spike_arr, 'red','filled');
        title('a vs first spike')
    grid
    figure(11)
        scatter(-0.10:0.01:0.1, num_of_spikes_arr,'red', 'filled');
        title('a vs num of spikes')
    grid

    disp('bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb')
    first_spike_arr = [];
    num_of_spikes_arr = [];
    for i=-0.10:0.01:0.1
        voltage_array = run_izhikevich(standard_a, standard_b+i, standard_c, standard_d);
        first_spike_arr = [first_spike_arr, get_first_spike_time(voltage_array, 1)];
        num_of_spikes_arr = [num_of_spikes_arr, get_num_of_spikes(voltage_array)];
    end
    figure(20)
        scatter(-0.10:0.01:0.1, first_spike_arr, 'red','filled');
        title('b vs first spike')
    grid
    figure(21)
        scatter(-0.10:0.01:0.1, num_of_spikes_arr,'red', 'filled');
        title('b vs num of spikes')
    grid

    disp('ccccccccccccccccccccccccccccccc')
    first_spike_arr = [];
    num_of_spikes_arr = [];
    for i=-10:1:10
        voltage_array = run_izhikevich(standard_a, standard_b, standard_c+i, standard_d);
        first_spike_arr = [first_spike_arr, get_first_spike_time(voltage_array, 1)];
        num_of_spikes_arr = [num_of_spikes_arr, get_num_of_spikes(voltage_array)];
    end
    figure(30)
        scatter(-10:1:10, first_spike_arr, 'red','filled');
        title('c vs first spike')
    grid
    figure(31)
        scatter(-10:1:10, num_of_spikes_arr,'red', 'filled');
        title('c vs num of spikes')
    grid
    
    disp('ddddddddddddddddddddddddddddddd')
    first_spike_arr = [];
    num_of_spikes_arr = [];
    for i=-10:1:10
        voltage_array = run_izhikevich(standard_a, standard_b, standard_c, standard_d+i);
        first_spike_arr = [first_spike_arr, get_first_spike_time(voltage_array, 1)];
        num_of_spikes_arr = [num_of_spikes_arr, get_num_of_spikes(voltage_array)];
    end
    figure(40)
        scatter(-10:1:10, first_spike_arr, 'red','filled');
        title('d vs first spike')
    grid
    figure(41)
        scatter(-10:1:10, num_of_spikes_arr,'red', 'filled');
        title('d vs num of spikes')
    grid
end
function voltage_values = run_izhikevich(a,b,c,d)
    % a=0.03; b=0.25; c=-52;  d=0;
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
    
    fprintf("\n a=%f b=%f c=%f d=%f \n", a, b, c, d);
    num_of_spikes = get_num_of_spikes(VV);
    first_spike = get_first_spike_time(VV,1);
    first_spike_after_negative_i_cuts = get_first_spike_time(VV, 25000);
    fprintf("%d spikes \n", num_of_spikes);
    fprintf(" %f, %f \n", first_spike, first_spike_after_negative_i_cuts);

    voltage_values = VV;
    % plot(tspan,VV,[0 T1 T1 (T1+5) (T1+5) max(tspan)],-85+[0 0 -5 -5 0 0]);
    % axis([0 max(tspan) -90 30])
    % title([' value ', num2str(i)]);
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

function return_value = get_first_spike_time(voltage_arr, from_time)
    len = length(voltage_arr);
    return_value = 0;
    spike_time = 0;
    for i=from_time:len
        if voltage_arr(i) >= 30
            spike_time = i/1000;
            break;
        end
    end
    return_value = spike_time;
end


