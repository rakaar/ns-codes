iter = 1;
for t=2:length(tspan)
    fprintf("\n t is %d \n",t)
    oldest_impact_time = t-20;
    if oldest_impact_time < 1
        oldest_impact_time = 1;
    end 

    for n1=1:num_network_neurons
        for n2=1:num_network_neurons
            c1 = floor(n1/n_total_neurons) + 1;
            c2 = floor(n2/n_total_neurons) + 1;
            
            if c1 == 6
                c1 = 5;
            end
            if c2 == 6
                c2 = 5;
            end

            n1_index_in_column = mod(n1,n_total_neurons);
            if n1_index_in_column == 0
                n1_index_in_column = 25;
            end
    
            n2_index_in_column = mod(n2,n_total_neurons);
            if n2_index_in_column == 0
                n2_index_in_column = 25;
            end 

%             if n1_index_in_column > n_excitatory || n2_index_in_column > n_excitatory
%                 continue; % plastic connections only between excitatory columns
%             end

            % no learning rule checked
            if spikes(iter,c1,n1_index_in_column,t) == 0 && spikes(iter,c2,n2_index_in_column,t) == 0 
                if network_weight_matrix(iter,t,n1,n2) == network_weight_matrix(iter,t-1,n1,n2)
                    disp("okkkk")
                else
                    disp("XXXXXXXXXXXXXXXXX PROBLEM no learning XXXXXXXXXXXXXXXXX")
                    return       
                end
            end

            % LTD check
            if network_weight_matrix(iter,t,n1,n2) < network_weight_matrix(iter,t-1,n1,n2) 
                % 1. only pre has a spike at t, post has no spike at t,
                % 2.  post has spike in past
                pre_has_spike = (spikes(iter,c1,n1_index_in_column,t) == 1) && (spikes(iter,c2,n2_index_in_column,t) == 0);

                post_has_spike_in_past = 0;
                for ttt=t-1:-1:oldest_impact_time
                    if spikes(iter,c2,n2_index_in_column,ttt) == 1
                        post_has_spike_in_past = 1;
                        break;
                    end
                end

                if post_has_spike_in_past && pre_has_spike
                    disp("✅ LTD")
                else
                    fprintf("\n post_has_spike_in_past- %d, pre_has_spike - %d \n", post_has_spike_in_past,pre_has_spike)
                    disp("XXXXXXXXXXXXXXXXX PROBLEM XXXXXXXXXXXXXXXXXXXXXXXXXX")
                    return
                end
                
            end % end of LTD check

            % LTP check
            if network_weight_matrix(iter,t,n1,n2) > network_weight_matrix(iter,t-1,n1,n2) 
                % post has a spike at t, pre has no spike at t
                % pre has spike in past

                post_has_spike = (spikes(iter,c2,n2_index_in_column,t) == 1) && (spikes(iter,c1,n1_index_in_column,t) == 0);

                pre_has_spike_in_past = 0;
                for ttt=t-1:-1:oldest_impact_time
                    if spikes(iter,c1,n1_index_in_column,ttt) == 1
                        pre_has_spike_in_past = 1;
                        break;
                    end
                end

                if post_has_spike && pre_has_spike_in_past
                    disp("✅ LTP")
                else
                    fprintf("\n pre_has_spike_in_past  - %d, post_has_spike- %d \n", pre_has_spike_in_past,post_has_spike)
                    disp("xxxxxxxxxxxxxxxxxxxPROBLEMxxxxxxxxxxxxxxxxx")
                    return
                end
            end % end of LTP check
        end
    end
end

%% checking for a particular connection
close all
% 3 -> 5 LTD y no?, 3 spike 5 no spike and 5 past spike

col = 3;
pre_index_in_col = 3;
post_index_in_col = 5;

pre_index_in_network = (col-1)*n_total_neurons + pre_index_in_col;
post_index_in_network = (col-1)*n_total_neurons + post_index_in_col;

iter = 1;
for t=2:length(tspan)
    oldest_impact_time = t-20;
    if oldest_impact_time < 1
        oldest_impact_time = 1;
    end

    if spikes(iter,col,pre_index_in_col,t) == 1 && spikes(iter,col,post_index_in_col,t) == 0 
        post_has_spike_in_past = 0;
        for ttt=t-1:-1:oldest_impact_time
            if spikes(iter,col,post_index_in_col,ttt) == 1
                post_has_spike_in_past = 1;
                break;
            end
        end

        if post_has_spike_in_past == 1
            disp(" xxxxxxxxxxxxxxx PROBLEM - LTD has to be there xxxxxxxxxxxxxxx")
            fprintf("\n past weight %f , now weight %f \n", network_weight_matrix(iter,t-1,pre_index_in_network,post_index_in_network),network_weight_matrix(iter,t,pre_index_in_network,post_index_in_network))
            return
        else
            disp("okkkk")
        end
    end   
end