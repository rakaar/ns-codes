iter = 1;
 for cc=1:num_connected_pairs
    c1 = connected_columns_arr(cc,1);
    c2 = connected_columns_arr(cc,2);
    
    for t=2:length(tspan)
        fprintf("\n c1 c2 %d %d, t  %d \n", c1,c2,t);
             oldest_impact_time = t-20;
             if oldest_impact_time < 1
                oldest_impact_time = 1;
             end 


             for pre=1:20
                for post=1:20
                    
                    % LTD check
                    if across_columns_exc_to_exc_weight_matrix(iter,c1,c2,t,pre,post) < across_columns_exc_to_exc_weight_matrix(iter,c1,c2,t-1,pre,post)
                         % 1. pre has a have spike
                         % 2. post has a spike in past 20

                         pre_has_spike = (spikes(iter,c1,pre,t) == 1);
                          
                         post_has_past_spike = 0;
                         for i=t-1:-1:oldest_impact_time
                            if spikes(iter,c2,post,i) == 1
                                post_has_past_spike = 1;
                                break
                            end
                         end

                         if pre_has_spike && post_has_past_spike
                            disp("✅ LTD")
                         else
                            disp("xxxxxxxxxxxxxxxxxxxPROBLEM LTD xxxxxxxxxxxxxxxxx")
                            return
                         end
                    end



                    % LTP check
                    if across_columns_exc_to_exc_weight_matrix(iter,c1,c2,t,pre,post) > across_columns_exc_to_exc_weight_matrix(iter,c1,c2,t-1,pre,post)
                         % 1. post has a spike
                         % 2. pre has a spike in past 20

                         post_has_spike = (spikes(iter,c2,post,t) == 1);
                          
                         pre_has_past_spike = 0;
                         for i=t-1:-1:oldest_impact_time
                            if spikes(iter,c1,pre,i) == 1
                                pre_has_past_spike = 1;
                                break
                            end
                         end

                         if post_has_spike && pre_has_past_spike
                            disp("✅ LTP")
                         else
                            disp("xxxxxxxxxxxxxxxxxxxPROBLEM- LTP xxxxxxxxxxxxxxxxx")
                            return
                         end
                    end

                    
                end
             end
    end
 
end