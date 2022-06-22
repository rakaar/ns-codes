close all;clear all

batches = 200;
cols = 5;
batch_data_path = "D:\7_multi_col_big_clip_range";

for b=1:batches
    fprintf("\n \n bbbbbbbbbbbbbbbbath % d bbbbbbbbbbbbbbbbb\n ", b);
    batch_file_name = batch_data_path + "\batch_" + num2str(b) + ".mat";
        
    batch_spikes_all = load(batch_file_name, "spikes").spikes;
    t_simulate = load(batch_file_name, "t_simulate").t_simulate;
    batch_weight_matrix_all = load(batch_file_name, "exc_to_exc_weight_matrix").exc_to_exc_weight_matrix;
        
     for c=1:cols
        batch_spikes = squeeze(batch_spikes_all(1,c,:,:));
        batch_weight_matrix = squeeze(batch_weight_matrix_all(1,c,:,:,:));
      
        for t=2:t_simulate
            fprintf("\n col %d, t  %d \n", c,t);
             oldest_impact_time = t-20;
             if oldest_impact_time < 1
                oldest_impact_time = 1;
             end 
            
             for pre=1:20
                for post=1:20
                    % LTD check
                    if batch_weight_matrix(t,pre,post) < batch_weight_matrix(t-1,pre,post) 
                        % 1. pre has a have spike
                        % 2. post has a spike in past 20

                        spike_in_pre = (batch_spikes(pre,t) == 1 && batch_spikes(post,t) == 0);
                       
                        
                        spike_in_post_past = 0;
                        for ttt=t-1:-1:oldest_impact_time
                            if batch_spikes(post,ttt) == 1
                                spike_in_post_past = 1;
                                break
                            end
                        end
    
                        if (spike_in_post_past && spike_in_pre)
                            disp("✅ LTD")
                        else
                            disp("xxxxxxxxxxxxxxxxxxxPROBLEMxxxxxxxxxxxxxxxxx")
                            return
                        end 

                    end
                                
                    % LTP check
                    if batch_weight_matrix(t,pre,post) > batch_weight_matrix(t-1,pre,post) 
                        % 1. post doesnt have spike
                        % 2. pre has a spike in past 20

                        spike_in_post = (batch_spikes(post,t) == 1 && batch_spikes(pre,t) == 0);

                        spike_in_pre_past = 0;
                        for ttt=t-1:-1:oldest_impact_time
                            if batch_spikes(pre,ttt) == 1
                                spike_in_pre_past = 1;
                                break;
                            end
                        end
    
                        if (spike_in_pre_past && spike_in_post)
                            disp("✅ LTP")
                        else
                            disp("xxxxxxxxxxxxxxxxxxxPROBLEMxxxxxxxxxxxxxxxxx")
                            return
                        end 

                    end

                    
                end
            end
        end
     end
end