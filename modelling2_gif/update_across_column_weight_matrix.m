function updated_weight_matrix = update_across_column_weight_matrix(old_weight_matrix, spikes_pre_syn_column, spikes_post_syn_column, t, n_excitatory)
    % old_weight_matrix- old_weight_matrix of a pair of columns
    % spikes_pre_syn_column - spike trains of pre syn column
    % spikes_post_syn_column - spike trains of post syn column
    % t - time

    % returns updated_weight_matrix
            either_LTP_or_LTD_occured = zeros(n_excitatory, n_excitatory);
            for N=1:n_excitatory
                % N -> postyn : LTD
                for postsyn_neuron=1:n_excitatory
                    if spikes(iter,col_stdp,N,i) == 1 % if there is a spike
                        found_spike_in_window_LTD = 0;

                        for postsyn_spike_time=i-1:-1:i-20
                            if postsyn_spike_time >= 1 && spikes(iter,col_stdp,postsyn_neuron,postsyn_spike_time) == 1 && spikes(iter,col_stdp,postsyn_neuron,i) == 0 
                                exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron) = exc_to_exc_weight_matrix(iter,col_stdp,i-1,N,postsyn_neuron)*(1 - Amp_weak*exp(-abs(i-postsyn_spike_time)/tau_weak));
                                % clipping weights
                                if exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron) < minimum_weight_exc_to_exc
                                    exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron) = minimum_weight_exc_to_exc;
                                end
                                if exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron) > maximum_weight_exc_to_exc
                                    exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron) = maximum_weight_exc_to_exc;
                                end
                                num_of_LTDs(iter,col_stdp,i) = num_of_LTDs(iter,col_stdp,i) + 1;
                                
                                either_LTP_or_LTD_occured(N,postsyn_neuron) = 1;
                                found_spike_in_window_LTD = 1;
                                break;
                            end
                        end

                        if found_spike_in_window_LTD == 0 && either_LTP_or_LTD_occured(N,postsyn_neuron) == 0
                            exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron) = exc_to_exc_weight_matrix(iter,col_stdp,i-1,N,postsyn_neuron);

                        end
                    else % if there is no spike

                        if either_LTP_or_LTD_occured(N,postsyn_neuron) == 0
                               exc_to_exc_weight_matrix(iter,col_stdp,i,N,postsyn_neuron) = exc_to_exc_weight_matrix(iter,col_stdp,i-1,N,postsyn_neuron);

                        end
                        
                    end
                end


                % presyn -> N : LTP
                for presyn_neuron=1:n_excitatory
                    if spikes(iter,col_stdp,N,i) == 1 % if there is a spike
                        found_spike_in_window_LTP = 0;

                        for presyn_spike_time=i-1:-1:i-20
                            if presyn_spike_time >= 1 && spikes(iter,col_stdp,presyn_neuron,presyn_spike_time) == 1 && spikes(iter,col_stdp,presyn_neuron,i) == 0
                                    exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N) = exc_to_exc_weight_matrix(iter,col_stdp,i-1,presyn_neuron,N)*(1 + Amp_strength*exp(-abs(i-presyn_spike_time)/tau_strength));
                                    % clipping weights 
                                    if exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N) < minimum_weight_exc_to_exc
                                        exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N) = minimum_weight_exc_to_exc;
                                    end
                                    if exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N) > maximum_weight_exc_to_exc
                                        exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N) = maximum_weight_exc_to_exc;
                                    end

                                    num_of_LTPs(iter,col_stdp,i) = num_of_LTPs(iter,col_stdp,i) + 1; 
                                either_LTP_or_LTD_occured(presyn_neuron,N) = 1;
                                found_spike_in_window_LTP = 1;
                                break
                            end
                        end

                        if found_spike_in_window_LTP == 0 && either_LTP_or_LTD_occured(presyn_neuron,N) == 0
                            exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N) = exc_to_exc_weight_matrix(iter,col_stdp,i-1,presyn_neuron,N);

                        end
                    else % if there is no spike
                        if either_LTP_or_LTD_occured(presyn_neuron,N) == 0
                            exc_to_exc_weight_matrix(iter,col_stdp,i,presyn_neuron,N) = exc_to_exc_weight_matrix(iter,col_stdp,i-1,presyn_neuron,N);

                        end
                     end
                end

                

            end % end of for all excitatory neurons

end