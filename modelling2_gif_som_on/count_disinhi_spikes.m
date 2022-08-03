images_path = "D:\disinhi_spikes\";
iter = 1;
token_length = 2*single_stimulus_duration + gap_duration + pre_token_silence + post_token_silence;
num_of_disinhi_spikes = zeros(n_columns, n_excitatory, n_tokens);
for c=1:n_columns
    for n=1:n_excitatory
        for tok=1:n_tokens
           time_at_end_of_a = (tok-1)*token_length + 61;
           threshold_at_end_of_a = theta_tensor(iter,c,n,time_at_end_of_a);

           time_end_of_b = (tok-1)*token_length + token_length;
           for t=time_at_end_of_a+5:time_end_of_b
                threshold_at_t = theta_tensor(iter,c,n,t);
                if spikes(iter,c,n,t) == 1
                    if threshold_at_t < threshold_at_end_of_a
                        num_of_disinhi_spikes(c,n,tok) = sum(spikes(iter,c,n,time_at_end_of_a+5:time_end_of_b));
                    end

                    break
                end
           end
        end
     end
end

for c=1:n_columns
    figure
%         mean_num_of_neurons_giving_disinhi = mean(squeeze(num_of_disinhi_spikes(c,:,:)), 1);
%         plot(mean_num_of_neurons_giving_disinhi);
    imagesc(squeeze(num_of_disinhi_spikes(c,:,:)))
        title(['num of disinhi spikes-col-',num2str(c)])
    grid
end