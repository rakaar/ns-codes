function [mean_for_iters, mean_for_neurons] = get_mean(var_tensor, n_iters, n_neurons, time_length,column_index)
    % params
    %   var_tensor: a 4 dim tensor,
    %              num_iters x num_columns x num_neurons x time
    %  n_iters: num of iterations
    %  n_neurons: num of neurons
    %  time_length: time period
    %   column_index: column in the model

    % returns
    % mean_for_iters: a 2D matrix
    %                 n_neurons x time_length
    %  mean_for_neurons: a vector
    %                    1 x time_length

    mean_for_iters = zeros(n_neurons, time_length);
    mean_for_neurons = zeros(1, time_length);
  
    for t=1:time_length
        for n=1:n_neurons
            mean_for_iters(n,t) = sum(var_tensor(:,column_index,n,t))/n_iters;
        end
    end

    for t=1:time_length
        mean_for_neurons(1,t) = sum(mean_for_iters(:,t))/n_neurons;
    end

end