function result_vec = non_linear_gain_vector(vec)
    result_vec = zeros(1, length(vec));
    
    for i=1:length(vec)
        result_vec(1, i) = non_linear_gain(vec(1,i));
    end
end