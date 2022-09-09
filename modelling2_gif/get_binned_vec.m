function binned_vec = get_binned_vec(vec, bin_size)
    vec_length = length(vec);
    vec_reshaped = reshape(vec, bin_size, vec_length/bin_size);
    binned_vec = mean(vec_reshaped, 1);
end