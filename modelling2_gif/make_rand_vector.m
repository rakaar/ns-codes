function rand_vec = make_rand_vector(num,err,vec_size)
    % returns a random vec of number of size-vec_size btn num-err% and num+err%
   % num = number about which fluctuation is allowed
   % error = percentage of fluctuation
    a = num*(1 - (err/100));
    b = num*(1 + (err/100));
    
    rand_vec = a + (b-a)*rand(vec_size);
end