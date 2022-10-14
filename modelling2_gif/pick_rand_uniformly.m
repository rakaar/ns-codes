function rand_num = pick_rand_uniformly(num, error)
  % returns a random number btn num-err% and num+err%
  % num = number about which fluctuation is allowed
  % error = percentage of fluctuation

    a = num*(1 - (error/100));
    b = num*(1 + (error/100));
    
    rand_num = a + (b-a)*rand;
end