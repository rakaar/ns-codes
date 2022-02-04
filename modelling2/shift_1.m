function shifted_arr = shift_1(arr)
    len = length(arr);
    new_arr = zeros(1, len);
    for i=2:len
        new_arr(1,i) = arr(1, i-1); 
    end

    shifted_arr = new_arr;
end
