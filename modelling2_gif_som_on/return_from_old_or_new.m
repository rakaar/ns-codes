function value = return_from_old_or_new(old_vector, new_vector, t)
% takes old and new tensors
% based on value of t,quantity returns the value from appropriate tensor
% params
% old_vector and new_vector are old and new tensors of any shape
% t is the time or array index, 
%       if its non-positive, it will fetch from old vector

    if t >= 1
        value = new_vector(t);
    else
        t_new = length(old_vector) + t;
        value = old_vector(t_new);
    end

end