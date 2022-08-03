function new_tensor = remove_diagonal_elements(old_tensor)
    % given old_tensor of size time x 20 x 20
    % returns new_tensor of size time x 20 x 19
    
    n_exc = size(old_tensor, 2);
    time_period = size(old_tensor,1);
    new_tensor = zeros(time_period,n_exc,n_exc-1);
    for t=1:time_period
        matrix_at_t = zeros(n_exc, n_exc-1);
        
        for n=1:n_exc % rows
          
            if n == 1
                row_without_dia = reshape(old_tensor(t,n,2:n_exc),  1,n_exc-1); 
            elseif n == n_exc
                row_without_dia = reshape(old_tensor(t,n,1:n_exc-1),  1,n_exc-1);
            else
                pre_dia = reshape(old_tensor(t,n,1:n-1),  1,n-1);
                post_dia = reshape(old_tensor(t,n,n+1:n_exc), 1,n_exc-n);
                row_without_dia = [pre_dia post_dia];
            end

            
            new_tensor(t,n,:) = reshape(row_without_dia, 1,1,n_exc-1);
        end % end of rows
   end
end