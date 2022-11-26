function x = limit_x(old_x)
    if old_x < 0
        x = 0;
    elseif old_x > 1
        x = 1;
    else
        x = old_x;
    end
end