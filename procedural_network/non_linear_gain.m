function gain = non_linear_gain(w)
  %   E_max = 300; % original in paper
   E_max = 0.3; 
    if w < 0
        gain = 0;
    elseif (w >= 0) & (w < E_max)
        gain = w;
    else 
        gain = E_max;
    end
end