close all
iter_to_see = 2;
col = 1;

        for n=1:n_excitatory
            figure
                x = squeeze(I_background_tensor(iter_to_see,col,n,:));
                x = x(1:500);
                plot(x)

            grid
        end

