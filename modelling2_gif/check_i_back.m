close all
iter_to_see = 1;
col = 1;
figure
    hold on
        for n=1:2
                x = squeeze(I_background_tensor(iter_to_see,col,n,:));
                x = x(1:100);
                plot(x)
            
        end
    hold off
grid