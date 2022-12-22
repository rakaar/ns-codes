%% 1
close all
weights = zeros(300,1);
ee = 100;
for i=1:300
    x = randi(100);
    if x <= 33
        weights(i) = ee;
    elseif x > 34 && x < 67
        weights(i) = ee/2;
    else
        weights(i) = ee/4;
    end
end

figure
    hist(weights)
grid

%% 2 
close all


weights = zeros(3000,1);
ee = 100;
for i=1:3000
    x = rand;
    if x <= 1/3
        weights(i) = exp(unifrnd(log(ee)-(log(2)/2),log(ee)+(log(2)/2)));
        
    elseif x > 1/3 && x < 2/3
        weights(i) = exp(unifrnd(log(ee/2)-(log(2)/2),log(ee/2)+(log(2)/2)));
    else
        weights(i) = exp(unifrnd(log(ee/4)-(log(2)/2),log(ee/4)+(log(2)/2)));
    end
end

figure
    hist(weights)
grid

figure
    hist(log(weights))
grid

figure
    hist(log(weights),100)
grid
