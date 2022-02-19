function project
    [t r] = ode15s(@fitznagumo, [0 200], [-2.8; -1.8]);
    T1=20;
    figure(10)
        plot(t, r(:,1),[0 T1 T1 (T1+5) (T1+5) 200], -5 + [0 0 -0.5 -0.5 0 0]);
        title("fitz nagumo- v vs t");
    grid
end

function result = fitznagumo(t,r)
    a = 0.7; b = 0.8; c = 12.5;

    if t < 25 & t > 20
        i = -0.5;
    else
        i = 0;
    end

    result = zeros(2,1);
    v = r(1); w = r(2);
    result(1) = v - ((v^3)/3) - w + i;
    result(2) = (1/c)*(v + a - (b*w));
end