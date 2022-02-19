N_excitatory = 100;
N_inhibitory = 25;
Total_neurons = N_excitatory + N_inhibitory;
weight_matrix = zeros(Total_neurons, Total_neurons);

% generate matrix
Jee = 6;
Jie = 0.5;
Jei = -4;
Jii = -0.5;

for i=1:1:N_excitatory
    for j=1:1:N_excitatory
        weight_matrix(i,j) = Jee;
    end
end

for i=1:1:N_inhibitory
    for j=1:1:N_inhibitory
        weight_matrix(i,j) = Jii;
    end
end

for i=1:1:N_excitatory
    for j=1:1:N_inhibitory
        weight_matrix(i,j) = Jei;
    end
end


for i=1:1:N_inhibitory
    for j=1:1:N_excitatory
        weight_matrix(i,j) = Jie;
    end
end

sample_column = column(weight_matrix);

disp(sample_column.weights(10,10))
