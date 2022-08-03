function gini_coeff = calculate_gini_coeff(matrix, row_wise_or_column_wise)
    % calculates gini coeff across row or column on a matrix
    % params
    %        matrix : a square matrix
    %        row_wise_or_column_wise = 0 or 1; 0 if row 1 if column

    
    % get the sum vector
    if row_wise_or_column_wise == 0 % if wanted to calculate row wise
        sum_vector = zeros(1, size(matrix,1));
        for r=1:size(matrix,1)
            sum_vector(1,r) = sum(matrix(r,:));
        end
    else
        sum_vector = zeros(1, size(matrix,2));
        for c=1:size(matrix,2)
            sum_vector(1,c) = sum(matrix(:,c));
        end    
    end
    
    % normalise the sum vector
    sum_vector_normalised = zeros(1, length(sum_vector));
    sum_all_elements = sum(sum_vector(:));
    for s=1:length(sum_vector)
        sum_vector_normalised(1,s) = sum_vector(1,s)/sum_all_elements;
    end

    sorted_sum_vector_normalised = sort(sum_vector_normalised);
    cummulative_sum_sorted_sum_vector_normalised = cumsum(sorted_sum_vector_normalised);
    x_axis_normalised = (1:length(sum_vector)) / length(sum_vector);
    area_under_curve = trapz(x_axis_normalised,cummulative_sum_sorted_sum_vector_normalised);
    gini_coeff = 1 - 2*area_under_curve;

end