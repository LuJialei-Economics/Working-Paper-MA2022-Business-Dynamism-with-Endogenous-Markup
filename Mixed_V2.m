% Mixed Duopolist's Value Function. 
function V2 = Mixed_V2(z_grid,V2M, V2D, V20)

    global theta num_state;

    %------------------------------------------------------------------------
    %
    % Goal: the RHS of value functional equations have a weighted avergage
    % values of V2. This function aims to combine different parts 
    % of V2M, V2D, and V20 to prepare the weighted V2
    % The base of mixture is the value of z
    %
    %------------------------------------------------------------------------

    V2 = zeros(num_state, num_state);

    for i = 1:num_state
        upper = z_grid(i)*(theta)/(theta-1);
        lower = z_grid(i)*(theta-1)/theta;

        % if lower > z_grid(end), V2(i,:) = V2M(i,:)
        % if upper > z_grid(end) & z_grid(1) < lower < z_grid(end), 
        %                                       V2(i,:) = V2M(i,:):V2D(i,:)
        % if upper > z_grid(end) & lower < z_grid(1), V2(i,:) = V2D(i,:)
        % if upper < z_grid(1), V2(i,:) = V20(i,:)
        % if upper < z_grid(end) & lower < z_grid(1), 
        %                                       V2(i,:) = V2D(i,:):V20(i,:)
        % if upper < z_grid(end) & lower > z_grid(1), 
        %                              V2(i,:) = V2M(i,:):V2D(i,:):V20(i,:)
        U = upper>z_grid;
        L = lower>z_grid;
        V2(i,:) = ((1-L).*U)'.*V2D(i,:) + (L.*U)'.*V2M(i,:) + ...
            ((1-L).*(1-U))'.*V20(i,:);
    end
    
end