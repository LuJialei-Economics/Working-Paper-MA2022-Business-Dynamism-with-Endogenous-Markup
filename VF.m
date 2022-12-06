% Value function iteration given values of w and lambda

function [V1_pri,V2_pri] = VF(w,lambda, z_grid, Zprob,V2M, V2D, V20,V1)

    global beta a theta delta num_state h; 

    % value function iteraiton based on the HJBs 

    iter = 0;
    tor = 1; % torlerance of error

    pi_M = w/a * (theta/(theta-1)*w./z_grid).^((1-theta))/theta; 
                                                    % Monopoly Profit
    pi_D = zeros(num_state,num_state);              % Duopoly Profit
    for i = 1:num_state
        for j = 1: num_state
            % market share
            s = ((1-theta)*z_grid(j) + theta * z_grid(i))/(z_grid(i) + z_grid(j));
            pi_D(i,j) = w/a * s^2/theta * (theta/(theta-s)*w/z_grid(i))^(1-theta);
        end
    end
    

    while tor > 10^(-10)
    
        iter = iter+1;
        V2 = Mixed_V2(z_grid,V2M, V2D, V20);
        
        V2M_pri = repmat(pi_M,1, num_state) + beta * (1-delta) * ...
            Zprob * (Zprob * V2');
        V2D_pri = pi_D + beta * (1-delta) * Zprob * (Zprob * V2');
        V20_pri = beta * (1-delta) * Zprob * (Zprob * V2');

        V1_pri = pi_M + beta * (1-delta-lambda*delta/(1-lambda))*Zprob*V1 +...
            (beta * lambda*delta/(1-lambda))*Zprob*(h'*V2)';

        V2_pri = Mixed_V2(z_grid,V2M_pri, V2D_pri, V20_pri);

        tor1 = max(abs(V1-V1_pri));
        tor2 = max(max(abs(V2-V2_pri)));
        tor = max(tor1,tor2);

        V1 = V1_pri;
        V2M = V2M_pri;
        V2D = V2D_pri;
        V20 = V20_pri;

    
    end

end