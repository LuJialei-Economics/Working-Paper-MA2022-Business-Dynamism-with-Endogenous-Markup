% Tauchen's Method
function [Z,Zprob]= tauchen(mu,sig)

    global rho num_state;

    %------------------------------------------------------------------------
    % y_t = \rho * y_{t-1} + epsilon, epsilon ~ Normal(mu,sig^2)
    % num_state: number of states
    % m: numbers of standard errors to simulate
    %------------------------------------------------------------------------
    
    m = 4;

    % initialize the initial values and transition matrix
    Z = zeros(num_state,1); %column vector of state values in state space, N states
    Zprob = zeros(num_state,num_state); %transition matrix
    
    % set the boundary points
    Z(1) = mu/(1-rho) - m*sqrt(sig^2/(1-rho^2));
    Z(num_state) = mu/(1-rho) + m*sqrt(sig^2/(1-rho^2));

    % set other points
    zstep=(Z(num_state)-Z(1))/(num_state-1);
    for i=2:num_state-1
        Z(i)=Z(i-1)+zstep;
    end

    % set the transition matrix
    for j = 1:num_state
        for k = 1:num_state
            if k == 1
                Zprob(j,k) = normcdf((Z(1)-mu-rho*Z(j)+zstep/2)/sig);
            elseif k == num_state
                Zprob(j,k) = 1 - normcdf((Z(num_state)-mu-rho*Z(j)-zstep/2)/sig);
            else
                Zprob(j,k) = normcdf((Z(k)-mu-rho*Z(j)+zstep/2)/sig)-...
                          normcdf((Z(k)-mu-rho*Z(j)-zstep/2)/sig);
            end
        end
    end
end