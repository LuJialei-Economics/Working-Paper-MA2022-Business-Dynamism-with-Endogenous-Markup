% Tauchen's Method
function [Z,Zprob]= Tauchen(param,varargin)

    % TAUCHEN 
    % 
    % Tauchen method to discretize AR(1) process based on the following
    % formula:
    %
    % y_t = \rho * y_{t-1} + epsilon, epsilon ~ Normal(mu,sig^2)
    % num_state: number of states
    % m: numbers of standard errors to simulate
    % --------------------------------------------------------------------
    
    m = 5;
    Nz = param.Nz;
    mu = param.mu;

    par = inputParser;
    addOptional(par,'sig',param.sig);
    addOptional(par,'rho',param.rho);
    parse(par,varargin{:});

    rho = par.Results.rho;
    sig = par.Results.sig;


    % initialize the initial values and transition matrix
    Z = zeros(Nz,1); 
                %column vector of state values in state space, N states
    Zprob = zeros(Nz,Nz); %transition matrix
    
    % set the boundary points
    Z(1) = mu/(1-rho) - m*sqrt(sig^2/(1-rho^2));
    Z(Nz) = mu/(1-rho) + m*sqrt(sig^2/(1-rho^2));

    % set other points
    zstep=(Z(Nz)-Z(1))/(Nz-1);
    for i=2:Nz-1
        Z(i)=Z(i-1)+zstep;
    end

    % set the transition matrix
    for j = 1:Nz
        for k = 1:Nz
            if k == 1
                Zprob(j,k) = normcdf((Z(1)-mu-rho*Z(j)...
                    +zstep/2)/sig);
            elseif k == Nz
                Zprob(j,k) = 1 - normcdf((Z(Nz) ...
                    -mu-rho*Z(j)-zstep/2)/sig);
            else
                Zprob(j,k) = normcdf((Z(k)-mu-rho*Z(j) ...
                    +zstep/2)/sig)... 
                    -normcdf((Z(k)-mu-rho*Z(j) ...
                    -zstep/2)/sig);
            end
        end
    end
end

