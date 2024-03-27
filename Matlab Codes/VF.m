function [V1_pri,V2_pri,pi1,pi2] = VF(SMatrix, glob, param, varargin)

    %VF
    %
    % Value function iteration given values of w and lambda
    %

    par = inputParser;
    addOptional(par,'lambda',glob.lambda);
    addOptional(par,'w',glob.w);
    parse(par, varargin{:});

    lambda = par.Results.lambda;
    w = par.Results.w;

    beta = param.beta;
    a = param.a;
    theta = param.theta;
    gamma = param.gamma;
    Nz = param.Nz;
    h = glob.h;
    z_grid = glob.z_grid;
    Zprob = glob.Zprob;
    delta = glob.delta; 
    s = SMatrix;

    iter = 0;
    tor = 1; % torlerance of error
    maxiter = 500;

    % initialize the value functions
    V1 = zeros(Nz, 1);
    V2 = zeros(Nz, Nz);

    % profits
    xi = gamma + theta - gamma*theta;
    pi1 = (theta/(theta-1)/gamma - 1)*(theta/(theta-1)/gamma)^(-theta/xi)*...
        w^((xi-theta+1)/xi)*a^(-1/xi).*z_grid.^((theta-1)/xi);    
                                                    % Monopoly Profit
    zz_grid = kron(z_grid,ones(1,Nz));                                                 
    pi2 = s.^(1/gamma) .* (theta./(theta-s)./gamma - 1).*...
        (theta./(theta-s)./gamma).^(-theta/xi).*w^((xi-theta+1)/xi).*...
        a^(-1/xi).*zz_grid.^((theta-1)/xi);                                 
                                                    % Duopoly Profit
    
    % parameters
    parT = (1-lambda)*(1-delta)+lambda*delta;
    parA = ((1-lambda)*(1-delta)^2)/parT;
    parB = ((1-delta)*lambda*delta)/parT;
    parC = 0.5*((1-lambda)*(1-delta)*delta)/parT;
    parD = 0.5*(lambda*delta^2)/parT;

    % iteration
    while tor > 10^(-10)
    
        iter = iter+1;
        
        V2_pri = pi2 + beta*((1-delta)*(Zprob*(Zprob*V2)')' +...
            parC*Zprob*V1 + parD*(kron(h,ones(1,Nz)))*(Zprob*V2)');
        V1_pri = pi1 + beta*parA*Zprob*V1 + beta*parB *Zprob*V2*h;

        tor1 = max(abs(V1-V1_pri));
        tor2 = max(max(abs(V2-V2_pri)));
        tor = max(tor1,tor2);

        V1 = V1_pri;
        V2 = V2_pri;

        disp(iter);
        disp(tor);

        if iter > maxiter
            disp('Reach the Maximum Iteration')
            break
        end

    end

end


