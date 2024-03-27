% Calculate the steady-state productivity distributions

function [gD, GD, gM, GM] = SteadyState_Dist(glob,param,varargin)

    %STEADYSTATE_DIST
    %
    % Calculate the steady-state distribution given the transition
    % relationsip below:
    %
    % gM(z')*gM(z_')' = delta * h + 
    %       (1-lambda)*(1-delta)^2)/((1-lambda)*(1-delta)+lambda*delta)
    %                 *(Zprob' * gM(z)) +
    %       ((1-delta)*lambda*delta)/((1-lambda)*(1-delta)+lambda*delta)
    %                 * (Zprob' * gD(z_)) 
    % gD(z')*gD(z_')' = (1-delta) [(Zprob' * gD(z) )*(Zprob' * gD(z_))'] + 
    %          (1/2)((1-lambda)*(1-delta)*delta)/((1-lambda)
    %           *(1-delta)+lambda*delta)*[(Zprob' * gM(z)) * h(z_')'] +
    %          (1/2)((1-lambda)*(1-delta)*delta)/((1-lambda)
    %           *(1-delta)+lambda*delta)*[h(z') * (Zprob' * gM(z_))'] + 
    %          (1/2)(lambda*delta^2)/((1-lambda)*(1-delta)+lambda*delta)
    %           *[(Zprob' * gD(z)) * h(z_')'] +
    %          (1/2)(lambda*delta^2)/((1-lambda)*(1-delta)+lambda*delta)
    %           *[h(z') * (Zprob' * gD(z_))']
    %
    
    delta = glob.delta;
    h = glob.h;
    Nz = param.Nz;

    par = inputParser;
    addOptional(par,'Zprob',glob.Zprob);
    addOptional(par,'lambda',glob.lambda);
    parse(par,varargin{:});

    Zprob = par.Results.Zprob;
    lambda = par.Results.lambda;

    % initilization
    gD = ones(Nz,1)/Nz;

    % parameters
    parA = ((1-lambda)*(1-delta)^2)/...
        ((1-lambda)*(1-delta)+lambda*delta);
    parB = ((1-delta)*lambda*delta)/...
        ((1-lambda)*(1-delta)+lambda*delta);
    parC = (1/2)*((1-lambda)*(1-delta)*delta)/...
        ((1-lambda)*(1-delta)+lambda*delta);
    parD = (1/2)*(lambda*delta^2)/...
        ((1-lambda)*(1-delta)+lambda*delta);
    
    % iteration
    err = 1;
    i=0;
    while err > 3*1e-5
        i = i + 1;
        gM = (eye(Nz)-parA * Zprob')\...
            (delta * h + parB * Zprob'*gD);
        gJ = (1-delta)*Zprob'*(gD*gD')*Zprob + parC *...
            (Zprob'*gM*h') + parC * (h*gM'*Zprob) +...
            parD * (Zprob'*gD*h') + parD * (h*gD'*Zprob);
        gD_tilde = sum(gJ,2);
        err = max(abs(gD_tilde - gD));
        gD = gD_tilde;
    end

    GD = cumsum(gD);
    GM = cumsum(gM);

end


