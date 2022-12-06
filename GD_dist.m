% Calculate the steady-state productivity distribution of Duopolists

function [gD, GD, gJ] = GD_dist(Zprob, g)

    global delta num_state h;
    
    %------------------------------------------------------------------------ 
    % gD(z')*gD(z_')' = (1-delta) [(Zprob' * gD(z) )*(Zprob' * gD(z_))'] + 
    %                 (delta/2) [(Zprob' * gM(z)) * h(z_')'] +
    %                 (delta/2) [h(z') * (Zprob' * gM(z_))']
    %------------------------------------------------------------------------
    
    L = eye(num_state*num_state) - (1-delta) * kron(Zprob', Zprob');
    R1 = kron(h, Zprob')*g;
    R2 = kron(Zprob', h)*g;
    R = (R1 + R2);

    gJ = L\(delta/2*R);
    gJ = reshape(gJ, [num_state,num_state]);
    gD = sum(gJ, 1);

    GD = cumsum(gD);    
end