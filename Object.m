% define the object function for Nelder-Mead

function Obj = Object(w, lambda, C_M, C_D, Zprob, z_grid,V2M, V2D, V20,V1)

    global num_state h;

    % transfer the root-finding problem to an optimization problem
    % define the object function as the L-2 norm
    [gM, GM] = GM_dist(Zprob, lambda);
    [V1_pri,V2_pri] = VF(w,lambda, z_grid, Zprob,V2M, V2D, V20,V1);
    
    
    EV1 = h'*V1_pri;  % expected value of entering as a monopolist
    EV2 = h'*V2_pri*gM; % expected value of entering as a duopolist

    Obj = (EV1 - C_M)^2 + (EV2 - C_D)^2;

end