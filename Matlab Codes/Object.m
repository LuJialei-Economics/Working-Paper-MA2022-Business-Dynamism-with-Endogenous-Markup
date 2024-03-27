function Obj = Object(x, C_M, C_D, glob, param)

    %OBJECT
    %
    % Define the object function for minimize, when Obj = 0, the results
    % are the corresponding wage (w) and proportion of duopolist markets
    % (lambda).
    %

    Zprob = glob.Zprob;
    zgrid = glob.z_grid;
    Nz = param.Nz;
    h = glob.h;

    lambda = x(1);
    w = x(2);

    load('SMatrix.mat','SMatrix');
    % transfer the root-finding problem to an optimization problem
    % define the object function as the L-2 norm
    [~, ~, gM, ~] = SteadyState_Dist(glob,param, 'lambda',lambda);
    [V1_pri,V2_pri] = VF(SMatrix, glob, param, 'lambda',lambda, 'w', w);
    
    
    EV1 = h'*V1_pri;  % expected value of entering as a monopolist
    EV2 = h'*V2_pri*gM; % expected value of entering as a duopolist

    Obj = (EV1 - C_M)^2 + (EV2 - C_D)^2;

end