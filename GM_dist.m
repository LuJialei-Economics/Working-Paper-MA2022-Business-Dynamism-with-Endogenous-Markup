% Calculate the steady-state productivity distribution of Monopolists

function [gM, GM] = GM_dist(Zprob, lambda)

    global delta num_state h;
    
    %------------------------------------------------------------------------ 
    % g_M(z') = delta/(1-lambda)h(z') 
    %                   + [1-delta-lambda*delta/(1-lambda)] Zprob(z',z)g_M(z)
    %------------------------------------------------------------------------
    
    B = eye(num_state) - (1-delta-lambda*delta/(1-lambda))*Zprob';
    gM = B\((delta)/(1-lambda)* h);

    GM = cumsum(gM);    
end