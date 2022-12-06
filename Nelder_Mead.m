% Nelder-Mead Algorithm to solve equilibrium w and lambda
function [w_e, lambda_e] = Nelder_Mead(W, Lambda, C_M, C_D, Zprob, z_grid,V2M, V2D, V20,V1)

    tic

    % given the original w and lambda, calculate the corresponding Obj's values
    diff = zeros(3,1);
    for i = 1:3
        
        diff(i) = Object(W(i),Lambda(i), C_M, C_D, Zprob, z_grid,V2M, V2D, V20,V1);
    end

    epsi = 1;

    while epsi > 10^(-2)
        
        % sort the difference
        [val, indexmax] = max(diff);
        [val, indexmin] = min(diff);
        D = val < diff;
        D(indexmax) = 0;
        [val, indexmid] = max(D);

        w_1 = W(indexmin) + W(indexmid) - W(indexmax);
        lambda_1 = Lambda(indexmin) + Lambda(indexmid) - Lambda(indexmax);

        if w_1 < 0
            w_1 = 0.01;
        end
        if lambda_1 < 0
            lambda_1 = 0.01;
        elseif lambda_1 > 0.95
            lambda_1 = 0.95;
        end

        diff_1 = Object(w_1,lambda_1, C_M, C_D, Zprob, z_grid,V2M, V2D, V20,V1);

        if diff_1 < diff(indexmin)
            w_2 = 1.5*(W(indexmin) + W(indexmid)) - 2*W(indexmax);
            lambda_2 = 1.5*(Lambda(indexmin) + Lambda(indexmid)) - ...
                2*Lambda(indexmax);
            
            if w_2 < 0
                w_2 = 0.01;
            end
            if lambda_2 < 0
                lambda_2 = 0.01;
            elseif lambda_2 > 0.95
                lambda_2 = 0.95;
            end
                
            diff_2 = Object(w_2,lambda_2, C_M, C_D, Zprob, z_grid,V2M, V2D, V20,V1);

            if diff_1 < diff_2
                diff(indexmax) = diff_1;
                W(indexmax) = w_1;
                Lambda(indexmax) = lambda_1;
            else
                diff(indexmax) = diff_2;
                W(indexmax) = w_2;
                Lambda(indexmax) = lambda_2;
            end
            
        elseif diff_1 > diff(indexmid)
            w_3 = 0.75*(W(indexmin) + W(indexmid)) - 0.5*W(indexmax);
            lambda_3 = 0.75*(Lambda(indexmin) + Lambda(indexmid)) - ...
                0.5*Lambda(indexmax);

            if w_3 < 0
                w_3 = 0.01;
            end
            if lambda_3 < 0
                lambda_3 = 0.01;
            elseif lambda_3 > 0.95
                lambda_3 = 0.95;
            end
            
            diff_3 = Object(w_3,lambda_3, C_M, C_D, Zprob, z_grid,V2M, V2D, V20,V1);
            
            if diff_3 < diff(indexmax)
                diff(indexmax) = diff_3;
                W(indexmax) = w_3;
                Lambda(indexmax) = lambda_3;
            else
                w_4 = 0.5*(W(indexmin) + W(indexmid));
                lambda_4 = 0.5*(Lambda(indexmin) + Lambda(indexmid));
                
                w_5 = 0.5*(W(indexmin) + W(indexmax));
                lambda_5 = 0.5*(Lambda(indexmin) + Lambda(indexmax));

                if w_4 < 0
                    w_4 = 0.01;
                end
                if lambda_4 < 0
                    lambda_4 = 0.01;
                elseif lambda_4 > 0.95
                    lambda_4 = 0.95;
                end

                if w_5 < 0
                    w_5 = 0.01;
                end
                if lambda_5 < 0
                    lambda_5 = 0.01;
                elseif lambda_5 > 0.95
                    lambda_5 = 0.95;
                end
                
                diff(indexmid) = Object(w_4,lambda_4, C_M, C_D, Zprob, ...
                    z_grid,V2M, V2D, V20,V1);
                diff(indexmax) =  Object(w_5,lambda_5, C_M, C_D, Zprob, ...
                    z_grid,V2M, V2D, V20,V1);
                W(indexmid) = w_4;
                W(indexmax) = w_5;
                Lambda(indexmid) = lambda_4;
                Lambda(indexmax) = lambda_5;
            end
        
        else
            diff(indexmax) = diff_1;
            W(indexmax) = w_1;
            Lambda(indexmax) = lambda_1;
        end
    
    epsi = diff(indexmax);
    disp(['W=',num2str(W(1)),num2str(W(2)),num2str(W(3))])
    disp(['Lambda=',num2str(Lambda(1)),num2str(Lambda(2)),num2str(Lambda(3))])
    disp(['epsi=', num2str(epsi)])
    time = toc;
    disp(['Iteration Time:',num2str(time)])

    if max(Lambda)-min(Lambda) < 10^(-5) & max(W)-min(W) < 10^-5 & epsi > 10^-2 

        Lambda = [rand(1,3)*0.95];
        W = [rand(1,3)];

    end
    
    end
end