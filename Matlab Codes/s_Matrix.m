function s_Matrix(glob, param)

    %S_MATRIX
    %
    % Given the productivity matrix [z_ijt, z_-ijt], genrate a
    % corresponding s-matrix and save
    %

    gamma = param.gamma;
    theta = param.theta;
    Nz = param.Nz;
    z_grid = glob.z_grid;

    SMatrix = zeros(Nz, Nz);

    for i = 1:Nz
        for j = 1:Nz
            fun = @(s) log(z_grid(i))-log(z_grid(j))+gamma*log(theta-s)-...
                gamma*log(theta-1+s)-(1-gamma)*log(s) + (1-gamma)*log(1-s);
            ini = z_grid(i)/(z_grid(i)+z_grid(j));
            x = fsolve(fun,ini);
            SMatrix(i,j) = x;
        end
    end

    for i = 1:Nz
        for j  = 1:Nz
            if imag(SMatrix(i,j)) < 10^(-10)
                SMatrix(i,j) = real(SMatrix(i,j));
            end
        end
    end
    
    save('SMatrix.mat',"SMatrix");
end

