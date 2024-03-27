% ========================================================================
% Goal: 
%    This part of codes aims to solve the steady-state produtivity 
%    distribution and the value functions.
%    
%    The numerial results provide the illustration of models/ mechanism
%
% ========================================================================

%% Parameters Settings
clc, clear;
dbstop if error

%% options
options.distplot = 'Y'; % 'Y': plot graphs of productivity distribution, 
                        % 'N': not plot
options.vfplot = 'Y';   % 'Y': plot graphs of value functions, 
                        % 'N': not plot
options.calib  = 'N';   % 'Y': run calibration

% parameters
param.beta  = 0.92;     % discount factor
param.a     = 1/3000;   % disutility of labor supply
param.theta = 4;        % elasticity of substitution of goods on the 
                        % global markets
param.rho   = 0.9506;   % productivity persistence
param.mu    = 0;        % mean of innovations of productivity
param.sig   = 0.0865;   % standard error of innovations of productivity
                        % when running calibration, will be overwrite by 
                        % calibrated value
param.Nz    = 300;      % number of z's states
param.gamma = 0.80;     % Decreasing return to scale of labor

% global variables
glob.delta  = 0.09;    % exit rate of firms
glob.lambda = 0.5565;  % proportion of duopolistic markets
                       % when running calibration, will be overwrite by 
                       % calibrated value
glob.w      = 10;      % wage     

% generate the productivity distribution of entrants
% assume H(z) and h(z) is CDF and PDF of uniform distribution
glob.h = ones(param.Nz,1)/param.Nz;


%% Calibration
switch options.calib
    case 'Y'
        % Calibrate sigma, the standard error
        Phi1 = (param.theta-1)/(param.gamma+param.theta-param.theta...
            *param.gamma);
        std_gsale = sqrt(2)*0.7293;       % from Bottazzi and Secchi (2003b)
        param.sig = std_gsale /(Phi1*sqrt(2/(1-param.rho)));
end

% generate the productivity grid, since the method requires sigma, it
% should be run after calibrating sigma
% use Tauchen's Method to transfer the AR(1) process into Markov
[Z, Zprob]  = Tauchen(param);
glob.Z      = Z;
glob.Zprob  = Zprob;
glob.z_grid = exp(glob.Z);

% generate or load the market share matrix (S matrix)
if exist('SMatrix.mat','file')==2
    load('SMatrix.mat', 'SMatrix');
else
    s_Matrix(glob, param);
    load('SMatrix.mat', 'SMatrix');
end

switch options.calib
    case 'Y'
        % Calibrate lambda
        EMarkup = 1.28;         % from Christopoulou, Vermeulen (2008)
       
        % give an intuition about the reasonable range of lambda
        x = linspace(0.1,0.9,9);
        Mkup_duo = zeros(9,1);
        Mkup = zeros(9,1);
        gD_cali = zeros(param.Nz,9);
        gJ = zeros(param.Nz,param.Nz,9);
        for i = 1:9
            [gD_cali(:,i), ~, ~, ~] = SteadyState_Dist(glob,param,'lambda',x(i));
            gJ(:,:,i) = gD_cali(:,i) .* kron(gD_cali(:,i)',ones(param.Nz,1));
            Mkup_duo(i) = sum((param.theta*SMatrix./(param.theta-SMatrix)+...
                param.theta*(1-SMatrix)./(param.theta-1+SMatrix)).*...
                gJ(:,:,i),'all'); 
            Mkup(i) = (1-x(i))*param.theta/(param.theta-1) + x(i)*Mkup_duo(i);                
        end
        
        % calculate lambda by Bisection Method
        x1 = 0.2;
        x2 = 0.8;
        while abs(x1-x2)>1e-5
            [gD1, ~, ~, ~] = SteadyState_Dist(glob,param,'lambda',x1);
            gJ1 = gD1 .* kron(gD1',ones(param.Nz,1));
            Mkup_duo1 = sum((param.theta*SMatrix./(param.theta-SMatrix)+...
                param.theta*(1-SMatrix)./(param.theta-1+SMatrix)).*...
                gJ1,'all'); 
            Mkup1 = (1-x1)*param.theta/(param.theta-1) + x1*Mkup_duo1;
    
            [gD2, ~, ~, ~] = SteadyState_Dist(glob,param,'lambda',x2);
            gJ2 = gD2 .* kron(gD2',ones(param.Nz,1));
            Mkup_duo2 = sum((param.theta*SMatrix./(param.theta-SMatrix)+...
                param.theta*(1-SMatrix)./(param.theta-1+SMatrix)).*...
                gJ2,'all'); 
            Mkup2 = (1-x2)*param.theta/(param.theta-1) + x2*Mkup_duo2;
    
            xm = (x1+x2)/2;
            [gDm, ~, ~, ~] = SteadyState_Dist(glob,param,'lambda',xm);
            gJm = gDm .* kron(gDm',ones(param.Nz,1));
            Mkup_duom = sum((param.theta*SMatrix./(param.theta-SMatrix)+...
                param.theta*(1-SMatrix)./(param.theta-1+SMatrix)).*...
                gJm,'all'); 
            Mkupm = (1-xm)*param.theta/(param.theta-1) + xm*Mkup_duom;        
    
            if (Mkup1 - EMarkup) * (Mkupm - EMarkup) < 0
                x2 = xm;
            else
                x1 = xm;
            end
        end
        glob.lambda = xm;
end

%% Steady-State Distribution
% Calculate the Steady-state productivity distribution function of monopoly/
% duopoly markets GM(z) and gM(z)/ GD(z) and gD(z)

[result.gD, result.GD, result.gM, result.GM] = SteadyState_Dist(glob,param);

switch options.distplot
    case 'Y'
%       PDF of Monopolists' and Duopolists' Productivity
%       Figure (4) in Paper
        figure(4)
        set(gca,'FontSize',14, 'FontName','Times',...
            'LooseInset', max(get(gca,'TightInset'), 1.02));
        plot(glob.z_grid,result.gM,'LineWidth',1)
        hold on
        plot(glob.z_grid,result.gD,'LineWidth',1)    
        hold off
        legend('Monopolists','Duopolists')
        xlabel("Productivity",Interpreter='latex')
        ylabel('PDF',Interpreter='latex')

%       CDF of Monopolists' and Duopolists' Productivity
        figure(99)
        set(gca,'FontSize',14, 'FontName','Times',...
            'LooseInset', max(get(gca,'TightInset'), 0.02));
        plot(glob.z_grid,result.GM,'LineWidth',1)
        hold on
        plot(glob.z_grid,result.GD,'LineWidth',1)
        hold off
        legend('Monopolists','Duopolists')
        xlabel("Productivity",Interpreter='latex')
        ylabel('CDF',Interpreter='latex')
end



%% Value Function iteration
% given the values of wage (w) and proportion of Duopoly (lambda), 
% give out value functions
% Notions: 
% when as input, lambda and w will be independent  variables, 
% when as output, they will be included in the structure: result

% iterate to get value function
[result.V1_pri,result.V2_pri,result.pi1, result.pi2] = VF(SMatrix, glob, param);

result.EV2_cond = result.V2_pri*result.gM;


switch options.vfplot
    case 'Y'
%       Value Function of Monopolists Before the Shock
%       Figure (1) in paper
        figure(1)
        plot(glob.z_grid,result.V1_pri)
        xlabel('Productivity',Interpreter='latex')
        ylabel('Values of Monopolists $V_1(z)$',Interpreter='latex')
%       Value Function of Duopolist i Before the Shock
%       Figure (2) in paper      
        figure(2)        
        [X,Y] = meshgrid(glob.z_grid, glob.z_grid);
        surfc(X,Y, result.V2_pri,'FaceAlpha',0.5)
        xlabel('Productivity of firm -i',Interpreter='latex')
        ylabel('Productivity of firm i',Interpreter='latex')
        zlabel('Values of Duopolists $V_2(z,z_{-})$',Interpreter='latex')
        shading interp
%       Comparison Value Function of Monopolists and Expected Value 
%       Function of Duopolist i Before the Shock
%       Figure (3) in paper
        figure(3)
        plot(glob.z_grid,result.EV2_cond)
        hold on
        plot(glob.z_grid,result.V1_pri)
        xlabel('Productivity',Interpreter='latex')
        ylabel('Values $V_1$ v.s. Expected Values $E_{z{-}}V_2$','interpreter','latex');
        legend('Expected Value, Duopolists','Value, Monopolists', ...
            Location='northwest')
end

%% Original Entry Costs and Number of Industries 'N'
EV1 = glob.h'*result.V1_pri;  % expected value of entering as a monopolist
EV2 = glob.h'*result.V2_pri*result.gM; % expected value of entering as a duopolist
result.C_D = EV2;
result.C_M = EV1;

xi = param.gamma + param.theta - param.theta*param.gamma;
gJ = result.gD .* kron(result.gD',ones(param.Nz,1));
z_pair = zeros(param.Nz, param.Nz, 2);
z_pair(:,:,2) = kron(glob.z_grid, ones(1, param.Nz));
z_pair(:,:,1) = kron(glob.z_grid', ones(param.Nz,1));
temp1 = z_pair(:,:,1).^((param.theta-xi)/xi/param.gamma) .* SMatrix.^(1/xi).*...
    (param.theta./(param.theta - SMatrix)).^(-param.theta/xi);
temp2 = z_pair(:,:,2).^((param.theta-xi)/xi/param.gamma) .* (1-SMatrix).^(1/xi).*...
    (param.theta./(param.theta - 1 + SMatrix)).^(-param.theta/xi);
rhs = param.gamma^(param.theta/xi)*(1/param.a)^(1/xi)*glob.w^((1-param.theta)/xi)*...
    ((1-glob.lambda)*(param.theta/(param.theta-1))^(-param.theta/xi)*...
    glob.z_grid'.^((param.theta-xi)/xi/param.gamma)*result.gM +...
    glob.lambda * sum((temp1+temp2).*gJ,'all'));
N = 1/param.a/rhs;

%% Markup of Duopolists
result.gJ = result.gD .* kron(result.gD',ones(param.Nz,1));
result.Mkup_duo = sum((param.theta*SMatrix./(param.theta-SMatrix)+...
            param.theta*(1-SMatrix)./(param.theta-1+SMatrix)).*...
            result.gJ,'all'); 
result.Mkup = (1-glob.lambda)*param.theta/(param.theta-1) +...
    glob.lambda*result.Mkup_duo; 

