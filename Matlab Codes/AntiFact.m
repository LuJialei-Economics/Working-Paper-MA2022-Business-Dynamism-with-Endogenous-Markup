%% Anti-Fact Experiments
% Before runing this program, run Main.m first
% Benchmark
% result.C_M = 45.2019; result.C_D = 44.6412;
% Mkup = 1.28; Mkup_duo = 1.2375 
% glob.w = 10;  glob.lambda = 0.5565

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Decrease of C_D:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

calflg = 'N';

% new equilibrium lambda and wage
C_D = 0.99*result.C_D;    % entry cost of duopolistic market decrease 1%

switch calflg
    case 'Y'  % run non-linear solver to solve the new w and lambda
        x0 = [0.2; 10];  % x = [lambda w]
        option = optimset('Display','iter','PlotFcns',@optimplotfval);
        xDD = fminsearch(@Object,x0,option, result.C_M, C_D, glob,...
            param);
        lambdaDD = xDD(1);
        wDD = xDD(2);
    case 'N'
        lambdaDD = 0.3598;
        wDD = 10.1887;
end

%% Check the result
[gD, GD, gM, GM] = SteadyState_Dist(glob,param,'lambda',lambdaDD);
[V1_DD, V2_DD,~,~] = VF(SMatrix, glob, param,'lambda',lambdaDD,'w',wDD);
EV1DD = glob.h'*V1_DD;  % expected value of entering as a monopolist
EV2DD = glob.h'*V2_DD*gM; % expected value of entering as a duopolist

%% Markup distribution
gJ = gD .* kron(gD',ones(param.Nz,1));
Mkup_duoDD = sum((param.theta*SMatrix./(param.theta-SMatrix)+...
    param.theta*(1-SMatrix)./(param.theta-1+SMatrix)).*gJ,'all');  
MkupDD = (1-lambdaDD)*param.theta/(param.theta-1) + lambdaDD*Mkup_duoDD; 

% Number of markets after schocks
xi = param.gamma + param.theta - param.theta*param.gamma;
z_pair = zeros(param.Nz, param.Nz, 2);
z_pair(:,:,2) = kron(glob.z_grid, ones(1, param.Nz));
z_pair(:,:,1) = kron(glob.z_grid', ones(param.Nz,1));
temp1 = z_pair(:,:,1).^((param.theta-xi)/xi/param.gamma) .* SMatrix.^(1/xi).*...
    (param.theta./(param.theta - SMatrix)).^(-param.theta/xi);
temp2 = z_pair(:,:,2).^((param.theta-xi)/xi/param.gamma) .* (1-SMatrix).^(1/xi).*...
    (param.theta./(param.theta - 1 + SMatrix)).^(-param.theta/xi);
rhs = param.gamma^(param.theta/xi)*(1/param.a)^(1/xi)*wDD^((1-param.theta)/xi)*...
    ((1-lambdaDD)*(param.theta/(param.theta-1))^(-param.theta/xi)*...
    glob.z_grid'.^((param.theta-xi)/xi/param.gamma)*result.gM +...
    lambdaDD * sum((temp1+temp2).*gJ,'all'));
NDD = 1/param.a/rhs;  


%% plot
% Value Function of Monopolists After the Shock
% Figure (15) in paper
figure(15)
plot(glob.z_grid,V1_DD)
xlabel('Productivity',Interpreter='latex')
ylabel('Values of Monopolists $V_1(z)$',Interpreter='latex')

% Value Function of Duopolist i After the Shock
% Figure (16) in paper   
figure(16)        
[X,Y] = meshgrid(glob.z_grid, glob.z_grid);
surfc(X,Y, V2_DD,'FaceAlpha',0.5)
xlabel('Productivity of firm -i',Interpreter='latex')
ylabel('Productivity of firm i',Interpreter='latex')
zlabel('Values of Duopolists $V_2(z,z_{-})$',Interpreter='latex')
shading interp

EV2DD_cond = V2_DD*gM;  % expected duopolistic values after shock


% Comparison Value Function of Monopolists and Expected Value 
% Function of Duopolist i Before and After the Shock
% Figure (17) in paper
figure(17)
plot(glob.z_grid,EV2DD_cond,color = 'red')
hold on
plot(glob.z_grid,V1_DD,color = 'blue')
plot(glob.z_grid,result.EV2_cond, 'LineStyle','--', color = 'red')
plot(glob.z_grid,result.V1_pri,'LineStyle','--', color = 'blue')
xlabel('Productivity',Interpreter='latex')
ylabel('Values $V_1$ v.s. Expected Values $E_{z{-}}V_2$','interpreter','latex');
legend('Expected Value, Duopolists After','Value, Monopolists After', ...
    'Expected Value, Duopolists Before','Value, Monopolists Before', ...
    Location='northwest')

% PDF of Monopolists' and Duopolists' Productivity Before and After
% Figure (18) in paper
figure(18)
set(gca,'FontSize',14, 'FontName','Times',...
    'LooseInset', max(get(gca,'TightInset'), 1.02));
plot(glob.z_grid,result.gM,'LineWidth',1,'LineStyle','--', color = 'blue')
hold on
plot(glob.z_grid,result.gD,'LineWidth',1,'LineStyle','--', color = 'red')    
plot(glob.z_grid,gM,'LineWidth',1, color = 'blue')
hold on
plot(glob.z_grid,gD,'LineWidth',1, color = 'red')  
hold off
legend('PDF of Monopolists Before Shocks','Marginal PDF of Duopolists Before Shocks',...
    'PDF of Monopolists After Shocks','Marginal PDF of Duopolists  After Shocks')
xlabel("Productivity",Interpreter='latex')
ylabel('PDF',Interpreter='latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Increase of C_D:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calflg = 'N';

% new equilibrium lambda and wage
C_D = 1.01*result.C_D;    % entry cost of duopolistic market increase 1%

switch calflg
    case 'Y'  % run non-linear solver to solve the new w and lambda
        x0 = [0.2; 10];  % x = [lambda w]
        option = optimset('Display','iter','PlotFcns',@optimplotfval);
        xID = fminsearch(@Object,x0,option, result.C_M, C_D, glob,...
            param);
        lambdaID = xID(1);
        wID = xID(2);
    case 'N'
        lambdaID = 0.7403;
        wID = 9.7863;
end

%% Check the result
[gD, GD, gM, GM] = SteadyState_Dist(glob,param,'lambda',lambdaID);
[V1_ID, V2_ID,~,~] = VF(SMatrix, glob, param,'lambda',lambdaID,'w',wID);
EV1ID = glob.h'*V1_ID;  % expected value of entering as a monopolist
EV2ID = glob.h'*V2_ID*gM; % expected value of entering as a duopolist

%% Markup distribution
gJ = gD .* kron(gD',ones(param.Nz,1));
Mkup_duoID = sum((param.theta*SMatrix./(param.theta-SMatrix)+...
    param.theta*(1-SMatrix)./(param.theta-1+SMatrix)).*gJ,'all');  %1.2367
MkupID = (1-lambdaID)*param.theta/(param.theta-1) + lambdaID*Mkup_duoID; %1.2618

% Number of markets after schocks
rhs = param.gamma^(param.theta/xi)*(1/param.a)^(1/xi)*wID^((1-param.theta)/xi)*...
    ((1-lambdaID)*(param.theta/(param.theta-1))^(-param.theta/xi)*...
    glob.z_grid'.^((param.theta-xi)/xi/param.gamma)*result.gM +...
    lambdaID * sum((temp1+temp2).*gJ,'all'));
NID = 1/param.a/rhs;  %2926

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Decrease of C_M:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calflg = 'N';

% new equilibrium lambda and wage
C_M = 0.99*result.C_M;    % entry cost of duopolistic market increase 1%

switch calflg
    case 'Y'  % run non-linear solver to solve the new w and lambda
        x0 = [0.2; 10];  % x = [lambda w]
        option = optimset('Display','iter','PlotFcns',@optimplotfval);
        xDM = fminsearch(@Object,x0,option, C_M, result.C_D, glob,...
            param);
        lambdaDM = xDM(1);
        wDM = xDM(2);
    case 'N'
        lambdaDM = 0.7424;
        wDM = 9.8967;
end

%% Check the result
[gD, GD, gM, GM] = SteadyState_Dist(glob,param,'lambda',lambdaDM);
[V1_DM, V2_DM,~,~] = VF(SMatrix, glob, param,'lambda',lambdaDM,'w',wDM);
EV1DM = glob.h'*V1_DM;  % expected value of entering as a monopolist
EV2DM = glob.h'*V2_DM*gM; % expected value of entering as a duopolist

%% Markup distribution
gJ = gD .* kron(gD',ones(param.Nz,1));
Mkup_duoDM = sum((param.theta*SMatrix./(param.theta-SMatrix)+...
    param.theta*(1-SMatrix)./(param.theta-1+SMatrix)).*gJ,'all');  %1.2367
MkupDM = (1-lambdaID)*param.theta/(param.theta-1) + lambdaDM*Mkup_duoDM; %1.2644

% Number of markets after schocks
rhs = param.gamma^(param.theta/xi)*(1/param.a)^(1/xi)*wID^((1-param.theta)/xi)*...
    ((1-lambdaDM)*(param.theta/(param.theta-1))^(-param.theta/xi)*...
    glob.z_grid'.^((param.theta-xi)/xi/param.gamma)*result.gM +...
    lambdaDM * sum((temp1+temp2).*gJ,'all'));
NDM = 1/param.a/rhs;  %2926

