% Analysis and Mechanism
% run after runing Main.m and AntiFact.m

%% ----------------------------------------------------------------------
% Relationship between wages and other variables
% -----------------------------------------------------------------------

%% Detect the relationship between wages and difference between V1 and V2
Wage = 6:1:14;
V1_storew = zeros(param.Nz, length(Wage));
V2_storew = zeros(param.Nz, param.Nz, length(Wage));
EV2_storew = zeros(param.Nz, length(Wage));
[gDAM, GDAM, gMAM, GMAM] = SteadyState_Dist(glob,param);
gJAM = gDAM.*gDAM';
for i = 1:length(Wage)
    [V1_storew(:,i),V2_storew(:,:,i),~,~] =...
        VF(SMatrix, glob, param, 'w', Wage(i));
    EV2_storew(:,i) = V2_storew(:,:,i)*gMAM;
end

Labelw = [];
for j=1:2:length(Wage)
    i = int64(ceil(j/2));
    Labelw{i} = [append('w=',num2str(Wage(j)))];
end

% Difference between Value of Monopolists and Expected Value of Duopolists
% given Different wage
% Figure (7) in paper
diff_VFw12 = V1_storew - EV2_storew;
figure(7)
plot(glob.z_grid,diff_VFw12(:,1:2:end))
xlabel('Productivity',Interpreter='latex')
ylabel('Value Differences ($V_1-E_{z_{-}}V_2$) Given $\lambda$'...
    ,'interpreter','latex')
hl = refline([0,0]);
hl.Color = 'k';
legend(Labelw)


%% Detect the relationship between wages and EXPECTED difference between 
%% V1 and V2

Ediff_VFw12 = glob.h' * diff_VFw12;
% Expected Value Differences w.r.t different wage levels
% Figure (10) in paper
figure(10)
plot(Wage,Ediff_VFw12, LineWidth=1)
xlabel('Wage Level',Interpreter='latex')
ylabel('EXPECTED Value Differences ($E_z[V_1-E_{z_{-}}V_2]$) Given $\lambda$'...
    ,'interpreter','latex')
legend('EV1-EEV2')

%% Detect the relationship between wages and V1 and EV2

% Monopolistic and Expected Duopolistic Values w.r.t different wage levels
% Figure (5) in paper
figure(5)
figureSize = [488, 242, 1120, 420]; % [left, bottom, width, height]
set(gcf, 'Position', figureSize);
subplot(1,2,1)
plot(glob.z_grid, V1_storew(:,1:2:end),'LineWidth',1);
xlabel('Productivity',Interpreter='latex')
ylabel('Values of Monopolists $V_1$','interpreter','latex')
legend(Labelw,Location='northwest')
colors = get(get(gca, 'Children'), 'Color');
subplot(1,2,2)
plot(glob.z_grid, EV2_storew(:,1:2:end),'LineWidth',1);
xlabel('Productivity',Interpreter='latex')
ylabel('Expected Values of Duopolists $E_{z_{-}}V_2$','interpreter','latex')
legend(Labelw,Location='northwest')

%% ----------------------------------------------------------------------
% Relationship between lambda and other variables
% -----------------------------------------------------------------------

%% Detect the relationship between lambda and distributions

Lambda = 0.2:0.15:0.8;
gMAM_store = zeros(param.Nz,length(Lambda));
EM_store = zeros(1,length(Lambda));
gDAM_store = zeros(param.Nz,length(Lambda));
for j = 1:length(Lambda)
    [gDAM, ~, gMAM, ~] = SteadyState_Dist(glob,param, 'lambda',Lambda(j));
    gMAM_store(:,j) = gMAM;
    gDAM_store(:,j) = gDAM;
    EM_store(j) = glob.z_grid' * gMAM;
end

Label = [];
for j=1:length(Lambda)
    Label{j} = [append('\lambda=',num2str(Lambda(j)))];
end

% Productivity Distribution of Monopolists and corresponding gravity
% centers after the shock
% Figure (9) in paper
figure(12)
figureSize = [488, 242, 1120, 420]; % [left, bottom, width, height]
set(gcf, 'Position', figureSize);
subplot(1,2,1)
plot(glob.z_grid, gMAM_store,'LineWidth',1);
xlabel('Productivity',Interpreter='latex')
ylabel('PDF',Interpreter='latex')
legend(Label)
colors = get(get(gca, 'Children'), 'Color');
subplot(1,2,2)
for j = 1:length(Lambda)
    xline(EM_store(j),'Color', colors{6-j},'LineWidth',1)
end
xlabel('Mean Productivity',Interpreter='latex')
legend(Label)


%% Detect the relationship between lambdas and difference between V1 and V2

V1_store = zeros(param.Nz, length(Lambda));
V2_store = zeros(param.Nz, param.Nz, length(Lambda));
EV2_store = zeros(param.Nz, length(Lambda));

for i = 1:length(Lambda)
    [V1_store(:,i),V2_store(:,:,i),~,~] =...
        VF(SMatrix, glob, param, 'lambda', Lambda(i));
    EV2_store(:,i) = V2_store(:,:,i)*gMAM_store(:,i);
end

% Difference between Value of Monopolists and Expected Value of Duopolists
% given Different lambda
% Figure (8) in paper
diff_VF12 = V1_store - EV2_store;
figure(14)
plot(glob.z_grid,diff_VF12)
xlabel('Productivity',Interpreter='latex')
ylabel('Value Differences ($V_1-E_{z_{-}}V_2$) Given w' ...
    ,'interpreter','latex')
hl = refline([0,0]);
hl.Color = 'k';
legend(Label)


%% Detect the relationship between lambdas and V1 and EV2

% Monopolistic and Expected Duopolistic Values w.r.t different lambda levels
% Figure (6) in paper
figure(6)
figureSize = [488, 242, 1120, 420]; % [left, bottom, width, height]
set(gcf, 'Position', figureSize);
subplot(1,2,1)
plot(glob.z_grid, V1_store./V1_store(:,1),'LineWidth',1);
xlabel('Productivity',Interpreter='latex')
label2 = sprintf('Proportion Changes of Values of Moopolists\n $V_1/V_1(\\lambda=0.2)$')
ylabel(label2,Interpreter='latex')
legend(Label,Location='northwest')
colors = get(get(gca, 'Children'), 'Color');
subplot(1,2,2)
plot(glob.z_grid, EV2_store./EV2_store(:,1),'LineWidth',1);
xlabel('Productivity',Interpreter='latex')
labels = sprintf('Proportion Changes of Expected Values of Duopolists \n $EV_2/EV_2(\\lambda=0.2)$')
ylabel(labels,Interpreter='latex');
legend(Label,Location='northwest')


%% Detect the relationship between lambdas and EXPECTED difference between 
%% V1 and V2

% Expected Value Differences w.r.t different lambda levels
% Figure (11) in paper
Ediff_VF12 = glob.h' * diff_VF12;
figure(11)
plot(Lambda,Ediff_VF12, LineWidth=1)
xlabel('$\lambda$',Interpreter='latex')
ylabel('EXPECTED Value Differences ($E_z[V_1-E_{z_{-}}V_2]$) Given $\lambda$',...
    Interpreter='latex')

%% ----------------------------------------------------------------------
% Relationship between lambda and EEV1, EEV2
% -----------------------------------------------------------------------

LambdaM = 0.3:0.01:0.6;
WageM = [10, 10.1887];

gMM_store = zeros(param.Nz,length(LambdaM));
for i = 1:length(LambdaM)
    [~, ~, gMM, ~] = SteadyState_Dist(glob,param, 'lambda',LambdaM(i));
    gMM_store(:,i) = gMM;
end

V1M_store = zeros(param.Nz, length(LambdaM),length(WageM));
V2M_store = zeros(param.Nz, param.Nz, length(LambdaM),length(WageM));
EEV1_store = zeros(length(LambdaM),length(WageM));
EEV2_store = zeros(length(LambdaM),length(WageM));

for j = 1:length(WageM)
    for i = 1:length(LambdaM)
        [V1M_store(:,i,j),V2M_store(:,:,i,j),~,~] =...
        VF(SMatrix, glob, param, 'lambda', LambdaM(i),'w',WageM(j));
    end
end

for j = 1:length(WageM)
    for i = 1:length(LambdaM)
        EEV1_store(i,j) = glob.h'*V1M_store(:,i,j);
        EEV2_store(i,j) = glob.h'*V2M_store(:,:,i,j)*gMM_store(:,i);
    end
end

%%
% Expected monopolistic values and expected duopolistic values' expectation
% faced by entrants with repsect to different lambda and wage levels 
% before and after the shock
% Figure (12) in paper
figure(12)
figureSize = [488, 0, 560, 840]; % [left, bottom, width, height]
set(gcf, 'Position', figureSize);
subplot(2,1,1)
plot(LambdaM, EEV1_store(:,1),'LineWidth',1);
hold on
plot(LambdaM, EEV1_store(:,2),'LineWidth',1);
yline(result.C_M, 'LineWidth',1);
xline(glob.lambda,'LineWidth',1)
xline(lambdaDD,'LineWidth',1,LineStyle='-.')
xlabel('$\lambda$',Interpreter='latex')
ylabel('Expected Value of Monopolists $E_zV_1(z)$',Interpreter='latex')
legend('EV1 before Shock','EV1 after Shock')
colors = get(get(gca, 'Children'), 'Color');

subplot(2,1,2)
plot(LambdaM, EEV2_store(:,1),'LineWidth',1);
hold on
plot(LambdaM, EEV2_store(:,2),'LineWidth',1);
yline(result.C_D, 'LineWidth',1)
yline(C_D, 'LineWidth',1, LineStyle='--')
xline(glob.lambda,'LineWidth',1)
xline(lambdaDD,'LineWidth',1,LineStyle='-.')
xlabel('$\lambda$',Interpreter='latex')
ylabel('Expected Values of Duopolists $E_zE_{z_{-}}V_2(z,z_{-})$',Interpreter='latex')
legend('EEV2 before Shock','EEV2 after Shock')
%title('Productivity Distribution of Monopolists')


% Expected value difference before and after shock
% Figure (13) in paper
figure(13)
plot(LambdaM, EEV1_store(:,1) - EEV2_store(:,1),'LineWidth',1);
hold on
plot(LambdaM, EEV1_store(:,2) - EEV2_store(:,2),'LineWidth',1);
xline(glob.lambda,'LineWidth',1)
xline(lambdaDD,'LineWidth',1,LineStyle='-.')
yline(result.C_M-result.C_D, 'LineWidth',1)
yline(result.C_M-C_D, 'LineWidth',1, LineStyle='--')
xlabel('$\lambda$',Interpreter='latex')
ylabel('Expected Value Differneces $E_z(V_1-E_{z_{-}}V_2)$',Interpreter='latex')
legend('Expected Value Differneces before Shock',...
    'Expected Value Differneces after Shock')

%% Market share distributions w.r.t different lambda
dmarketshare = zeros(param.Nz, param.Nz);
dmarketshare_GJ = zeros(param.Nz*param.Nz,3);

for i = 1:param.Nz
    for j = 1:param.Nz
        dmarketshare(i,j) = abs(SMatrix(i) - SMatrix(j)); 
               % differences between maket shares in the same local market
    end
end

dz = zeros(param.Nz, param.Nz);
Dz = zeros(param.Nz*param.Nz,length(Lambda)+1);

for i = 1:param.Nz
    for j = 1:param.Nz
        dz(i,j) = abs(SMatrix(i) - SMatrix(j));
    end
end

dmarketshare_GJ(:,1) = dz(:);
dmarketshare_GJ(:,2) = result.gJ(:);
dmarketshare_GJ(:,3) = gJ(:);
dshare_unique = unique(dmarketshare_GJ(:,1));

pdf0 = zeros(length(dshare_unique),2);

for j = 1:2
    for i = 1:length(dshare_unique)
        pdf0(i,j) = sum(dmarketshare_GJ(dmarketshare_GJ(:,1) ==...
            dshare_unique(i), 1+j));
    end
end

%%
% Changes on PDF of market share difference before and after shock
% Figure (14) in paper
figure(14)
plot(dshare_unique,pdf0(:,1)-pdf0(:,2))
xlabel('Market Share Difference $|s_i-s_{-i}|$',Interpreter='latex')
ylabel("PDF Differences $\psi-\psi'$",Interpreter='latex')

