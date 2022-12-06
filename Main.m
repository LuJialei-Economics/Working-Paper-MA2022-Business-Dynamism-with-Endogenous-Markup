%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Goal: 
%    This part of codes aims to solve the steady-state produtivity 
%    distribution and the value functions The numerial results provide the 
%    illustration of models/ mechanism
%
% Structure:
%    Part I : Set the values of parameters
%    Part II: Given lambda, calculate the steady-state distribution
%             The method is the discretization. 
%    Part III : Given the steady-state distribution, and the guess of 
%             wage and lambda, calculate the value fucntions.
%    Part IV  : Use the entry conditions to pin down the value of 
%             wage and lambda
%    Part V : Calculate the number of markets, N
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ========================== Part I: Parameters Settings ===================

clear all; clc;

global beta a theta delta rho num_state h;

% parameters
beta = 0.95;
a = 0.9;
theta = 1.3;
delta = 0.05;


% productivity grid
rho = 0.9;
num_state = 50;

% assume H(z) and h(z) is CDF and PDF of uniform distribution
h = ones(num_state,1)/num_state;

% generate the productivity grid
% use Tauchen's Method to transfer the AR(1) process into Markov
[Z, Zprob] = tauchen(0, 1);
z_grid = exp(Z);

%% ========================== Part II: Steady-State Distribution ============

tic;

% Calculate the Steady-state productivity distribution function of monopoly 
% markets, GM(z) and gM(z)

Lambda = 0.1:0.1:0.9;
gM_store = zeros(num_state,9);
GM_store = zeros(num_state,9);
for i=1:9
    [gM, GM] = GM_dist(Zprob, Lambda(i));
    gM_store(:,i) = gM;
    GM_store(:,i) = GM;
end

% Calculate the Steady-state joint productivity distribution function of 
% duopoly markets, GD(z) and gD(z)

gD_store = zeros(num_state,9);
GD_store = zeros(num_state,9);
for i=1:9
    [gD, GD, gJ] = GD_dist(Zprob, gM_store(:,i));
    gD_store(:,i) = gD;
    GD_store(:,i) = GD;
end

% plot

states = 1:1:num_state;

Label = [];
for i=1:9
    Label{i} = [append('\lambda=',num2str(Lambda(i)))];
end

figure
set(gca,'Fontsize',14);
subplot(2,2,1)
for i=1:9
    plot(states,gM_store(:,i),'LineWidth',1)
    xlabel("Monopolists' Productivity states")
    ylabel('PDF')
    title('Steady-state Productivity of Monopolists')
    hold on
end
legend(Label)

subplot(2,2,2)
for i=1:9
    plot(states,GM_store(:,i),'LineWidth',1)
    xlabel("Monopolists' Productivity states")
    ylabel('CDF')
    title('Steady-state Productivity of Monopolists')
    hold on
end
legend(Label)

subplot(2,2,3)
for i=1:9
    plot(states,gD_store(:,i),'LineWidth',1)
    xlabel("Duopolists' Productivity states")
    ylabel('PDF')
    title('Steady-state Productivity of Duopolists')
    hold on
end
legend(Label)

subplot(2,2,4)
for i=1:9
    plot(states,GD_store(:,i),'LineWidth',1)
    xlabel("Duopolists' Productivity states")
    ylabel('CDF')
    title('Steady-state Productivity of Duopolists')
    hold on
end
legend(Label)

% joint distribution plot
figure
X = 1:1:num_state;
Y = 1:1:num_state;
[X, Y] = meshgrid(X,Y);

Color = X.*sqrt(Y);
surf(X,Y, gJ, Color)
shading interp
title('Joint Distribution')

toc

%% ========================== Part III: Value Function Iteration ============

lambda = 0.5;
w = 0.7;

% Value fuunction iteration

% initialize the value functions
V1 = zeros(num_state, 1);
V2M = zeros(num_state, num_state);
V2D = zeros(num_state, num_state);
V20 = zeros(num_state, num_state);
itera = 0;

[V1_pri,V2_pri] = VF(w, lambda, z_grid, Zprob,V2M, V2D, V20,V1);

figure
subplot(1,2,1)
plot(z_grid,V1_pri)
xlabel('Asset')
ylabel('Value of Monopolists')
title('Value Function of Monopolists')

subplot(1,2,2)
[X,Y] = meshgrid(z_grid, z_grid);
Color = X.*sqrt(Y);
surfc(X,Y, V2_pri, Color)
shading interp
title('Value Function of Duopolists')

%% ========= Part IV: Calculate the equilibrium wage and lambda =============

% assume the entry coos of duopolists and monopolists
C_D = 5.6149;
C_M = 8.7646;
W = [0.7,0.73,0.75];
Lambda = [0.40,0.55,0.65];

Nelder_Mead(W, Lambda, C_M, C_D, Zprob, z_grid,V2M, V2D, V20,V1)


%% ========================= Subroutines ====================================




