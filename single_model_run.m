%% Single model run
% This script simulates the model once.

% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 23/11/2022


clear; 
close all;

%% parameters
tmax = 100; % end time
tvec = [0,tmax];
options = odeset('MaxStep',1e-2); % options for ode solver
% select nutrient medium
% A0 = 0.5; N0 = 0; %MSgg
A0 = 0.005; N0 =50; %MSbg
% A0 = 0.5; N0 =50; %MS-GA + BSA
run("parameters.m") % load parameters from file

%% initial conditions
ic_tot = 0.01; % total cell density in IC
wt_ic = ic_tot*0.5; % fraction of initial WT
ic_coex = [wt_ic,ic_tot - wt_ic,A0,N0,0,0,0,0]; % IC vector

%% solve the system
[t, sol_coex] = ode15s(@(t,y) odesys(t,y,param), tvec, ic_coex, options);

%% visualisation
col = lines;
f1 = figure;
subplot(3,2,1)

wt = sol_coex(:,1); d8 = sol_coex(:,2);
wtod  = wt + sol_coex(:,6);
d8od = d8 + sol_coex(:,7);

semilogy(t,wtod)
hold on
grid on
semilogy(t,d8od)
legend("WT", "D8", 'location', 'southeast')
xlabel('Time, $t$', 'interpreter','latex')
ylabel('OD$_{600}$', 'interpreter','latex')
pbaspect([1 1 1])

subplot(3,2,2)
hold on
grid on

plot(t,sol_coex(:,6)./wtod)
plot(t,sol_coex(:,7)./d8od)
legend("WT", "D8 ", 'location', 'southeast')
xlabel('Time, $t$', 'interpreter','latex')
ylabel('$\%$ spores', 'interpreter','latex')
pbaspect([1 1 1])

subplot(3,2,3)
hold on
grid on
plot(t,sol_coex(:,3))
xlabel('Time, $t$', 'interpreter','latex')
ylabel('Available nutrient, $A$', 'interpreter','latex')
pbaspect([1 1 1])

subplot(3,2,4)
hold on
grid on
plot(t,sol_coex(:,4))
xlabel('Time, $t$', 'interpreter','latex')
ylabel('Base nutrient, $N$', 'interpreter','latex')
pbaspect([1 1 1])


subplot(3,2,5)
hold on
grid on
plot(t,sol_coex(:,5))
xlabel('Time, $t$', 'interpreter','latex')
ylabel('Exoproteases, $E$', 'interpreter','latex')
pbaspect([1 1 1])

subplot(3,2,6)
hold on
grid on
plot(t,sol_coex(:,8))
xlabel('Time, $t$', 'interpreter','latex')
ylabel('Intermediate nutrient, $A_0$', 'interpreter','latex')
pbaspect([1 1 1])



    
set(f1,'Windowstyle','normal')
set(findall(f1,'-property','FontSize'),'FontSize',11)
set(f1,'Units','centimeters')
set(f1,'Position',[18 1 15 16])

