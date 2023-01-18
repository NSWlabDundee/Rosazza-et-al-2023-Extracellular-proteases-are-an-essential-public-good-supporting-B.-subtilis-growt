%% Growth curves for different initial WT populations
% This script simulates co-cultures for different initial fractions of the
% WT population and visualises the final WT pop and time dynamics.

% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 23/11/2022



clear; 
close all;

%% parameters
run("parameters.m") % load parameters from file

A0_col = [0.5,0.005]; % growth media definition
N0_col = [0,50];
title_col = ["MS-Ga+Gly", "MS-BSA+Gly"];
icwt = 0:5:100; % range of initial WT percentage
% icwt = [5,20,50,70,90,95]; % range of initial WT percentage

tmax = 100;
tt = linspace(0,tmax);
tvec = 0:1:tmax;
options = odeset('MaxStep',1e-2);

f1 = figure; f = figure;
col1 = 1/255*[66, 150, 255; 133, 255, 23];
for rr = 1:length(A0_col) % loop through all growth media
    figure(f)
    wt_frac = NaN*ones(1,length(icwt));
    for ii = 1:length(icwt)
    
        ic = [0.01*icwt(ii)/100,0.01*(100-icwt(ii))/100,A0_col(rr),N0_col(rr),0,0,0,0]; % initial condition
        [~, sol] = ode15s(@(t,y) odesys(t,y,param), tt, ic, options);
        wt_frac(ii) = (sol(end,1)+sol(end,6))./(sol(end,1)+sol(end,6)+sol(end,7)+sol(end,2));
    
    %% plots
        col = lines;
        noslots = length(wt_frac);
        subplot(2,ceil(noslots/2),ii)
        
        hold on
        grid on
        p(1) = plot(tt, 100*(sol(:,1)+sol(:,6))./(sol(:,1)+sol(:,6)+sol(:,7)+sol(:,2)), 'color',col(5,:), 'displayname','WT');
        p(2) = plot(tt, 100*(sol(:,2)+sol(:,7))./(sol(:,1)+sol(:,6)+sol(:,7)+sol(:,2)),  'color',col(3,:),'displayname','$\Delta 8$');
        xlabel('Time', 'interpreter','latex')
        if ii == 1 || ii == ceil(ceil(noslots/2))+1
            ylabel('$\%$ population', 'interpreter','latex')
        end
        title([num2str(icwt(ii)),'$\%$'], 'interpreter','latex')
        if ii == 1
            leg = legend(p,'location','east','NumColumns',1, 'Orientation', 'vertical');
            set(leg,'interpreter','latex')
        end
        ylim([0,100])
        pbaspect([1 1 1])
        set(leg,'interpreter','latex')
    end
    
    
    
    set(f,'Windowstyle','normal')
    set(findall(f,'-property','FontSize'),'FontSize',10.5)
    set(f,'Units','centimeters')
    set(f,'Position',[10 5 12.5 17/2])
    
    %% figure final WT % vs initial WT %
    figure(f1)
    hold on
    p1(rr) = plot(100 - icwt,100 - 100*wt_frac, '--o','color', col1(rr,:), 'DisplayName', title_col(rr));
    
end
figure(f1)
% legend(p1,'location', 'northwest')
xlabel('Initial \Delta8 %')
ylabel('Final \Delta8 %')
grid on
xticks([0,25,50,75,100])
yticks([0,25,50,75,100])
set(f1,'Windowstyle','normal')
set(findall(f1,'-property','FontSize'),'FontSize',10.5)
set(f1,'Units','centimeters')
set(f1,'Position',[10 5 10 8.4])
