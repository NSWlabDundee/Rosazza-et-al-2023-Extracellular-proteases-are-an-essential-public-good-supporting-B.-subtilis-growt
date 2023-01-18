%% Initial cell ratio vs yield
% This script calculates the total yield in MSgg and MSbg for a range of
% different initial ratios of WT and D8 cells and visualises the output on
% a single plot.

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
col = 1/255*[66, 150, 255; 133, 255, 23];
tmax = 100; % output time
tt = linspace(0,tmax,1e3);
tvec = 0:1:tmax;
options = odeset('MaxStep',1e-2); % options for ode solver

ic_tot = 0.01; % total initial cell pop
wt_ic = linspace(0,ic_tot,21); % initial WT pop

f1 = figure;
grid on

for rr = 1:length(A0_col) % loop through all growth media
    yield_od = zeros(1,length(wt_ic));
    for ii = 1:length(wt_ic) % loop through all initial ratios
        ic_coex = [wt_ic(ii),ic_tot - wt_ic(ii),A0_col(rr),N0_col(rr),0,0,0,0]; % define IC
        [~, sol_coex(ii,:,:)] = ode15s(@(t,y) odesys(t,y,param), tt, ic_coex, options); % solve system
        wt = sol_coex(ii,end,1); d8 = sol_coex(ii,end,2);
        wtod  = wt + sol_coex(ii,end,6);
        d8od = d8 + sol_coex(ii,end,7);
        yield_od(ii) = (wt+wtod+d8+d8od)/norm; %calculate yield
    end
    
    p(rr) = plot(100-100*wt_ic/ic_tot,yield_od, '--o', 'color', col(rr,:), 'DisplayName', title_col(rr)); % plot
    hold on
    grid on
    
end


%% finalise visualisation
% legend(p,'location', 'northwest')
xlabel('Initial \Delta8 %')
ylabel('Yield')
xticks([0,25,50,75,100])
set(f1,'Windowstyle','normal')
set(findall(f1,'-property','FontSize'),'FontSize',11)
set(f1,'Units','centimeters')
set(f1,'Position',[10 5 10 8.4])

