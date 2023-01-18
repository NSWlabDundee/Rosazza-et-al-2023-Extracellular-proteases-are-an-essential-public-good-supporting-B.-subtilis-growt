%% Single strain, single model runs
% This script simulates the model for single-strain cultures across all
% three nutrient media and visualises the outcome by overlaying the growth
% curves of WT and D8 in all media.

% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 23/11/2022


clear; 
close all;

%% parameters
f1 = figure;
tmax = 100; % end time
tvec = [0,tmax];
options = odeset('MaxStep',1e-2); % options for ODE solver
run("parameters.m") % load parameters from file
% define nurtrient conditions
A0_col = [0.5,0.5,0.005]; 
N0_col = [0,50,50];
title_col = ["MS-Ga+Gly", "MS-Ga+BSA", "MS-BSA+Gly"];
col = 1/255*[64, 122, 255; 150, 252 ,242; 255, 46, 15; 255, 135, 125; 36, 199, 0; 133, 255, 23]; % colours for plots


%% initial conditions
ic_tot = 0.01; % total cell density in IC

for rr = 1:length(A0_col)
    ic_wt = [ic_tot,0,A0_col(rr),N0_col(rr),0,0,0,0]; % IC vector
    ic_d8 = [0,ic_tot,A0_col(rr),N0_col(rr),0,0,0,0]; % IC vector

    %% solve the system
    [t_wt, sol_wt] = ode15s(@(t,y) odesys(t,y,param), tvec, ic_wt, options);
    [t_d8, sol_d8] = ode15s(@(t,y) odesys(t,y,param), tvec, ic_d8, options);
%% visualisation

    wt = sol_wt(:,1); d8 = sol_d8(:,2);
    wtod  = (wt + sol_wt(:,6))/norm;
    d8od = (d8 + sol_d8(:,7))/norm;

    if rr == 1
        axes('Position',[.1 .1 0.25 0.9])
    elseif rr == 2
        axes('Position',[.4 .1 0.25 0.9])
    elseif rr == 3
        axes('Position',[.7 .1 0.25 0.9])
    end
    subplot(1,length(A0_col),rr)
    semilogy(t_wt,wtod, 'color', 'magenta')
    hold on
    grid on
    semilogy(t_d8,d8od, 'color', 'green')
    if rr <3
        legend("WT", "D8", 'location', 'southwest')
        legend boxoff
    elseif rr == 3
        legend("WT", "D8", 'location', 'northwest')
        legend boxoff
    end
    
    xlabel('Time, t')
    ylabel('Cell density')
    pbaspect([1 1 1])
    ylim([1e-3,1e1])
    boxwidth = 1;
    title(title_col(rr))
    if rr<3
        rectangle('Position',[0 1e-3 2*boxwidth 1e1-1e-3], 'EdgeColor','r')
        if rr == 1
            axes('Position',[.25 .2 .08 .4])
        elseif rr == 2
            axes('Position',[.53 .2 .08 .4])
        end
        box on
        semilogy(t_wt,wtod, 'color', 'magenta')
        hold on
        grid on
        semilogy(t_d8,d8od, 'color', 'green')
        set(gca, 'XTickLabel', [])
        set(gca, 'yTickLabel', [])
        xlim([0,boxwidth])
        ylim([1e-3,1e1])
        ax = gca;
        ax.XColor = 'red';
        ax.YColor = 'red';
    end
    
end    

set(f1,'Windowstyle','normal')
set(findall(f1,'-property','FontSize'),'FontSize',11)
set(f1,'Units','centimeters')
set(f1,'Position',[18 1 21 6])

