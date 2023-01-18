%% Single model run for co-culture
% This script performs one single run for a co-culture experiment across
% MSgg and MSbg and visualises the outcome.

% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 23/11/2022


clear; 
close all;

%% parameters

wt_frac = 0.5; % WT fraction at start
f1 = figure;
tmax = 100; % end time
tvec = [0,tmax];

options = odeset('MaxStep',1e-2); % options for ode solver
run("parameters.m") % load parameters from file

% initialise nutrient conditions
A0_col = [0.5,0.005];
N0_col = [0,50];
title_col = ["MS-Ga+Gly", "MS-BSA+Gly"];

%% initial conditions
ic_tot = 0.01; % total cell density in IC

for rr = 1:length(A0_col)
    ic = [ic_tot*wt_frac,ic_tot-ic_tot*wt_frac,A0_col(rr),N0_col(rr),0,0,0,0]; % IC vector
    %% solve the system
    [t, sol] = ode15s(@(t,y) odesys(t,y,param), tvec, ic, options);
%% visualisation
    col = lines;

    wt = sol(:,1); d8 = sol(:,2);
    wtod  = (wt + sol(:,6))/norm;
    d8od = (d8 + sol(:,7))/norm;

    subplot(1,length(A0_col),rr)
    semilogy(t,wtod)
    hold on
    grid on
    semilogy(t,d8od)
    
    xlabel('Time, t')
    ylabel('Cell density')
    pbaspect([1 1 1])
    ylim([1e-3,1e1])
    
    title(title_col(rr))
    
    
    
end    
legend("WT", "D8", 'location', 'southeast')
set(f1,'Windowstyle','normal')
set(findall(f1,'-property','FontSize'),'FontSize',11)
set(f1,'Units','centimeters')
set(f1,'Position',[18 1 13.5 6])

