%% Co-cultures with different EP production costs
% This script simulates co-cultures in MSbg for a range of different EP
% production costs and initial ratios and visualises the outcome on a
% single plot

% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 23/11/2022


clear; 
close all;

%% parameters

wt_frac_col = [0.05,.2,.5,.7,.9,.95, 1]; % WT fraction at start
tmax = 100; % end time
tvec = linspace(0,tmax,100000);
options = odeset('MaxStep',1e-2);
run("parameters.m") % load parameters from file
EP_cost_col = 0:0.05:2; %range of EP production costs

f1 = figure(1);
f2 = figure(2);
f3 = figure(3);
f4 = figure(4);
f5 = figure(5);



%% initial conditions
ic_tot = 0.01; % total cell density in IC
Aic = 0.005; Nic = 50; % MSbg
% Aic = 0.5; Nic = 0; % MSgg

for ww = 1:length(wt_frac_col)
    wt_frac = wt_frac_col(ww);
    ic = [ic_tot*wt_frac,ic_tot-ic_tot*wt_frac,Aic,Nic,0,0,0,0]; % IC vector
    
    WT_frac_end = NaN*ones(1,length(EP_cost_col));
    yield = NaN*ones(1,length(EP_cost_col));
    fm1 = figure;
    f0 = figure;
    for rr = 1:length(EP_cost_col)
        param(13) = EP_cost_col(rr);
        %% solve the system
        [t, sol] = ode15s(@(t,y) odesys(t,y,param), tvec, ic, options);
    %% visualisation
        col = lines;

        wt = sol(:,1); d8 = sol(:,2);
        wtod  = (wt + sol(:,6))/norm;
        d8od = (d8 + sol(:,7))/norm;
        
        WT_frac_end(rr) = wtod(end)/(wtod(end)+d8od(end)); % wild type fraction at end
        yield(rr) = wtod(end) + d8od(end);
       
        

%         calculate relative EP production cost 


        k1 = param(1);
        k2 = param(3);
        gamma = param(5);
        Ecap = param(9);
        chi = param(13);
        gamma1 = param(23);
        k5 = param(12);
        W = sol(:,1); C = sol(:,2); A = sol(:,3); N = sol(:,4); E = sol(:,5); Ws = sol(:,6); Cs = sol(:,7); A0 = sol(:,8);
        g = k1*A.^2./(k2^2+A.^2);
        gA0 = k1*A0.^2./(k2^2+A0.^2);
        f = k5*(1-E/Ecap);
       
       
        growth_ind = find(EP_cost_col(rr)*f./(gamma1*gA0 + gamma*g)<1); %find indices where growth occurs
        t_growth = t(growth_ind);  % times where growth occurs
        tot_growth(rr) = trapz(t_growth,gamma1*gA0(growth_ind) + gamma*g(growth_ind)); % calculate total growth over time
        tot_cost(rr) = trapz(t_growth,chi*f(growth_ind)); % calculate total cost over time
        rel_cost_int(rr) = tot_cost(rr)/tot_growth(rr); % calculate relative cost as fraction of total cost vs total growth

    
    end    

    %% visualisation of relative cost vs final D8 pop
    figure(f2)
    hold on
    grid on
    p(ww) = plot(rel_cost_int,100 - 100*WT_frac_end, '--o', 'color', col(ww,:), 'Displayname', "Initial \Delta8 = "+ num2str(100 - 100*wt_frac) + "%");
    plot(rel_cost_int(EP_cost_col == 1),100 - 100*WT_frac_end(EP_cost_col == 1), 'o','color', col(ww,:), 'MarkerFaceColor',  'r')
    xlabel('avg. relative EP production cost, $C_{rel}$', 'interpreter', 'latex')
    ylabel('Final $\Delta 8$', 'interpreter', 'latex')
   
%     xlim([0,0.3])

    
    %% visualisation of total cost vs yield
    figure(f1);
    hold on
    grid on
    plot(tot_cost,yield, '--o')
    xlabel('avg. total EP production cost, $C_{tot}$', 'interpreter', 'latex')
    ylabel('Yield', 'interpreter', 'latex') 
    pbaspect([1 1 1])
    set(f1,'Windowstyle','normal')
    set(findall(f1,'-property','FontSize'),'FontSize',11)
    set(f1,'Units','centimeters')
    set(f1,'Position',[18 1 8 8])

    %% visualisation total cost vs chi
    figure(f3);
    hold on
    grid on
    plot(EP_cost_col, tot_cost, '--o')
    xlabel('EP production cost, $\chi$', 'interpreter', 'latex')
    ylabel('avg. total EP production cost, $C_{tot}$', 'interpreter', 'latex')
    pbaspect([1 1 1])
    set(f3,'Windowstyle','normal')
    set(findall(f3,'-property','FontSize'),'FontSize',11)
    set(f3,'Units','centimeters')
    set(f3,'Position',[18 1 8 8])

    %% visualisation total growth vs chi
    figure(f4);
    hold on
    grid on
    plot(EP_cost_col, tot_growth, '--o')
    xlabel('EP production cost, $\chi$', 'interpreter', 'latex')
    ylabel('avg. total growth of $\Delta 8$ (per unit), $G_{tot}$', 'interpreter', 'latex')
    pbaspect([1 1 1])
    set(f4,'Windowstyle','normal')
    set(findall(f4,'-property','FontSize'),'FontSize',11)
    set(f4,'Units','centimeters')
    set(f4,'Position',[18 1 8 8])

    %% visualisation of rel cost vs yield
    figure(f5);
    hold on
    grid on
    plot(rel_cost_int,yield, '--o')
    xlabel('avg. relative EP production cost, $C_{rel}$', 'interpreter', 'latex')
    ylabel('Yield', 'interpreter', 'latex') 
    pbaspect([1 1 1])
    set(f5,'Windowstyle','normal')
    set(findall(f5,'-property','FontSize'),'FontSize',11)
    set(f5,'Units','centimeters')
    set(f5,'Position',[18 1 8 8])

end

figure(f2)
leg = legend(p,'location','southeast');
% set(leg, 'Interpreter', 'latex')
set(f2,'Windowstyle','normal')
set(findall(f2,'-property','FontSize'),'FontSize',11)
set(f2,'Units','centimeters')
set(f2,'Position',[18 1 11.4 13])