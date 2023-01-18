function dydt = odesys(t,y,param)

%% function called during ode solver
% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 23/11/2022

k1 = param(1);
k2 = param(3);
gamma = param(5);
Ecap = param(9);
k6 = param(11);
k5 = param(12);
chi = param(13);
k7 = param(14);
k3  = param(18);
k4 = param(19);
gamma1 = param(23);

W = y(1); C = y(2); A = y(3); N = y(4); E = y(5); Ws = y(6); Cs = y(7); A0 = y(8);
dydt = zeros(4,1);

g = k1*A^2/(k2^2+A^2);
gA0 = k1*A0^2/(k2^2+A0^2);
f = k5*(1-E/Ecap);
h = k6*N/(k7+N)*E;
d = k3*(1-(A+A0)^2/(k4^2 + (A+A0)^2));

dydt(1) = max(0,(gamma*g+gamma1*gA0 - chi*f))*W - d*W; % W eq
dydt(2) = (gamma*g+gamma1*gA0)*C - d*C; % C eq
dydt(3) =  - g*(W+C); % A eq
dydt(4) =  -h; % N eq
dydt(5) = f*W; % E eq
dydt(6) = d*W; % Ws eq
dydt(7) = d*C; % Cs eq
dydt(8) = h - gA0*(W+C); %A0 eq