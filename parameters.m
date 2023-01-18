%% Parameter values for the model

param(3) = 0.02; % growth half saturation constant
param(1) = 3; %growth rate in response to A
param(18) = 0.4; %max sporulation rate
param(19) = 0.003; %sporulation half saturation constant
param(11) = 0.5; % max nutrient conversion rate
param(9) = 0.2; % carrying capacity of E
param(13) = 1; % EP production cost
param(12) = 0.03; %max EP production rate
param(5) = 4; %gamma
param(23) = 1; %gamma1
param(14) = 0.45; % nutrient conversion half saturation constant

norm = 2.01; % normalisation constant for visualisations (stationary phase in WT MSgg) - could automate by running the model here once.