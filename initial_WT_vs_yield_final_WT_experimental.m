%% Plot final yield and final WT % vs initial WT %

clear; close all

%% load data
data_msgg = readtable("Data_repo/CFU_coculture_MSgg_merge.csv");
data_msbg = readtable("Data_repo/CFU_coculture_MSbg_merge.csv");

% restrict msbg to 48h and exlcude heat treatments

data_msbg = data_msbg(data_msbg.Timepoint == 48,:);
data_msbg = data_msbg(data_msbg.Heat_treatment == 0,:);

title_col = ["MS-Ga+Gly", "MS-BSA+Gly"];
col = 1/255*[66, 150, 255; 133, 255, 23];

% %% Benford's law test
% firstDigit_msbg = testBenfordConformity(data_msbg.Count_Tot);

%% final CFU
data_msbg.CFU_tot = 10*data_msbg.Count_Tot./data_msbg.Dilution;
starting_wt_msbg = unique(data_msbg.x__NCIB3610); % find initial WT %
for ss = 1:length(starting_wt_msbg) % calc mean for each initial WT %
    msbg_mean(ss) = mean(data_msbg.CFU_tot(data_msbg.x__NCIB3610 == starting_wt_msbg(ss)));
end

data_msgg.CFU_tot = 10*data_msgg.Count_Tot./data_msgg.Dilution;
starting_wt_msgg = unique(data_msgg.x__NCIB3610); % find initial WT %
for ss = 1:length(starting_wt_msgg) % calc mean for each initial WT %
    msgg_mean(ss) = mean(data_msgg.CFU_tot(data_msgg.x__NCIB3610 == starting_wt_msgg(ss)));
end


f1 = figure;
scatter(100 - data_msgg.x__NCIB3610', data_msgg.CFU_tot',36, col(1,:),'LineWidth',2)
hold on
grid on
scatter(100 - data_msbg.x__NCIB3610', data_msbg.CFU_tot',36, col(2,:),'LineWidth',2)
scatter(100 - starting_wt_msgg, msgg_mean, 100, 'k', '_','LineWidth',3)
scatter(100 - starting_wt_msbg, msbg_mean, 100, 'k', '_','LineWidth',3)
legend(title_col, 'location', 'southeast')
set(gca,'yscale','log')


xlabel('Initial \Delta8 %')
ylabel('CFU/ml')
xticks([0,25,50,75,100])
set(f1,'Windowstyle','normal')
set(findall(f1,'-property','FontSize'),'FontSize',11)
set(f1,'Units','centimeters')
set(f1,'Position',[10 5 10 8.4])

%% final WT 
data_msgg.WTCFU_perc = 100*data_msgg.Count_NCIB3610./data_msgg.Count_Tot;
data_msbg.WTCFU_perc = 100*data_msbg.Count_NCIB3610./data_msbg.Count_Tot;

for ss = 1:length(starting_wt_msbg) % calc mean for each initial WT %
    msbg_mean_finalWT(ss) = mean(data_msbg.WTCFU_perc(data_msbg.x__NCIB3610 == starting_wt_msbg(ss)));
end

for ss = 1:length(starting_wt_msgg) % calc mean for each initial WT %
    msgg_mean_finalWT(ss) = mean(data_msgg.WTCFU_perc(data_msgg.x__NCIB3610 == starting_wt_msgg(ss)));
end

f2 = figure;
scatter(100 - data_msgg.x__NCIB3610, (100 - data_msgg.WTCFU_perc)./(100 - data_msgg.x__NCIB3610),36, col(1,:),'LineWidth',2)
hold on
grid on
scatter(100 - data_msbg.x__NCIB3610, (100 - data_msbg.WTCFU_perc)./(100 - data_msbg.x__NCIB3610),36, col(2,:),'LineWidth',2)
scatter(100 - starting_wt_msgg, (100 - msgg_mean_finalWT)./(100 - starting_wt_msgg'), 100, 'k', '_','LineWidth',3)
scatter(100 - starting_wt_msbg, (100 - msbg_mean_finalWT)./(100 - starting_wt_msbg'), 100, 'k', '_','LineWidth',3)
% legend(title_col, 'location', 'northwest')



xlabel('Initial \Delta8 %')
ylabel('\Delta8 relative fitness')
xticks([0,25,50,75,100])
set(f2,'Windowstyle','normal')
set(findall(f2,'-property','FontSize'),'FontSize',11)
set(f2,'Units','centimeters')
set(f2,'Position',[10 5 10 8.4])