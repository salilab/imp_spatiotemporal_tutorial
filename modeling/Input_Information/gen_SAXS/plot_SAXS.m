clear all;
close all;

% Set colors ----------------------------------------------------

color1=[0, 0.4470, 0.7410];
color2=[0.8500, 0.3250, 0.0980];
color3=[0.9290, 0.6940, 0.1250];
color4=[0.4940, 0.1840, 0.5560];
color5=[0.4660, 0.6740, 0.1880];
color6=[0.3010, 0.7450, 0.9330];
color7=[0.6350, 0.0780, 0.1840];
color8=[160,  82,  45]/255;

blue1=[0,0,128]/255;
blue2=[65,105,225]/255;
blue3=[135,206,250]/255;

% Import SAXS data from files ---------------------------------------------

min0=importdata("0min_exp.dat",' ',2);
min0=min0.data;
min1=importdata("1min_exp.dat",' ',2);
min1=min1.data;
min2=importdata("2min_exp.dat",' ',2);
min2=min2.data;

% Plot SAXS data ----------------------------------------------------------

figure;
hold on;
plot(min0(:,1),min0(:,2),'Linewidth',4,'Color',color5);
plot(min1(:,1),min1(:,2),'Linewidth',4,'Color',color6);
plot(min2(:,1),min2(:,2),'Linewidth',4,'Color',color7);
errorbar(min0(:,1),min0(:,2),min0(:,3),'Linewidth',4,'Color',color5);
errorbar(min1(:,1),min1(:,2),min1(:,3),'Linewidth',4,'Color',color6);
errorbar(min2(:,1),min2(:,2),min2(:,3),'Linewidth',4,'Color',color7);
legend({' 0 min',' 1 min',' 2 min'},'location','northeast');
set(gca,'FontSize',56,'FontName','Helvetica','Linewidth',4,'YScale', 'log');
% axis([-1 2.5 0 1])
% box on;
legend boxoff;

