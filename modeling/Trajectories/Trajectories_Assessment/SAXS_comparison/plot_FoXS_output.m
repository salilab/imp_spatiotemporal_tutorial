clear all;
close all;

% Set colors ----------------------------------------------------

grey=[188,  188,  188]/255;

blue1=[0,0,128]/255;
blue2=[65,105,225]/255;
blue3=[135,206,250]/255;

% Import SAXS data from files ---------------------------------------------
% 0 min
SS1_0=importdata("snapshot1_0min_0min_exp.fit",' ',3);
SS1_0=SS1_0.data;
SS2_0=importdata("snapshot2_0min_0min_exp.fit",' ',3);
SS2_0=SS2_0.data;
SS3_0=importdata("snapshot3_0min_0min_exp.fit",' ',3);
SS3_0=SS3_0.data;

% 1 min
SS1_1=importdata("snapshot1_1min_1min_exp.fit",' ',3);
SS1_1=SS1_1.data;
SS2_1=importdata("snapshot2_1min_1min_exp.fit",' ',3);
SS2_1=SS2_1.data;
SS3_1=importdata("snapshot3_1min_1min_exp.fit",' ',3);
SS3_1=SS3_1.data;

% 2min
SS1_2=importdata("snapshot1_2min_2min_exp.fit",' ',3);
SS1_2=SS1_2.data;

% Plot SAXS data ----------------------------------------------------------
% 0 min
fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
plot(SS1_0(:,1),SS1_0(:,2),'Linewidth',3,'Color',grey);
plot(SS1_0(:,1),SS1_0(:,4),'Linewidth',3,'Color',blue1);
plot(SS2_0(:,1),SS2_0(:,4),'Linewidth',3,'Color',blue2);
plot(SS3_0(:,1),SS3_0(:,4),'Linewidth',3,'Color',blue3);
errorbar(SS1_0(:,1),SS1_0(:,2),SS1_0(:,3),'Linewidth',3,'Color',grey);
plot(SS1_0(:,1),SS1_0(:,4),'Linewidth',3,'Color',blue1);
plot(SS2_0(:,1),SS2_0(:,4),'Linewidth',3,'Color',blue2);
plot(SS3_0(:,1),SS3_0(:,4),'Linewidth',3,'Color',blue3);
legend({' experiment',' state 1',' state 2', ' state 3'},'location','northeast');
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3,'YScale', 'log');
axis([0 0.5 10^4 3*10^6]);
xticks([0;0.1;0.2;0.3;0.4;0.5]);
yticks([10^4;10^5;10^6]);
yticklabels(['10^4';'10^5';'10^6']);
% box on;
legend boxoff;
print(fig,'0min_FoXS.png','-dpng');

% 1 min
fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
plot(SS1_1(:,1),SS1_1(:,2),'Linewidth',3,'Color',grey);
plot(SS1_1(:,1),SS1_1(:,4),'Linewidth',3,'Color',blue1);
plot(SS2_1(:,1),SS2_1(:,4),'Linewidth',3,'Color',blue2);
plot(SS3_1(:,1),SS3_1(:,4),'Linewidth',3,'Color',blue3);
errorbar(SS1_1(:,1),SS1_1(:,2),SS1_1(:,3),'Linewidth',3,'Color',grey);
plot(SS1_1(:,1),SS1_1(:,4),'Linewidth',3,'Color',blue1);
plot(SS2_1(:,1),SS2_1(:,4),'Linewidth',3,'Color',blue2);
plot(SS3_1(:,1),SS3_1(:,4),'Linewidth',3,'Color',blue3);
legend({' experiment',' state 1',' state 2', ' state 3'},'location','northeast');
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3,'YScale', 'log');
axis([0 0.5 10^4 10^7]);
xticks([0;0.1;0.2;0.3;0.4;0.5]);
yticks([10^4;10^5;10^6;10^7]);
yticklabels(['10^4';'10^5';'10^6';'10^7']);
% box on;
legend boxoff;
print(fig,'1min_FoXS.png','-dpng');

% 2 min
fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
plot(SS1_2(:,1),SS1_2(:,2),'Linewidth',3,'Color',grey);
plot(SS1_2(:,1),SS1_2(:,4),'Linewidth',3,'Color',blue1);
errorbar(SS1_2(:,1),SS1_2(:,2),SS1_2(:,3),'Linewidth',3,'Color',grey);
plot(SS1_2(:,1),SS1_2(:,4),'Linewidth',3,'Color',blue1);
legend({' experiment',' state 1'},'location','northeast');
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3,'YScale', 'log');
axis([0 0.5 10^4 4*10^7]);
xticks([0;0.1;0.2;0.3;0.4;0.5]);
yticks([10^4;10^5;10^6;10^7]);
yticklabels(['10^4';'10^5';'10^6';'10^7']);
% box on;
legend boxoff;
print(fig,'2min_FoXS.png','-dpng');

 