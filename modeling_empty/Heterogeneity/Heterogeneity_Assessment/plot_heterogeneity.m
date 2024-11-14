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

grey=[188,  188,  188]/255;

blue1=[0,0,128]/255;
blue2=[65,105,225]/255;
blue3=[135,206,250]/255;

% Import experimental data from files -------------------------------------

expA=importdata("../../Input_Information/gen_FCS/A_v2.txt");
expA=expA.data;
expB=importdata("../../Input_Information/gen_FCS/B_v2.txt");
expB=expB.data;
expC=importdata("../../Input_Information/gen_FCS/C_v2.txt");
expC=expC.data;

% Import forward copy numbers as a function of time -----------------------

model0=importdata("../Heterogeneity_Modeling/0min.txt",' ',1);
model0=model0.data;
model1=importdata("../Heterogeneity_Modeling/1min.txt",' ',1);
model1=model1.data;
model2=importdata("../Heterogeneity_Modeling/2min.txt",' ',1);
model2=model2.data;

model_time=[0,1,2];

% Convert copy number as a function of time to copy number as a function of 
% protein -----------------------------------------------------------------

modelA=zeros(7,2);
modelA(1:3,1)=model0(:,1);
modelA(1:3,2)=0;
modelA(4:6,1)=model1(:,1);
modelA(4:6,2)=1;
modelA(7,1)=model2(:,1);
modelA(7,2)=2;

modelB=zeros(7,2);
modelB(1:3,1)=model0(:,2);
modelB(1:3,2)=0;
modelB(4:6,1)=model1(:,2);
modelB(4:6,2)=1;
modelB(7,1)=model2(:,2);
modelB(7,2)=2;

modelC=zeros(7,2);
modelC(1:3,1)=model0(:,3);
modelC(1:3,2)=0;
modelC(4:6,1)=model1(:,3);
modelC(4:6,2)=1;
modelC(7,1)=model2(:,3);
modelC(7,2)=2;


% Plot A ---------------------------------------

fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
plot(expA(:,1),expA(:,2),'Linewidth',3,'Color',grey);
plot(modelA(:,2),modelA(:,1),'o','MarkerFaceColor',blue1,'MarkerSize',16,'Linewidth',3,'Color',blue1);
errorbar(expA(:,1),expA(:,2),expA(:,3),'Linewidth',3,'Color',grey);
legend({' experiment',' model'},'location','north');
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3);
axis([-0.1 2.6 -1 3])
legend boxoff;
print(fig,'heterogeneity_compA.png','-dpng');



% Plot B ---------------------------------------

fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
plot(expB(:,1),expB(:,2),'Linewidth',3,'Color',grey);
plot(modelB(:,2),modelB(:,1),'o','MarkerFaceColor',blue1,'MarkerSize',16,'Linewidth',3,'Color',blue1);
errorbar(expB(:,1),expB(:,2),expB(:,3),'Linewidth',3,'Color',grey);
legend({' experiment',' model'},'location','north');
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3);
axis([-0.1 2.6 -1 3])
legend boxoff;
print(fig,'heterogeneity_compB.png','-dpng');

% Plot C ---------------------------------------

fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
plot(expC(:,1),expC(:,2),'Linewidth',3,'Color',grey);
plot(modelC(:,2),modelC(:,1),'o','MarkerFaceColor',blue1,'MarkerSize',16,'Linewidth',3,'Color',blue1);
errorbar(expC(:,1),expC(:,2),expC(:,3),'Linewidth',3,'Color',grey);
legend({' experiment',' model'},'location','north');
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3);
axis([-0.1 2.6 -1 3])
legend boxoff;
print(fig,'heterogeneity_compC.png','-dpng');

