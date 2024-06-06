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

% Import sampled data from files ----------------------------------------------------

expA=importdata("A_v2.txt");
expA=expA.data;
expB=importdata("B_v2.txt");
expB=expB.data;
expC=importdata("C_v2.txt");
expC=expC.data;

time=[-1:0.1:2.5];

% Compute true distribution ------------------------------------------------

eta=5;
shift=-0.2;
CNA=0.5+0.5*tanh(eta*(time-(2+shift)));
CNB=0.5+0.5*tanh(eta*(time-(shift)));
CNC=0.5+0.5*tanh(eta*(time-(1+shift)));

% Plot true distribution ---------------------------------

figure;
hold on;
plot(time,CNA,'Linewidth',4,'Color',color1);
plot(time,CNB,'Linewidth',4,'Color',color2);
plot(time,CNC,'Linewidth',4,'Color',color4);
legend({' A',' B',' C'},'location','northwest');
set(gca,'FontSize',56,'FontName','Helvetica','Linewidth',4);
axis([-1 2.5 0 1])
% box on;
legend boxoff;

% Plot sampled data ---------------------------------------

figure;
hold on;
plot(expA(:,1),expA(:,2),'Linewidth',4,'Color',color1);
plot(expB(:,1),expB(:,2),'Linewidth',4,'Color',color2);
plot(expC(:,1),expC(:,2),'Linewidth',4,'Color',color4);
errorbar(expA(:,1),expA(:,2),expA(:,3),'Linewidth',4,'Color',color1);
errorbar(expB(:,1),expB(:,2),expB(:,3),'Linewidth',4,'Color',color2);
errorbar(expC(:,1),expC(:,2),expC(:,3),'Linewidth',4,'Color',color4);
legend({' A',' B',' C'},'location','north');
set(gca,'FontSize',56,'FontName','Helvetica','Linewidth',4);
axis([-0.1 2.6 -1 3])
legend boxoff;

