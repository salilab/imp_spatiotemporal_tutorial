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

% ChiSquare ----------------------------------------------------------

Chi=importdata('snapshot_1_2min.ChiSquare_Grid_Stats.txt');
prec=importdata('snapshot_1_2min.Sampling_Precision_Stats.txt',' ',3);
prec=cell2mat(prec(3));
prec=str2num(prec);

x=[0;max(Chi(:,1))];
y_p=[0.05;0.05];
y_V=[0.1;0.1];
y_pop=[0.8;0.8];

y2=[0;1];
x2=[prec(1);prec(1)];

fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
plot(Chi(:,1),Chi(:,2),'o','MarkerFaceColor',color5,'MarkerSize',8,'Linewidth',3,'Color',color5);
plot(Chi(:,1),Chi(:,3),'o','MarkerFaceColor',color3,'MarkerSize',8,'Linewidth',3,'Color',color3);
plot(Chi(:,1),Chi(:,4)./100,'o','MarkerFaceColor',color4,'MarkerSize',8,'Linewidth',3,'Color',color4);
plot(x,y_p,'--','MarkerFaceColor',color5,'MarkerSize',8,'Linewidth',3,'Color',color5);
plot(x,y_V,'--','MarkerFaceColor',color3,'MarkerSize',8,'Linewidth',3,'Color',color3);
plot(x,y_pop,'--','MarkerFaceColor',color4,'MarkerSize',8,'Linewidth',3,'Color',color4);
plot(x2,y2,'--k','MarkerSize',8,'Linewidth',3);
legend({' \it\chi^2\rm p-value',' Cramer''s \itV\rm',' % clustered'},'location','west');
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3);
axis([0 24 0 1]);
legend boxoff;
print(fig,'ChiSquare_v2.png','-dpng');

% Top score ----------------------------------------------------------

TopScore=importdata('snapshot_1_2min.Top_Score_Conv.txt');

y=[TopScore(1,2);TopScore(1,2)];
x=[0;13000];

fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
plot(TopScore(:,1),TopScore(:,2),'o','MarkerFaceColor',color7,'MarkerSize',8,'Linewidth',3,'Color',color7);
plot(x,y,'--','MarkerFaceColor',color6,'MarkerSize',8,'Linewidth',3,'Color',color6);
errorbar(TopScore(:,1),TopScore(:,2),TopScore(:,3),'o','MarkerFaceColor',color7,'MarkerSize',8,'Linewidth',3,'Color',color7);
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3);
% legend boxoff;
axis([0 13000 240 250]);
print(fig,'TopScore_v2.png','-dpng');

% Score histogram ----------------------------------------------------------

A=importdata('snapshot_1_2min.Score_Hist_A.txt');
B=importdata('snapshot_1_2min.Score_Hist_B.txt');

% Fill in 0 for preceeding and trailing entries in the histogram
N=max(size(A));
N2=max(size(B));
A2=zeros(N+2,2);
B2=zeros(N2+2,2);
dRA=A(2,1)-A(1,1);
A2(1,1)=A(1)-dRA;
dRB=B(2,1)-B(1,1);
B2(1,1)=B(1)-dRB;
A2(2:N+1,:)=A;
B2(2:N2+1,:)=B;
A2(N+2,1)=A(N)+dRA;
B2(N2+2,1)=B(N2)+dRB;



fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
plot(A2(:,1),A2(:,2),'MarkerFaceColor',color1,'MarkerSize',8,'Linewidth',3,'Color',color1);
plot(B2(:,1),B2(:,2),'MarkerFaceColor',color2,'MarkerSize',8,'Linewidth',3,'Color',color2);
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3);
legend({' Sample 1',' Sample 2'},'location','northeast');
legend boxoff;
axis([240 350 0 300]);
print(fig,'Score_dist_v2.png','-dpng');

% Cluster population bar graph --------------------------------------------

Pop=importdata('snapshot_1_2min.Cluster_Population.txt');

width=0.2;

fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
hold on;
bar(Pop(:,1)-width,Pop(:,2),2*width,'Linewidth',3,'FaceColor',color1,'EdgeColor','none');
bar(Pop(:,1)+width,Pop(:,3),2*width,'Linewidth',3,'FaceColor',color2,'EdgeColor','none');
set(gca,'FontSize',36,'FontName','Helvetica','Linewidth',3);
legend({' Sample 1',' Sample 2'},'location','northeast');
legend boxoff;
axis([-0.5 2.5 0 8000]);
xticks([0; 1; 2]);
xticklabels(['Cluster 0';'Cluster 1';'Cluster 2'])
print(fig,'Cluster_pop_v2.png','-dpng');
