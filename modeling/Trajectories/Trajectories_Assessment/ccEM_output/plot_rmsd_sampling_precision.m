% Clear workspace and close all figures
clear all;
close all;

% Import data from the txt file
ccEM = importdata('ccEM_calculations.txt');
ccEM=ccEM.data;

% Select most likely trajectory and make array
trj1_cc=[ccEM(1);ccEM(4);ccEM(7)];
% Array for time
time=[0;1;2];

% Define the colors
color1 = [217,73,49]/255;
color2 = [233,147,133]/255;

% Create the bar plot with the same pixels as for FoXS
fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
BAR = bar(time,trj1_cc, 'grouped', 'BarWidth', 0.5); % bar width 1, so bars are touching
hold on;
% Set colors of bars (based on the list from y_group)
BAR.FaceColor = color1;
% Customize the axes and labels
snapshots=["0 min";"1 min";"2 min"];
set(gca, 'XTickLabel', snapshots); % Set x-axis labels to snapshots
set(gca, 'FontSize',36,'FontName','Helvetica','Linewidth',3);
axis([-0.5 2.5 0 1]);
box off;
% Print the plot
print(fig,'ccEM_trj1.png','-dpng');