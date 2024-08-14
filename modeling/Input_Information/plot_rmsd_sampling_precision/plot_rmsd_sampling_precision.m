% Clear workspace and close all figures
clear all;
close all;

% Import data from the CSV file
data = importdata('rmsd_sampling_precision.csv');

% Extract relevant data from CSV files
snapshots = data.textdata(2:end); % skip first row for x-axis
sampling_precision = data.data(:, 1); % Second column
average_rmsd = data.data(:, 2); % Third column

% Create a matrix for bar plot
y_group = [sampling_precision, average_rmsd];

% Define the colors
color1 = [217,73,49]/255;
color2 = [233,147,133]/255;

% Create the bar plot with the same pixels as for FoXS
fig=figure('Renderer', 'painters', 'Position', [0 0 1000 525]);
BAR = bar(y_group, 'grouped', 'BarWidth', 1); % bar width 1, so bars are touching
hold on;

% Set colors of bars (based on the list from y_group)
BAR(1).FaceColor = color1;
BAR(2).FaceColor = color2;

% Customize the axes and labels
set(gca, 'XTickLabel', snapshots); % Set x-axis labels to snapshots
set(gca, 'FontSize',36,'FontName','Helvetica','Linewidth',3);

% Add legend without box
legend({' Precision', ' RMSD'}, 'Location', 'NorthWest','FontSize',36,'FontName','Helvetica','Linewidth',3);
legend boxoff

% Set limit for y-axis
ylim([0, 25]); 

% Print the plot
print(fig,'rmsd_sampling_precision.png','-dpng')