% Post-processing data results from .tb document in real-time
% rafael sanabria 2021

% clear memory
clear all; clc; close all;

data = load('LA5_NON_EC2.txt');

num_colours = 100;
colourmap = jet(100);
weights = linspace(0,1,num_colours);

coordinates = [data(:,1),data(:,2)];
shear_dist = -data(:,3);

figure
num_points = size(data,1);

min_shear = -5;
max_shear = 5;

    hold on
    column_points_1 = [0 150  0];
    column_points_2 = [300 150  0];
    column_points_3 = [300 -150   0];
    column_points_4 = [0 -150 0];

    xyz = vertcat(column_points_1,column_points_2,column_points_3,column_points_4,column_points_1);
    plot3(xyz(:,1),xyz(:,2),xyz(:,3),'k');
    
    plot3(coordinates(:,1),coordinates(:,2),coordinates(:,2)*0,'k')

for idx = 1:num_points
   
    w = (shear_dist(idx) - min_shear)/(max_shear - min_shear);%(shear_dist(idx) - min(shear_dist))/(max(shear_dist) - min(shear_dist));  
    [~,ind] = min(abs(weights-w));
    color = colourmap(ind,:);
    stem3(coordinates(idx,1),coordinates(idx,2),shear_dist(idx),'Color',color,'MarkerEdgeColor', color, 'MarkerFaceColor', color)
    zlim([min_shear max_shear])
    zticks([-5 -2.5 0  2.5  5]);
 
    %surf(xdata, ydata, 217*ones(2), 'FaceAlpha', 0.1, 'FaceColor', 'blue');
end

 cmap = colormap(jet(100)) ; %Create Colormap
 cbh = colorbar ; %Create Colorbar
 cbh.Ticks = linspace(0, 1, 11) ; %Create 8 ticks from zero to 1
 cbh.TickLabels = num2cell(min_shear:(max_shear-min_shear)/(10):max_shear); 
 view([115 33])
 grid on

 xlabel('x [mm]') 
 ylabel('y [mm]') 
 zlabel('Avergage shear stress[MPa]') 

 fig = gcf;
 fig.PaperPositionMode = 'auto'
 fig_pos = fig.PaperPosition;
 fig.PaperSize = [fig_pos(3) fig_pos(4)];
 print(fig,'LA5_NON_EC2','-dpdf')
 
