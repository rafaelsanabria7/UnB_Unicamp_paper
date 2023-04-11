% Correlation analysis among parameter and peak loads
% rafael sanabria 01-07-2022

% clear memory
clear all; clc; close all;

peaks = load('peak_values.txt');
input = load('input.txt');  

fc = input(:,1);
ft = input(:,2);
Ec = input(:,3);
Gf = input(:,4);
K = input(:,5);

%  histogram(fc)
%  fig = gcf;
%  fig.PaperPositionMode = 'auto'
%  fig_pos = fig.PaperPosition;
%  fig.PaperSize = [fig_pos(3) fig_pos(4)];
%  print(fig,'fc','-dpdf')
%  
%  histogram(ft)
%  fig = gcf;
%  fig.PaperPositionMode = 'auto'
%  fig_pos = fig.PaperPosition;
%  fig.PaperSize = [fig_pos(3) fig_pos(4)];
%  print(fig,'ft','-dpdf')
%  
%  histogram(Ec)
%  fig = gcf;
%  fig.PaperPositionMode = 'auto'
%  fig_pos = fig.PaperPosition;
%  fig.PaperSize = [fig_pos(3) fig_pos(4)];
%  print(fig,'Ec','-dpdf')
%  
%  histogram(Gf)
%  fig = gcf;
%  fig.PaperPositionMode = 'auto'
%  fig_pos = fig.PaperPosition;
%  fig.PaperSize = [fig_pos(3) fig_pos(4)];
%  print(fig,'Gf','-dpdf')
%  
%  histogram(K)
%  fig = gcf;
%  fig.PaperPositionMode = 'auto'
%  fig_pos = fig.PaperPosition;
%  fig.PaperSize = [fig_pos(3) fig_pos(4)];
%  print(fig,'K','-dpdf')

fc_corr =  corr(fc,peaks,'Type','Spearman','Rows','complete');
ft_corr =  corr(ft,peaks,'Type','Spearman','Rows','complete');
Ec_corr =  corr(Ec,peaks,'Type','Spearman','Rows','complete');
Gf_corr =  corr(Gf,peaks,'Type','Spearman','Rows','complete');
K_corr  =  corr(K,peaks,'Type','Spearman','Rows','complete');

bar([fc_corr ft_corr Ec_corr Gf_corr K_corr])

 fig = gcf;
 fig.PaperPositionMode = 'auto'
 fig_pos = fig.PaperPosition;
 fig.PaperSize = [fig_pos(3) fig_pos(4)];
 print(fig,'L5_correlation','-dpdf')