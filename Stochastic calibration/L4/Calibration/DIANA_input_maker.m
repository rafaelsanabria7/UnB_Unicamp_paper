% input data (.DAT) for DIANA 
% Rafael sanabria 2021

clc; clear all; close all;
tic
%LHS variables for fc, fct, E, Gf
LHS_sample = load('input.txt');

for j=1:length(LHS_sample) 
fc = LHS_sample(:,1);
fct = LHS_sample(:,2);
Ec = LHS_sample(:,3);
Gf = LHS_sample(:,4);
K = LHS_sample(:,5);

% open .dat file and export information
file_input = strcat('L4','.dat');  

% Open model
fid = fopen(file_input);
i = 1;

tline = fgetl(fid);
A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
    end

% capture specific lines
specific_lines = [35331 35336 35337 35341 35342 35399];

% replacing values - units [N mm]
A{specific_lines(1)} = cell2mat(strcat({'     YOUNG     '},num2str(Ec(j))));
A{specific_lines(2)} = cell2mat(strcat({'     TENSTR    '},num2str(fct(j))));
A{specific_lines(3)} = cell2mat(strcat({'     GF1   '},num2str(Gf(j)/1000)));
A{specific_lines(4)} = cell2mat(strcat({'     COMSTR    '},num2str(fc(j))));
A{specific_lines(5)} = cell2mat(strcat({'     GC    '},num2str(Gf(j)/1000*250)));
A{specific_lines(6)} = cell2mat(strcat({'     SPRING    '},num2str(K(j))));

% print new .datfile
inputname = strcat('input_',num2str(j));
fid = fopen([inputname,'.dat'], 'w');
for i=1:size(A,2)
fprintf(fid,'%s \n', A{i});
end
fclose('all');

%% Commands and output
% open .dcf file and export information
file_input = strcat('commands','.dcf');  

% Open model
fid = fopen(file_input);
i = 1;

tline = fgetl(fid);
B{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        B{i} = tline;
    end

% capture specific lines
specific_output = [109 131 279 301 373 395];

% replacing output names
B{specific_output(1)} = cell2mat(strcat({'    FILE  L4_'},num2str(j),'_FORCE'));
B{specific_output(2)} = cell2mat(strcat({'    FILE  L4_'},num2str(j),'_DISP'));
B{specific_output(3)} = cell2mat(strcat({'    FILE  L4_'},num2str(j),'_FORCE'));
B{specific_output(4)} = cell2mat(strcat({'    FILE  L4_'},num2str(j),'_DISP'));
B{specific_output(5)} = cell2mat(strcat({'    FILE  L4_'},num2str(j),'_FORCE'));
B{specific_output(6)} = cell2mat(strcat({'    FILE  L4_'},num2str(j),'_DISP'));


% print new .dcf file
inputname = strcat('command_',num2str(j));
fid = fopen([inputname,'.dcf'], 'w');
for i=1:size(B,2)-1
fprintf(fid,'%s \n', B{i});
end
fclose('all'); 
end

%% Runner
% Load step
k = 0;
model_path = 'C:/Users/labmem/Desktop/Rafael/UnB/';
for j=1:length(LHS_sample) 
C{1 + k} = cell2mat(strcat({'importModel( "'},model_path,'input_',num2str(j),{'.dat" )'}));
C{2 + k} = cell2mat(strcat({'addAnalysis( "command_'},num2str(j),{'" )'}));
C{3 + k} = cell2mat(strcat({'loadAnalysisCommands( "command_'},num2str(j),{'","'},model_path,'command_',num2str(j),{'.dcf" )'}));
C{4 + k} = cell2mat(strcat({'runSolver( "command_'},num2str(j),{'" )'}));
C{5 + k} = cell2mat({' '});
k = k + 5;
end

% print new .dcf file
inputname = strcat('batch_runner');
fid = fopen([inputname,'.txt'], 'w');
for i=1:size(C,2)
fprintf(fid,'%s \n', C{i});
end
fclose('all'); 


toc