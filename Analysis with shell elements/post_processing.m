% Post-processing of composed lines elements in DIANA for LOA III
clc; clear all; close all;
file_tab = strcat('MEdX','.tb');  


bsrx  = 477;

% Open model
fid = fopen(file_tab);
i = 1;

tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end

k=1;
for i=11:size(A,2) - 4
   B = regexp(A{i}, '\s', 'split');
   aux = B(~cellfun('isempty',B));
   Y_coordinate{k} = num2str(aux{end-1});
   MEdX{k} = num2str(aux{end-3});
   C{i}  = aux;
   k = k+1;
end
fclose(fid);

[Y_coordinate_sorted,I] = sort(str2double(Y_coordinate));
MEdX_sorted = str2double(MEdX(I));

hold on
plot(Y_coordinate_sorted, MEdX_sorted/180,'.-')
plot([bsrx/2 bsrx/2],[0 7E4])
plot([-bsrx/2 -bsrx/2],[0 7E4])

indexesInRange = Y_coordinate_sorted > -bsrx/2  & Y_coordinate_sorted < bsrx/2;
subVector = Y_coordinate_sorted(indexesInRange);

Msdx_ave = trapz(Y_coordinate_sorted(indexesInRange), MEdX_sorted(indexesInRange))/bsrx/180


%% Analysis for the other direction
file_tab = strcat('MEdY','.tb');  
bsry  =388.65;

% Open model
fid = fopen(file_tab);
i = 1;

tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end

k=1;
for i=11:size(A,2) - 4
   B = regexp(A{i}, '\s', 'split');
   aux = B(~cellfun('isempty',B));
   X_coordinate{k} = num2str(aux{end-2});
   MEdY{k} = num2str(aux{end-3});
   C{i}  = aux;
   k = k+1;
end
fclose(fid);

[X_coordinate_sorted,I] = sort(str2double(X_coordinate));
MEdY_sorted = str2double(MEdY(I));

figure
hold on
plot(X_coordinate_sorted, MEdY_sorted/180,'.-')
plot([150 150],[0 7E4])
plot([bsry bsry],[0 7E4])
% 
indexesInRange = X_coordinate_sorted > 150  & X_coordinate_sorted < bsry;
subVector = X_coordinate_sorted(indexesInRange);
% 
Msdy_ave = trapz(X_coordinate_sorted(indexesInRange), MEdY_sorted(indexesInRange))/bsry/180
