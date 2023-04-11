% Post-processing data results from .tb document in real-time
% rafael sanabria 2021
% clear memory
clear all; clc; close all;

% Load experimental curve
L1_exp = load('L1_LVDT.csv');

for step = 1:30
output = strcat('L1_',num2str(step),'_FORCE','.tb'); 
% Open model
fid = fopen(output);

        i = 1;
        j = 1;
        k = 1;
        flag=0;
        
       while flag==0
            readline=fgetl(fid);
            if readline==-1
                flag=1;
            end
            lineloc=strfind(readline,'7496');
            if isempty(lineloc)==0
                %stores line number
                srcline(j,1)=i;
                line = regexp(readline,'\s','split');
                aux = line(~cellfun('isempty',line));
                if strcmp(aux{:,1},'7496') == 1
                disp_1(k) =  str2double(aux{2});
                k = k + 1;
                end
            end
            i=i+1;
       end
 fclose(fid);
  
output = strcat('L1_',num2str(step),'_DISP','.tb');  
fid = fopen(output);      
        i = 1;
        j = 1;
        k = 1;
        flag=0;
        
       while flag==0 
            readline=fgetl(fid);
            if readline==-1
                flag=1;
            end
            lineloc=strfind(readline,'39185');
            if isempty(lineloc)==0
                %stores line number
                srcline(j,1)=i;
                line = regexp(readline,'\s','split');
                aux = line(~cellfun('isempty',line));
                if strcmp(aux{:,1},'39185') == 1
                   force(k) =  str2double(aux{2});
                   k = k + 1;
                end
               
            end
            i=i+1;
       end
 fclose(fid);

%  output_3 = strcat('N0_10_',num2str(step),'_DISP','.tb');  
%  fid = fopen(output_3);      
%         i = 1;
%         j = 1;
%         k = 1;
%         flag=0;
%         
%        while flag==0 
%             readline=fgetl(fid);
%             if readline==-1
%                 flag=1;
%             end
%             lineloc=strfind(readline,'383');
%             if isempty(lineloc)==0
%                 %stores line number
%                 srcline(j,1)=i;
%                 line = regexp(readline,'\s','split');
%                 aux = line(~cellfun('isempty',line));
%                 if strcmp(aux{:,1},'383') == 1
%                    disp_1(k) =  str2double(aux{2});
%                    k = k + 1;
%                 end
%                
%             end
%             i=i+1;
%        end
%  fclose(fid);
 
 
A = [abs(disp_1 - disp_1(1)); abs(force)*2/1000]'; % store each numerical response
 
% Load at specific displacement values
u_disp = [1 2 3 4.5 max(abs(L1_exp(:,5)))];
[new_disp index] = unique(abs(L1_exp(:,5)));

test_force = interp1(new_disp,abs(L1_exp(index,1)),u_disp);

% peaks loads
[peaks(step), indexmax] = max(abs(force)*2/1000);

% interpolation of P values for the specififc displacement:
[uniq_A, index] = unique(A(1:indexmax,1));   
P_values(step,:) = interp1(uniq_A, A(index,2),u_disp);
P_values(isnan(P_values))=0; % delete NaN results

hold on
if step == 28
plot(abs(disp_1- disp_1(1)), abs(force)*2/1000,'-','Color',[0 0 1],'LineWidth',1.8)  
else
plot(abs(disp_1- disp_1(1)), abs(force)*2/1000,'-','Color',[0.6010 0.7450 0.9330])
end

%plot(u_disp, P_values(step,:),'r.')
title('L1 Machine Learning')
ylim([0 700])
xlim([0 20])
xlabel('Displacement [mm]') 
ylabel('Load [kN]') 
% 
peaks(step) = max(abs(force)*2/1000);
convergence(step) = sum(peaks)./step;
% plot(convergence,'r.-')

disp_1 = [];
force = [];
end
plot(-L1_exp(:,5), L1_exp(:,1)+25,'r-.','LineWidth',1.5)
% plot(-L1_exp(:,5), L1_exp(:,1)+25,'k-.')
% plot([0 10], [437 437],'r-')

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'L1','-dpdf')
 
% 
% N = [1:9];
% legendStrings = "N = " + string(N);
% legend(legendStrings)

%Histogram
clear line
hold on
figure
P_exp = 437;
histogram(peaks)
hax=gca;
line([P_exp P_exp],get(hax,'YLim'),'LineWidth',2,'Color','r')
line([mean(peaks) mean(peaks)],get(hax,'YLim'),'LineWidth',2,'Color','b')
hold off;

[M, I] = min(abs(mean(peaks) - peaks));
best_curve = I;

%print new .txt file
save('load_values.txt', 'P_values', '-ascii', '-tabs')
peaks = peaks';
save('peak_values.txt', 'peaks', '-ascii', '-tabs')