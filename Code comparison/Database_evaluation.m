%% Database evaluation of flat slab-edge column connections 
% According to current EC2, fib Model Code 2010, FprEC2 (2022)
% code developed by: Rafael Diaz - University of Campinas

% Paper: Flat slab-edge column connections subjected to outward and inward 
% eccentricity: numerical modeling and comparison with codes
% DOI: (To be defined)

%% License 
%MIT License
%Copyright (c) [2023] [Rafael A S Diaz]

%Permission is hereby granted, free of charge, to any person obtaining a copy
%of this software and associated documentation files, to deal
%in the Software without restriction, including without limitation the rights
%to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%copies of the Software, and to permit persons to whom the Software is
%furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all
%copies or substantial portions of the Software.

%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%SOFTWARE.

% clear memory
clear all; clc; close all;

% load data
data =  readtable('Database.csv');

%% Input parameters
e = table2array(data(:,2));  % eccentricity
c1 = table2array(data(:,3)); % column size parallel to the free edge
c2 = table2array(data(:,4)); % column size perpendicular to the free edge
d  = table2array(data(:,5)); % slab effective depth
rho = table2array(data(:,6));% reinforcement ratio
fc = table2array(data(:,7)) - 4; % characteristic concrete strength
fy = table2array(data(:,8)); % steel yield strength
V_exp = table2array(data(:,9)); % experimental punching capacity
dg = table2array(data(:,10)); % aggregate size
rs_x = table2array(data(:,13)); % moment contraflexure radius in x-direction
rs_y = table2array(data(:,14)); % moment contraflexure radius in y-direction

tests_i = 1;
tests_j = 86;

%% Eurocode 2 (EN 1992-1-1:2004)
for i = tests_i:tests_j 
   
   k(i) = sqrt(200/d(i))+1; % size effect coefficient

   % for original values of C1 and C2
   if c1(i) == c2(i)
       K(i) = 0.6;
   else
     if c1(i)/c2(i) < 0.5
       K(i) = 0.45;
     else
       if c1(i)/c2(i) > 3.0
       K(i) = 0.8;
       else
       K(i) = interp1([0.5 1 2 3],[0.45 0.6 0.7 0.8],c1(i)/c2(i));
       end
     end
   end
   
% Two approaches considered: 
% (i) Fully plasitc shear stress distribution 
% (ii) Reduced perimeter

   u1(i) = 2*c1(i) + c2(i) + 2*d(i)*pi;
   u1_reduced(i) = c2(i) + 2*d(i)*pi + 2*min(1.5*d(i),0.5*c1(i));
   centroid(i) = (2*c1(i)*c1(i)/2 + 2*(2*d(i)*pi/2)*((c1(i)*2*d(i)*pi/2+(2*d(i)).^2)/(2*d(i)*pi/2)) +   c1(i)*(c1(i) + 2*d(i)))/u1(i);
   Wp1(i) = (c1(i) - centroid(i)).^2 + centroid(i).^2 + c1(i)*c2(i) + 2.*d(i)*c2(i) - c2(i).*centroid(i) + 8*d(i).^2  + 2*pi.*d(i).*(c1(i) - centroid(i));
   
   beta(i) = 1 + K(i)*abs(-e(i) - centroid(i) + c1(i)/2)*u1(i)/Wp1(i);
   beta_reduced(i) = u1(i)/u1_reduced(i);
   Mflex(i) = min(rho(i)/100*d(i)^2*fy(i)*(1 - rho(i)/100*fy(i)/(2*fc(i)))*(2*c1(i)+c2(i)),0.25*d(i)^2*fc(i)*(2*c1(i)+c2(i)));
   
   VEC2(i) = 0.18*k(i)*(rho(i)*fc(i))^(1/3)*d(i)*u1(i)/1000/beta(i);
   VEC2_reduced(i) = 0.18*k(i)*(rho(i)*fc(i))^(1/3)*d(i)*u1(i)/1000/beta_reduced(i);   
   ratio_moments_EC2(i) = V_exp(i)*abs(-e(i) - centroid(i) + c1(i)/2)/Mflex(i)*1E3;
end

ratio_EC2 = V_exp(tests_i:tests_j)./VEC2(tests_i:tests_j)';
Average_EC2 = mean(V_exp(tests_i:tests_j)./VEC2(tests_i:tests_j)')
CoV_EC2 = std(V_exp(tests_i:tests_j)./VEC2(tests_i:tests_j)')/Average_EC2
Min_EC2 = min(V_exp(tests_i:tests_j)./VEC2(tests_i:tests_j)')

ratio_EC2_reduced = V_exp(e<=0)./VEC2_reduced(e<=0)';
Average_EC2_reduced = mean(ratio_EC2_reduced)
CoV_EC2_reduced = std(V_exp(e<=0)./VEC2_reduced(e<=0)')/Average_EC2_reduced
Min_EC2_reduced = min(V_exp(e<=0)./VEC2_reduced(e<=0)')

%% fib Model Code 2010 (2013) (LoA II)
for i = tests_i:tests_j 

if dg(i) ~= 0 % Apply only with slabs with size aggregate information
    
% variables related with the eccentricity
Ac(i) = c1(i).*c2(i) + 2*c1(i).*d(i)/2 + c2(i).*d(i)/2 + d(i).^2/8*pi;
bu(i) = sqrt(4/pi*Ac(i));
kdg(i) = max(0.75,32/(16+dg(i)));

delta_ex(i) = 1/4*(2*c1(i)^2 + 6*c1(i)*d(i) + 3*d(i)^2)/(3*c1(i) + 2*d(i));
delta_ey(i) = 0;

eu(i) = abs(-e(i) - delta_ex(i)); 
ke(i) = 1/(1+ eu(i)/bu(i));
bo(i) = ke(i)*(2*c1(i) + c2(i) + d(i)*pi/2);

% variables related with moment capacity in the strip
bs = sqrt(rs_x(i)*rs_y(i))*1.5;
bsr_x(i) = min([3*c1(i) bs]);
bsr_y(i) =  min([c2(i)/2 + bs/2 bs]);

Vd = 0:1:1000000;
for j = 1:size(Vd,2)
mdx(j) = Vd(j)/8 + abs(-Vd(j)*e(i) - Vd(j)*delta_ex(i))/bsr_x(i);
mdy(j)  = min(Vd(j)/8 + abs(-Vd(j)*e(i) - Vd(j)*delta_ey(i))/(2*bsr_y(i)),Vd(j)/4);

MdRx(i) = rho(i)/100*d(i)^2*fy(i)*(1 - rho(i)/100*fy(i)/(2*fc(i)));
MdRy(i) = rho(i)/100*d(i)^2*fy(i)*(1 - rho(i)/100*fy(i)/(2*fc(i)));

E = 210000; % Assumed elasticity modulus for steel

% Rotation computation
psi_x(j) = 1.5*rs_x(i)*fy(i)./(d(i)*E)*(mdx(j)./MdRx(i)).^1.5;
psi_y(j) = 1.5*rs_y(i)*fy(i)./(d(i)*E)*(mdy(j)./MdRy(i)).^1.5;

% Load computation
Vx(j) = 0.75/(1 + 15*psi_x(j)*d(i)/(16 + dg(i)))*sqrt(fc(i))*bo(i)*d(i);
Vy(j) = 0.75/(1 + 15*psi_y(j)*d(i)/(16 + dg(i)))*sqrt(fc(i))*bo(i)*d(i);

VMC(j) = min(Vx(j),Vy(j));

if VMC(j) == Vx(j)
    direction_psi(j) = 1;
else
    direction_psi(j) = 2;
end

residual(j) = abs(VMC(j) - Vd(j));
end

% Intersection finding between L-psi and Criterion-psi curves
[M(i), I] = min(residual);
m_query(i) = mdx(I);
V_MC_II_kN(i) = Vd(I)/1000;
psi_max_direction(i) = direction_psi(I);
psi(i) = max([psi_x(I) psi_y(I)]);

if psi(i) == psi_x(I)
    Mflex_fib(i) = min(MdRx(i)*((2*c1(i)+c2(i))),0.25*d(i)^2*fc(i)*(2*c1(i)+c2(i)));
else
    Mflex_fib(i) = min(MdRy(i)*c2(i)/2);
end
ratio_V_MC_II(i) = V_exp(i)'/V_MC_II_kN(i); 
else
    V_MC_II_kN(i) = 0;
end
end

Average_V_MC_II = mean(ratio_V_MC_II(ratio_V_MC_II>0))
CoV_V_MC_II = std(V_exp(ratio_V_MC_II>0)./V_MC_II_kN(ratio_V_MC_II>0)')/Average_V_MC_II
Min_V_MC_II = min(V_exp(ratio_V_MC_II>0)./V_MC_II_kN(ratio_V_MC_II>0)')

%% Next generation of Eurocode 2 (FprEN1992-1-1:2022)
for i = tests_i:tests_j 
    
if dg(i) ~= 0 % Apply only with slabs with size aggregate information
    
% variables related with the eccentricity
Ac(i) = c1(i).*c2(i) + 2*c1(i).*d(i)/2 + c2(i).*d(i)/2 + d(i).^2/8*pi;
bb(i) = sqrt((c1(i) + d(i))*(c2(i) + d(i)));
delta_ex(i) = 1/4*(2*c1(i)^2 + 6*c1(i)*d(i) + 3*d(i)^2)/(3*c1(i) + 2*d(i));
delta_ey(i) = 0;

eb(i) = 0.5*abs(-e(i) - delta_ex(i)); 
beta_e(i) = max(1 + 1*eb(i)/bb(i),1.05);  
    
    b0(i) = (2*c1(i) + c2(i) + d(i)*pi/2);
    up = 4;
    kpb(i) = max(min(sqrt(5*up*d(i)/b0(i)),2.5),1);
    ddg(i) = min(16 + dg(i),40);
    tau_Rd(i) = min(0.6*kpb(i)*(rho(i)*fc(i)*ddg(i)/d(i))^(1/3),0.6*sqrt(fc(i)));
    pr_EC2_kN(i) =   tau_Rd(i) * d(i)* b0(i)/beta_e(i)/1000;
    ratio_prEC2(i) = V_exp(i)'/pr_EC2_kN(i);   

else
    pr_EC2_kN(i) = 0;  
end
    
end

Average_pr_EC2 = mean(ratio_prEC2(ratio_prEC2>0))
CoV_pr_EC2 = std(V_exp(ratio_prEC2>0)./pr_EC2_kN(ratio_prEC2>0)')/Average_pr_EC2
Min_pr_EC2 = min(V_exp(ratio_prEC2>0)./pr_EC2_kN(ratio_prEC2>0)')
