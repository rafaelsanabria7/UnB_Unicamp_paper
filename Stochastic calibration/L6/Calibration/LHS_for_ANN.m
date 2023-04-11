clc; clear all; close all
%LHS for punching shear model calibration 09-09-2021
%Slab L6 (Albuquerque,2012)
%====================================================================================================%
%                                              LHS                                                   %
%                                                                                                    %
%                                   latin Hypercube Sampling                                         %
%====================================================================================================%
%Units: kN and meters.
N_SIM=str2double(inputdlg('Insert the number of Samples'));

%PARAMETERS AND RANDOM VARIABLES:
%Assuming all variables with normal distribution
%Compression resistence of concrete, f'c MPa
mx_fc=52.1;
vx_fc=0.0632;
sigx_fc=mx_fc*vx_fc;

%Modulus of Eslasticity, Ec MPa
mx_Ec= 33600; %0.9*21500*((mx_fc/10)^(1/3));
vx_Ec=0.15;
sigx_Ec=mx_Ec*vx_Ec;

%Tensile Strength, fct Mpa
mx_fct= 3.85; %2.12*log(1 + 0.1*mx_fc);%0.3*((mx_fc - sigx_fc)^(2/3));
vx_fct= 0.2015;
sigx_fct=mx_fct*vx_fct;

%Fracture of Energy , Gf 
mx_Gf= 133; %0.9*73*((mx_fc)^0.18); %/sqrt(2);
vx_Gf=0.15;
sigx_Gf=mx_Gf*vx_Gf;

% Spring Support, K
mx_K= 120000;
vx_K=0.15;
sigx_K=mx_K*vx_K;

%SAMPLES OBTAINED BY LHS:
lhs_cdf = lhsdesign(N_SIM,5,'criterion','correlation');
fc = norminv(lhs_cdf(:,1),mx_fc,sigx_fc);
fct = norminv(lhs_cdf(:,2),mx_fct,sigx_fct);
Ec = norminv(lhs_cdf(:,3),mx_Ec,sigx_Ec);
Gf = norminv(lhs_cdf(:,4),mx_Gf,sigx_Gf);
K = norminv(lhs_cdf(:,5),mx_K,sigx_K);

%Variables to be correlated in standard normal space
fc_std = (fc-mx_fc)./sigx_fc;
Ec_std = (Ec-mx_Ec)./sigx_Ec;
fct_std = (fct-mx_fct)./sigx_fct;
Gf_std = (Gf-mx_Gf)./sigx_Gf;

%Correlationing the variables f'c, Ec, Fct, and Gf
%1. Create the random variables in the Standard Normal Space
%2. Covariance Matrix
rho_fc_fct = 0.8;
rho_fct_Ec = 0.7;
rho_fc_Ec = 0.9;
rho_fc_Gf = 0.6;
rho_fct_Gf = 0.9;
rho_Ec_Gf = 0.5;

Cov = corr2cov([std(fc_std);std(fct_std);std(Ec_std);std(Gf_std)],[1 rho_fc_fct rho_fc_Ec rho_fc_Gf;rho_fc_fct 1 rho_fct_Ec rho_fct_Gf;rho_fc_Ec rho_fct_Ec 1 rho_Ec_Gf;rho_fc_Gf rho_fct_Gf rho_Ec_Gf 1]);

%3. Calculate Cholesky
l = chol(Cov);
correlated_variables = [fc_std,fct_std,Ec_std,Gf_std]*l;

%Variáveis aleatorias correlacionadas
fc = correlated_variables(:,1).*sigx_fc + mx_fc;
fct = correlated_variables(:,2).*sigx_fct + mx_fct;
Ec = correlated_variables(:,3).*sigx_Ec + mx_Ec;
Gf = correlated_variables(:,4).*sigx_Gf + mx_Gf;

%Gerar .txt
var = [fc fct Ec Gf K];
fid=fopen('input.txt','w+t');%Abre um novo arquivo .txt para a escritura dos resultados
fprintf(fid,'%2.3f %1.3f %1.3f %3.3f %3.3f\n',var');
fclose(fid);



