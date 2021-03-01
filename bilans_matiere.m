clc
clear all

% Informations connues
n1 = 100;
y1_H2S = 0.25;
y1_H2O = 0.25;
y1_CH4 = 0.25;
y1_CO2 = 0.25;

% Variable à poser
e = 0.05;

% Masses molaires
% M_O2 = 0.032; % kg/mol
M_CH4 = 0.01604; % kg/mol
% M_CO2 = 0.04401; % kg/mol
% M_N2 = 0.02801; % kg/mol

% Bilan de matière sur R-03 (CH4, CO2, H2O, H2S)
n3 = y1_H2O*n1; % Toute l'eau est retirée 
n2 = n1 - n3;
y2_CH4 = y1_CH4*n1/n2;
y2_CO2 = y1_CO2*n1/n2;
y2_H2S = y1_H2S*n1/n2;

% Bilan de matière sur AD-01 (CH4, CO2, H2S)
n5 = y2_H2S*n2; % Tout le H2S est retiré
n4 = n2 - n5;
y4_CH4 = y2_CH4*n2/n4;
y4_CO2 = y2_CO2*n2/n4;

% Bilan de matière sur COMP-01
n6 = n4;
y6_CH4 = y4_CH4;
y6_CO2 = y4_CO2;

% Bilan de matière sur HX-03
n7 = n6;
y7_CH4 = y6_CH4;
y7_CO2 = y6_CO2;

% Composition par spécification de 97% molaire de CH4
y8_CH4 = 0.97;
y8_CO2 = 1 - y8_CH4;
n8 = y7_CH4*n7/y8_CH4;

% Bilan sur le point de séparation 
sp_10 = 0.10; % 10% massique du méthane produit est brulé
n10 = sp_10*n8;
n9 = (1 - sp_10)*n8;
y9_CH4 = y8_CH4;
y9_CO2 = y8_CO2;
y10_CH4 = y8_CH4;
y10_CO2 = y8_CO2;

% Par relation d'excès de O2 à la fournaise
n11_O2 = y10_CH4*n10/2*(e+1); % Par relation d'excès de O2
n11 = n11_O2/0.21; % Ratio N2/O2 dans l'air sec
n11_N2 = 0.79*n11;

% Bilan de la réaction dans la fournaise
eps = y10_CH4*n10; % Combustion complète (on pose)
n14_CH4 = y10_CH4*n10 - eps;
n14_CO2 = (1 - y10_CH4) + 2*eps;
n14_H2O = 2*eps;
n14_N2 = 0.79*n11;
n14_O2 = 0.21*n11-2*eps;









