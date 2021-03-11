clc 
clear all

% Masse molaires
M_h2O = 0.01802 %kg/mol

% Bilan de matière sur AB-01
y19_CH4 = 0.97;
y19_CO2 = 0.03;
n19 = n18*y18_CH4/y19_CH4; 
n26_CO2 = (n18_C02-n19_C02)/(1-0.01);
n28_CO2 = 0.01*n26_CO2;
n36_H2O = n28_H2O; 
n28_H2O = n26_H2O;
n28 = n28_CO2+n28_H2O;
n26 = n26_CO2+n26_H2O; 
y26_H2O = n26_H2O/n26; 
y26_CO2 = 1-y26_H2O; 
y28_H2O = n28_H2O/n28; 
y28_CO2 = 1-y28_H2O;
