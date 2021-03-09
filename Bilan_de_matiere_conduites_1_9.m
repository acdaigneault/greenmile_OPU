clear all;clc;

%% Résolution de l'équation de Buswell et Boyle

% Coefficients de l'équation de Buswell et Boyle
a = 1054;
b = 1916;
c = 397;
d = 24;
e = 1;

% Coefficients stoechiométriques de la réaction globale de RX-01
v_H2O = a-b/4-c/2+3*d/4+e/2;
v_CH4 = a/2+b/8-c/4-3*d/8-e/4;
v_CO2 = a/2-b/8+c/4+3*d/8+e/4;
v_NH3 = d;
v_H2S = e;

%% Bilans de masse sur les conduites 1 à 9

%Masses molaires
M_C = 12.0107; %g/mol
M_H = 1.0078; %g/mol
M_O = 15.9994; %g/mol
M_N = 14.0067; %g/mol
M_S = 32.065; %g/mol
M_MO = 1054*M_C+1916*M_H+397*M_O+24*M_N+M_S; %g/mol (M supposée de la MO pour les calculs de Buswell et Boyle)
M_H2O = 2*M_H+M_O; %g/mol
M_CH4 = M_C+4*M_H; %g/mol
M_CO2 = M_C+2*M_O; %g/mol
M_NH3 = M_N+3*M_H; %g/mol
M_H2S = 2*M_H+M_S; %g/mol

% Conduite 1
m1_MO = 75000*10^6/365/24; %g/h (Base de calcul)
n1_MO = m1_MO/M_MO; %mol/h (Débit molaire supposé pour les calculs de Buswell et Boyle)
m1_H2O = 75/25*m1_MO; %g/h (Taux de siccité de 25% de la matière qui entre dans l'usine)
n1_H2O = m1_H2O/M_H2O; %mol/h

% Conduite 3
m3_MO = m1_MO; %g/h
n3_MO = n1_MO; %mol/h
m3_H2O = 85/15*m3_MO; %g/h (Taux de siccité de 15% de la matière qui entre dans RX-01)
n3_H2O = m3_H2O/M_H2O; %mol/h

% Conduite 5
% Hypothèse: La totalité des produits de la réaction se trouvent dans le biogaz sorant (Conduite 5) de RX-01
nRx_MO = 0.6*n3_MO; %mol/h (40% de la MO est non digérable)
n5_CH4 = v_CH4*nRx_MO; %mol/h
n5_CO2 = v_CO2*nRx_MO; %mol/h
n5_NH3 = v_NH3*nRx_MO; %mol/h
n5_H2S = v_H2S*nRx_MO; %mol/h
Pv_H20 = 0.05562; %atm (Pv de l'eau à 35°C)(Yaws, 2009)
y5_H2O = 0.05562; % (Hypothèse: Le biogaz sortant de RX-01 est un gaz parfait à pression atmosphérique)
n5_H2O = y5_H2O*(n5_CH4+n5_CO2+n5_NH3+n5_H2S)/(1-y5_H2O); %mol/h

% Conduite 4
n4_MO = 0.4*n3_MO; %mol/h (40% de la MO est non digérable)
n4_H2O = n3_H2O-v_H2O*nRx_MO-n5_H2O; %mol/h

% Conduite 6
n6_MO = n4_MO; %mol/h (Hypothèse: La totalité de la MO sort par le bas de DC-01)
m6_MO = n6_MO*M_MO; %g/h
m6_H2O = 90/10*m6_MO; %g/h (Taux de siccité de 10% de la matière qui sort de DC-01)
n6_H2O = m6_H2O/M_H2O; %mol/h

% Conduite 7
n7_H2O = n4_H2O-n6_H2O; %mol/h

% Conduite 8
n8_H2O=0.9*n7_H2O; %mol/h (10% de l'eau de procédé peut être recirculée au digesteur)

% Conduite 9
n9_H2O=0.1*n7_H2O; %mol/h (10% de l'eau de procédé peut être recirculée au digesteur)

% Conduite 2
n2_H2O = n3_H2O-n1_H2O-n9_H2O; %mol/h
