% -----------------------------------------------------------------------
% Nom du fichier : bilans_mati�re_energie.m
% Description : Calculs des bilans de mati�re et d'�nergie du proc�d� de 
%               biom�thanisation servant de base de calcul pour le 
%               dimensionnement d'�quipements.
% Auteurs : GREEN MILE (Audrey Collard-Daigneault, ...)
% Date : mars 2021
% -----------------------------------------------------------------------
clc
clear all

% --- Informations connues --- %
n10 = 100;
y10_H2S = 0.25;
y10_H2O = 0.25;
y10_CH4 = 0.25;
y10_CO2 = 0.25;


% Masses molaires
% M_O2 = 0.032; % kg/mol
M_CH4 = 0.01604; % kg/mol
% M_CO2 = 0.04401; % kg/mol
% M_N2 = 0.02801; % kg/mol

% Bilan de mati�re sur R-03 (CH4, CO2, H2O, H2S)
n11 = y10_H2O*n10; % Toute l'eau est retir�e 
n12 = n10 - n11;
y12_CH4 = y10_CH4*n10/n12;
y12_CO2 = y10_CO2*n10/n12;
y12_H2S = y10_H2S*n10/n12;

% Bilan de mati�re sur AD-01 (CH4, CO2, H2S)
n13 = y12_H2S*n12; % Tout le H2S est retir�
n14 = n12 - n13;
y14_CH4 = y12_CH4*n12/n14;
y14_CO2 = y12_CO2*n12/n14;

% Bilan de mati�re sur COMP-01
n15 = n14;
y15_CH4 = y14_CH4;
y15_CO2 = y14_CO2;

% Bilan de mati�re sur HX-03
n18 = n15;
y18_CH4 = y15_CH4;
y18_CO2 = y15_CO2;

% Bilan d'�nergie sur HX-03
n_HX03 = [n15*y15_CH4 n15*y15_CO2]; % Courants �changeur in = out
T15 = 343.15; % K T op�ration AD-01
T18 = 298.15; % K T op�ration AB-01
T16 = 280.15; % K Moyenne temp�rature in de refroidissement (Hall, 2018, p.391)
T17 = 293.15; % K < Temp�rature max vers tour de refroissement (Hall, 2018, p.394) et T op�ration AB-01
T_HX03 = [T15 T18 T16 T17]; % Temp�rature in et out des courants
[n16, Qech_HX03] = bilan_energie_HX03(n_HX03, T_HX03);
n17 = n16;

% Composition par sp�cification de 97% molaire de CH4
y19_CH4 = 0.97;
y19_CO2 = 1 - y19_CH4;
n19 = y18_CH4*n18/y19_CH4;

% Bilan sur le point de s�paration 
sp_21 = 0.10; % 10% massique du m�thane produit est brul�
n21 = sp_21*n19;
n20 = (1 - sp_21)*n19;
y20_CH4 = y19_CH4;
y20_CO2 = y19_CO2;
y21_CH4 = y19_CH4;
y21_CO2 = y19_CO2;

% Par relation d'exc�s de O2 � la fournaise
e = 0.1; % 10% d'exc�s d'air (Perry & Green, 2008, p.24-22)
n22_O2 = 2*y21_CH4*n21*(e+1); % Par relation d'exc�s de O2
n22 = n22_O2/0.21; % Ratio N2/O2 dans l'air sec
n22_N2 = 0.79*n22;

% Bilan de la r�action dans la fournaise
eps = y21_CH4*n21; % Combustion compl�te (on pose)
n25_CH4 = y21_CH4*n21 - eps;
n25_CO2 = (1 - y21_CH4) + 2*eps;
n25_H2O = 2*eps;
n25_N2 = 0.79*n22;
n25_O2 = 0.21*n22 - 2*eps;

% Bilan d'�nergie sur F-01
n_F01 = [y21_CH4*n21 y21_CO2*n21 n25_CO2 n25_O2 n25_H2O n22_N2 n22_O2];
T21 = T18; % K T op�ration AB-01
T25 = 573.15; % K T gaz sortant (A DETERMINER)
T23 = 308.15; % K Eau � 35C minimum (T op�ration r�acteur)
T24 = 373.15; % K Vapeur satur�e � 1 atm
T_F01 = [T21 T25 T23 T24];
[n23, Qech_F01] = bilan_energie_F01(n_F01,T_F01);
