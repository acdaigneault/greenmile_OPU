% -----------------------------------------------------------------------
% Nom du fichier : bilans_mati�re_energie.m
% Description : Calculs des bilans de mati�re et d'�nergie du proc�d� de 
%               biom�thanisation servant de base de calcul pour le 
%               dimensionnement d'�quipements.
% Auteurs : GREEN MILE (Audrey Collard-Daigneault, Jesse Giroux)
% Date : mars 2021
% -----------------------------------------------------------------------
clc
clear all

%% R�solution de l'�quation de Buswell et Boyle

% Coefficients de l'�quation de Buswell et Boyle
a = 1054;
b = 1916;
c = 397;
d = 24;
e = 1;

% Coefficients stoechiom�triques de la r�action globale de RX-01
v_H2O = a-b/4-c/2+3*d/4+e/2;
v_CH4 = a/2+b/8-c/4-3*d/8-e/4;
v_CO2 = a/2-b/8+c/4+3*d/8+e/4;
v_NH3 = d;
v_H2S = e;


%% Masses molaires
M_C = 12.0107; %g/mol
M_H = 1.0078; %g/mol
M_O = 15.9994; %g/mol
M_N = 14.0067; %g/mol
M_S = 32.065; %g/mol
M_MO = 1054*M_C+1916*M_H+397*M_O+24*M_N+M_S; %g/mol (M suppos�e de la MO pour les calculs de Buswell et Boyle)
M_H2O = 2*M_H+M_O; %g/mol
M_CH4 = M_C+4*M_H; %g/mol
M_CO2 = M_C+2*M_O; %g/mol
M_NH3 = M_N+3*M_H; %g/mol
M_H2S = 2*M_H+M_S; %g/mol
M_O2 = 2*M_O; %g/mol
M_N2 = 2*M_N; % g/mol


%% Bilans de masse sur les conduites 1 � 9

% Conduite 1
m1_MO = 75000*10^6/365/24/3600; %g/s (Base de calcul)
n1_MO = m1_MO/M_MO; %mol/s (D�bit molaire suppos� pour les calculs de Buswell et Boyle)
m1_H2O = 75/25*m1_MO; %g/s (Taux de siccit� de 25% de la mati�re qui entre dans l'usine)
n1_H2O = m1_H2O/M_H2O; %mol/s

% Conduite 3
m3_MO = m1_MO; %g/s
n3_MO = n1_MO; %mol/s
m3_H2O = 85/15*m3_MO; %g/s (Taux de siccit� de 15% de la mati�re qui entre dans RX-01)
n3_H2O = m3_H2O/M_H2O; %mol/s

% Conduite 5
% Hypoth�se: La totalit� des produits de la r�action se trouvent dans le biogaz sorant (Conduite 5) de RX-01
nRx_MO = 0.6*n3_MO; %mol/s (40% de la MO est non dig�rable)
n5_CH4 = v_CH4*nRx_MO; %mol/s
n5_CO2 = v_CO2*nRx_MO; %mol/s
n5_NH3 = v_NH3*nRx_MO; %mol/s
n5_NH3 = 0; %mol/s (L'Ammoniac contenu dans le biogaz est n�gli�e)
n5_H2S = v_H2S*nRx_MO; %mol/s
Pv_H20 = 0.05562; %atm (Pv de l'eau � 35�C)(Yaws, 2009)
y5_H2O = 0.05562; % (Hypoth�se: Le biogaz sortant de RX-01 est un gaz parfait � pression atmosph�rique)
n5_H2O = y5_H2O*(n5_CH4+n5_CO2+n5_NH3+n5_H2S)/(1-y5_H2O); %mol/s

% Conduite 4
n4_MO = 0.4*n3_MO; %mol/s (40% de la MO est non dig�rable)
n4_H2O = n3_H2O-v_H2O*nRx_MO-n5_H2O; %mol/s

% Conduite 6
n6_MO = n4_MO; %mol/s (Hypoth�se: La totalit� de la MO sort par le bas de DC-01)
m6_MO = n6_MO*M_MO; %g/s
m6_H2O = 90/10*m6_MO; %g/s (Taux de siccit� de 10% de la mati�re qui sort de DC-01)
n6_H2O = m6_H2O/M_H2O; %mol/s

% Conduite 7
n7_H2O = n4_H2O-n6_H2O; %mol/s

% Conduite 8
n8_H2O=0.9*n7_H2O; %mol/s (10% de l'eau de proc�d� peut �tre recircul�e au digesteur)

% Conduite 9
n9_H2O=0.1*n7_H2O; %mol/s (10% de l'eau de proc�d� peut �tre recircul�e au digesteur)

% Conduite 2
n2_H2O = n3_H2O-n1_H2O-n9_H2O; %mol/s


%% Bilan sur les conduites 10 � ?

% Conduite 10
n10 = n5_CH4+n5_CO2+n5_H2S+n5_H2O;  %mol/s
y10_CH4 = n5_CH4/n10;
y10_CO2 = n5_CO2/n10;
y10_H2S = n5_H2S/n10;
y10_H2O = n5_H2O/n10;


% Bilan de mati�re sur R-03 (CH4, CO2, H2O, H2S)
n11= y10_H2O*n10; %mol/s (Toute l'eau est retir�e du biogaz par R-03) 
n12 = n10 - n11; %mol/s
y12_CH4 = y10_CH4*n10/n12;
y12_CO2 = y10_CO2*n10/n12;
y12_H2S = y10_H2S*n10/n12;

% Bilan de mati�re sur AD-01 (CH4, CO2, H2S)
n13 = y12_H2S*n12; %mol/s (Tout le H2S est retir� du biogaz par AD-01A/B)
n14 = n12 - n13; %mol/s
y14_CH4 = y12_CH4*n12/n14;
y14_CO2 = y12_CO2*n12/n14;

% Bilan de mati�re sur COMP-01
n15 = n14; %mol/s
y15_CH4 = y14_CH4;
y15_CO2 = y14_CO2;

% Bilan de mati�re sur HX-03
n18 = n15; %mol/s
y18_CH4 = y15_CH4;
y18_CO2 = y15_CO2;

% Bilan d'�nergie sur HX-03
n_HX03 = [n15*y15_CH4 n15*y15_CO2]; %mol/s (Courants �changeur in = out)
T15 = 343.15; %K (T op�ration AD-01)
T18 = 298.15; %K (T op�ration AB-01)
T16 = 280.15; %K (Moyenne temp�rature in de refroidissement)(Hall, 2018, p.391)
T17 = 293.15; %K (Temp�rature max vers tour de refroissement et T op�ration AB-01)(Hall, 2018, p.394)
T_HX03 = [T15 T18 T16 T17]; %K (Temp�rature in et out des courants)
[n16, Qech_HX03] = bilan_energie_HX03(n_HX03, T_HX03);
n17 = n16; %mol/s

% Composition par sp�cification de 97% molaire de CH4
y19_CH4 = 0.97;
y19_CO2 = 1 - y19_CH4;
n19 = y18_CH4*n18/y19_CH4; %mol/s

% Bilan sur le point de s�paration 
sp_21 = 0.10; % (10% massique du m�thane produit est brul�)
n21 = sp_21*n19; %mol/s
n20 = (1 - sp_21)*n19; %mol/s
y20_CH4 = y19_CH4;
y20_CO2 = y19_CO2;
y21_CH4 = y19_CH4;
y21_CO2 = y19_CO2;

% Par relation d'exc�s de O2 � la fournaise
e = 0.1; % (10% d'exc�s d'air)(Perry & Green, 2008, p.24-22)
n22_O2 = 2*y21_CH4*n21*(e+1); %mol/s (Par relation d'exc�s de O2)
n22 = n22_O2/0.21; %mol/s (Ratio N2/O2 dans l'air sec)
n22_N2 = 0.79*n22; %mol/s

% Bilan de la r�action dans la fournaise
eps = y21_CH4*n21; %mol/s (Hypoth�se: Combustion compl�te)
n25_CH4 = y21_CH4*n21 - eps; %mol/s
n25_CO2 = (1 - y21_CH4) + 2*eps; %mol/s
n25_H2O = 2*eps; %mol/s
n25_N2 = 0.79*n22; %mol/s
n25_O2 = 0.21*n22 - 2*eps; %mol/s

% Bilan d'�nergie sur F-01
n_F01 = [y21_CH4*n21 y21_CO2*n21 n25_CO2 n25_O2 n25_H2O n22_N2 n22_O2];
T21 = T18; %K (T op�ration AB-01)
T25 = 573.15; %K (T gaz sortant (A DETERMINER))
T23 = 308.15; %K (Eau � 35C minimum (T op�ration r�acteur))
T24 = 373.15; %K (Vapeur satur�e � 1 atm)
T_F01 = [T21 T25 T23 T24];
[n23, Qech_F01] = bilan_energie_F01(n_F01,T_F01);


%% Bilan sur les conduites ? � ? (Marieme et Nawshin)


