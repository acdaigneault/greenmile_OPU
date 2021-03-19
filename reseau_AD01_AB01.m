% -----------------------------------------------------------------------
% Nom du fichier : reseau_AD01_AB01.m
% Description : Calculs concernant le réseau de conduites entre le AD-01 
%               et AB-01 avec les équipements hydrauliques reliés.
% Auteurs : GREEN MILE (Audrey Collard-Daigneault)
% Date : mars 2021
% -----------------------------------------------------------------------
clc
clear all

% Informations sur le fluide
n = 70.576283371900020; % mol/s Débit molaire de biogaz
y_CH4 = 0.6243; % Fraction molaire de CH4
y_CO2 = 1 - y_CH4; % Fraction molaire de CO2
M_CH4 = 0.01604; % kg/mol
M_CO2 = 0.04401; % kg/mol
M_gaz = y_CH4*M_CH4 + y_CO2*M_CO2; % kg/mol
nu_gaz = 1.3; % Facteur de compressibilité du biogaz (Gallagher, 2006, p.52)
mu_gaz = 1.03e-5; % Pa*s Viscosité dynamique (Gallagher, 2006, p.51)
R = 8.314; % Pa*m3/molK

%% Conduite 14 du réseau (entre AD-01 et COMP-01)
% Conditions (Hypothèses : Isotherme et isobare)
T14 = 70 + 273.15; % K Température d'opération de AD-01
P14 = 450000; % Pa Pression d'opération de AD-01

% Vérification fluide incompressible
a14 = sqrt(nu_gaz*T14*R/M_gaz); % m/s vitesse du son 
v0 = 30; % m/s Vitesse posée pour un gaz entre 15 et 40 (Hall, 2018, p.51)

% Calcul du diamètre
rho14 = P14*M_gaz/(R*T14); % Masse volumique
V14 = n*M_gaz/rho14; % Débit volumique 
Dp14 = sqrt(4*V14/(v0*pi)); % m % Diamètre nécessaire pour la vitesse posée
D14 = 0.1282; % m Diamètre ajusté en fonction des standards (Nominal 5po schedule 40)
v14 = 4*V14/(pi*D14^2); % m Vitesse calculée

%% Conduite 15 du réseau (entre COMP-01 et HX-03)
% Conditions
T15 = T14; % K Température d'opération de AD-01 et compression isotherme
P15 = 800000; % Pa Pression d'opération de AB-01

% Calcul du diamètre
rho15 = P15*M_gaz/(R*T15); % Masse volumique
V15 = n*M_gaz/rho15; % Débit volumique 
Dp15 = sqrt(4*V15/(v0*pi)); % m % Diamètre nécessaire pour la vitesse posée
D15 = 0.1023; % m Diamètre ajusté en fonction des standards (Nominal 4po schedule 40)
v15 = 4*V15/(pi*D15^2); % m Vitesse calculée

%% Conduite 18 du réseau (entre HX-03 et AB-01)
% Conditions
T18 = 15 + 273.15; % K Température d'opération de AB-01 
P18 = 800000; % Pa Pression d'opération de AB-01 (dP négligée dans l'échangeur)************

% Calcul du diamètre
rho18 = P18*M_gaz/(R*T18); % Masse volumique
V18 = n*M_gaz/rho18; % Débit volumique 
Dp18 = sqrt(4*V18/(v0*pi)); % m % Diamètre nécessaire pour la vitesse posée
D18 = 0.0901; % m Diamètre ajusté en fonction des standards (Nominal 3 1/2po schedule 40)
v18 = 4*V18/(pi*D18^2); % m Vitesse calculée

