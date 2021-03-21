% -----------------------------------------------------------------------
% Nom du fichier : dimensionnement_HX03.m
% Description : Calculs concernant le dimensionnement de l'échangeur de
%               chaleur HX-03.
% Auteurs : GREEN MILE (Audrey Collard-Daigneault)
% Date : mars 2021
% -----------------------------------------------------------------------
clc
clear all

% --- Chaleur à transférer --- %
Q = 3.608429791771280e+04; % W Selon les bilans d'énergie

% --- Informations sur les débits --- %
n15 = 17.644070842975005; % mol/s
y15_CH4 = 0.6243;
y15_CO2 = 1 - y15_CH4;
n16 = 11.023229704450031; % mol/s

% --- Variable posée --- %
U0 = 150; % W/m2K U approximatif moyen pour transfert gaz-eau (Hall, 2018, p.221)

% --- Températures déterminées et pression --- %
T15 = 343.15; % K T opération AD-01 (Okoro & Sun, 2019)
T18 = 288.15; % K T opération AB-01, choix de l'équipe
T16 = 280.15; % K Moyenne température in de refroidissement (Hall, 2018, p.391)
T17 = 322.15; % K Température max vers tour de refroissement (Hall, 2018, p.394)
P = 800000; % Pa Pression d'opération de AB-01
R = 8.314; % m3*Pa/molK
Tm_gaz = (T15+T18)/2; % Température moyenne du gaz
Tm_H2O = (T16+T17)/2; % Température moyenne de l'eau de refroidissement

% --- Dimensions des tubes --- %
Dto = 16/1000; % m Diamètre suggéré (Frauss, 2015, p.23)
Dti = Dto - 2*1.7/1000; % m Diamètre interne (choix arbitraire d'épaisseur selon littérature) (Sinnott, 2020, p.786)

% --- Dimensions des ailettes (fins) --- %
hf = 14/1000; % m Hauteur d'une ailette (Sinnott, 2020, p.890; Thermofin, 2020)
Df = (Dto+2*hf); % m Diamètre ailette
Nf = 315; % Nombre d'ailette par metre de tube (Thermofin, 2020)
ef = 0.5/1000; % m Epaisseur d'une ailette (Frass, 2015, p.51)

% --- Calcul pertinent des aires --- %
At1m = pi*Dto; % m2 Aire d'un tube par 1m;
Ate1m = At1m - Nf*ef*pi*Dto; % m2 Aire tube sans ailette
Af1m = Nf*(2*pi*(Df^2 - Dto^2)/4 + ef*pi*Df); % m2 Aire d'une ailette
Atot1m = Ate1m + Af1m;
Nrangees = 8; % Nombre de rangées (Frass, 2015, p.53)
espace = 5/1000; % m Espace arbitraire de 5 mm entre les ailettes des tubes
Ntubes = 60; % 7-8 tubes par rangées
Ntprmax = 8; % Nombre max de tube pour une rangée
Hech = Ntprmax*(Df+espace); % Hauteur de l'échangeur avec ailettes collées
Lech = Nrangees*sqrt(3)/2*(Df+espace); % Longueur de l'échangeur

% --- Autre infos --- %
Rfo = 0.0003; % m2C/W Facteur de Fouling max pour de l'eau de refroidissement (Sinnott, 2020, p. 780)
Rfi = 0.0002; % m2C/W Facteur de Fouling pour les hydrocarbures légers (Sinnott, 2020, p. 780)
k_al = 177; % W/mK Conductivité thermique alliage d'aluminium

% --- Propriétés du CH4 et CO2 --- %
k_CO2 = -6.07829291578281E-03 + 7.53509863170161E-05*Tm_gaz + ...
    9.49276579202504E-09*Tm_gaz^2 + -1.12751601050646E-11*Tm_gaz^3; % W/mK  (Yaws, 2014, table 84)
k_CH4 = 5.37671485384522E-03 + 5.15550830467679E-05*Tm_gaz + ...
    1.66548977272723E-07*Tm_gaz^2 + -5.71678321678305E-11*Tm_gaz^3; % W/mK (Yaws, 2014, table 85)
Cp_CH4 = 46.241914499032198 + -0.15641465372*Tm_gaz + 0.00063340927021*Tm_gaz^2 ...
    + -9.8938057655e-007*Tm_gaz^3 + 8.9552744542e-010*Tm_gaz^4; % J/molK (Yaws, 2014, table 165)
Cp_CO2 = 22.870906792800799 + 0.053885041357*Tm_gaz + -3.0277722822e-006*Tm_gaz^2 ...
    + -5.9414717025e-008*Tm_gaz^3 + 5.758653574e-011*Tm_gaz^4; % J/molK (Yaws, 2014, table 165)

% --- Propriétés du biogaz et calculs intermédiaires --- %
M_CH4 = 16.04/1000; % kg/mol
M_CO2 = 44.01/1000; % kg/mol
M_gaz = y15_CH4*M_CH4 + y15_CO2*M_CO2; % kg/mol Masse molaire du biogaz
rho_gaz = P*M_gaz/(R*Tm_gaz); % kg/m3 Masse volumique (hypothèse de fluide incompressible)
mu_gaz = 1.03e-5; % Pa*s Viscosité dynamique (Gallagher, 2006, p.51)
k_gaz = 1/(y15_CH4/k_CH4 + y15_CO2/k_CO2); % W/mK Conductivité thermique du biogaz (Gilmore & Comings, 1966)
Cp_gaz = y15_CH4*Cp_CH4/M_CH4 + y15_CO2*Cp_CO2/M_CO2; % J/kgK 
Pr_gaz = mu_gaz*Cp_gaz/k_gaz; % Nombre de Prandtl
Q_gaz = n15*M_gaz/rho_gaz; % m3/s Débit volumétrique du biogaz

% --- Propriétés de l'eau de refroidissement et calculs intermédiaires --- %
M_H2O = 0.018015; % kg/mol Masse molaire
rho_H2O = 997; % kg/m3 à Tf = 300K Masse volumique (Bergman, 2018, Table A6) 
Cp_H2O = 4179; % J/KgK à Tf = 300K Cp (Bergman, 2018, Table A6) 
k_H2O = 613e-3; % W/mK à Tf = 300K Conductivité thermique (Bergman, 2018, Table A6) 
mu_H2O = 855e-6; % Pa*s à Tf = 300K Viscosité dynamique (Bergman, 2018, Table A6) 
Pr_H2O = 5.83; % à Tf = 300K  Nombre de Prandtl (Bergman, 2018, Table A6) 
Q_H2O = n16*M_H2O/rho_H2O; % m3/s Débit volumétrique 

% --- Calcul du facteur de correction courant croisé --- %
Tci = T16;
Tco = T17;
Thi = T15;
Tho = T18;
R_ = (Thi - Tho)/(Tco - Tci);
P_ = (Tco - Tci)/(Thi - Tci);
F = 0.595; % Facteur de correction (Bergman, 2018, fig 11S.3)

% --- Calcul de la différence moyenne de température logarithmique --- %
dT1 = Thi-Tco; 
dT2 = Tho-Tci;
dTlm = (dT2-dT1)/log(dT2/dT1);

% --- Début des intérations avec le U --- %
U = 0;
Ucal = U0;
tol = 5;

while abs(Ucal - U) > tol
    U = Ucal;
    % --- Calcul de l'aire et de rarrangement des tubes --- %
    A = Q/(F*Ucal*dTlm) % m2 Surface d'échange entre les fluides
    %Ntubes = floor(A/(Atot1m*Wech)); % Nombre de tube
    %Ntpr = floor(Ntubes/Nrangees); % Nombre de tube par rangée
    Wech = A/(Ntubes*Atot1m)
    
    % --- Calcul du coefficient de transfert thermique côté gaz --- %
    Alibre = Wech*(Hech - Ntprmax*(Dto+Nf*Df*ef)); % Aire de passage du gaz
    v_gaz = Q_gaz/Alibre; % m/s 
    L_gaz = Dto*Atot1m/At1m; % Longueur caractéristique 
    Re_gaz = rho_gaz*v_gaz*L_gaz/mu_gaz; % Nombre de Reynolds
    Nu_gaz = 0.45*Re_gaz^0.625*Pr_gaz^(1/3)*(Atot1m/At1m)^-0.375; % Nombre de Nusselt (Frauss, 2015, p.12)
    ho = k_gaz*Nu_gaz/L_gaz; % W/m2K h côté gaz
    
    % Nombre de Reynolds
    Q_H2O = n16*M_H2O/rho_H2O; % m3/s 
    v_H2O = Q_H2O/(Ntubes*pi*Dti^2/4); % m/s Vitesse de l'eau dans les tubes
    Re_H2O = rho_H2O*v_H2O*Dti/mu_H2O; 
    Nu_H2O = 0.023*Re_H2O^(4/5)*Pr_H2O^0.3; % (Bergman, 2018, p.519)
    hi = k_H2O*Nu_H2O/Dti;
    
    Ai = pi*Dti*Wech*Ntubes;
    Ao = Atot1m*Wech*Ntubes;
    nf = 0.9;
    
    Ucal = 1/(A*((1/(hi*Ai) + Rfi/Ai + Rfo/Ao + 1/(ho*Ao))/nf + log(Dto/Dti)/(2*pi*k_al*Wech*Ntubes)))

    
%     U0 = 1/(1/0.95*(1/ho+1/Rfo)*Atot1m*Wech/Af1m + Dto*log(Dto/Dti)/(2*k_al) + Dto/(Dti*Rfi) + Dto/(Dti*hi));
%     Ucal = Atot1m*Wech*U0/A
    
    

    
    
  
end
