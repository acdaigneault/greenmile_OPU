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
Q = 39512; % W Selon les bilans d'énergie
eff = 1.1; % Facteur de sécurité si perte vers environnement
Q = eff*Q; % W

% --- Informations sur les débits --- %
n15 = 19.3203; % mol/s
y15_CH4 = 0.6243;
y15_CO2 = 1 - y15_CH4;
n16 = 28.9135; % mol/s

% --- Variable posée --- %
U0 = 150; % W/m2K U approximatif pour un trasnfert gaz-eau (Hall, 2018, p.221)

% --- Températures déterminées et pression --- %
T15 = 343.15; % K T opération AD-01 (Okoro & Sun, 2019)
T18 = 288.15; % K T opération AB-01, choix de l'équipe
T16 = 280.15; % K Approximation température in de refroidissement (Hall, 2018, p.391)
T17 = 298.15; % K Température vers tour de refroissement ou autre unités
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

% --- Calculs pertinents des aires --- %
At1m = pi*Dto; % m2 Aire d'un tube par 1m;
Ate1m = At1m - Nf*ef*pi*Dto; % m2 Aire tube sans ailette pour 1m
Af1m = Nf*(2*pi*(Df^2 - Dto^2)/4 + ef*pi*Df); % m2 Aire d'une ailette pour 1m
Atot1m = Ate1m + Af1m; % m2 Aire totale extérieure pour 1m
Nrangees = 8; % Nombre de rangées (Frass, 2015, p.53)
espace = 5/1000; % m Espace horizontal arbitraire de 5 mm entre les ailettes des tubes
Ntubes = 36; % 5-4 tubes par rangées
Ntprmax = 5; % Nombre max de tube pour une rangée
Hech = Ntprmax*(Df+espace)+espace; % Hauteur de l'échangeur 
Lech = Nrangees*(sqrt(3)/2*Df+espace)+espace; % Longueur de l'échangeur
pf = 1/Nf; % m Pitch de l'ailette (distance entre deux ailettes)

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
rho_H2O = 999; % kg/m3 à Tf = 300K Masse volumique (Bergman, 2018, Table A6) 
Cp_H2O = 4184; % J/KgK à Tf = 300K Cp (Bergman, 2018, Table A6) 
k_H2O = 598e-3; % W/mK à Tf = 300K Conductivité thermique (Bergman, 2018, Table A6) 
mu_H2O = 1080e-6; % Pa*s à Tf = 300K Viscosité dynamique (Bergman, 2018, Table A6) 
Pr_H2O = 7.56; % à Tf = 300K  Nombre de Prandtl (Bergman, 2018, Table A6) 
Q_H2O = n16*M_H2O/rho_H2O; % m3/s Débit volumétrique 

% --- Méthode eNTU --- %
Tci = T16;
Tco = T17;
Thi = T15;
Tho = T18;

Cc = Q/(Tco-Tci); 
Ch = Q/(Thi-Tho);

if Cc < Ch
    Cmin = Cc;
    Cmax = Ch;
else 
    Cmin = Ch;
    Cmax = Cc;
end

Cr = Cmin/Cmax;
e = Q/(Cmin*(Thi-Tci));

% Calcul NTU selon C mixé
if abs(Cmin - Ch) < 1e-6
    NTU = -log(Cr*log(1-e)+1)/Cr; % Cmin mixé
else 
    NTU = -log(1+log(1-e*Cr)/Cr); % Cmax mixé
end

% --- Début des intérations avec le U --- %
U = 0;
Ucal = U0;
tol = 0.01;

while abs(Ucal - U) > tol
    U = Ucal;
    % --- Calcul de l'aire et de rarrangement des tubes --- %
    Ao = Cmin*NTU/U;
    Wech = Ao/(Ntubes*Atot1m);
    Lt = Wech;
    
    % --- Calcul du coefficient de transfert thermique côté gaz --- %
    Alibre = Wech*(Hech - Ntprmax*(Dto+Nf*2*hf*ef)); % Aire de passage du gaz
    v_gaz = Q_gaz/Alibre; % m/s 
    L_gaz = Dto; % Longueur caractéristique 
    Re_gaz = rho_gaz*v_gaz*L_gaz/mu_gaz; % Nombre de Reynolds
    Nu_gaz = 0.134*Re_gaz^0.681*Pr_gaz^0.33*((pf-ef)/hf)^0.2*(pf/ef)^0.1134; % (Sinnott, 2020, p.891)
    ho = k_gaz*Nu_gaz/L_gaz; % W/m2K h côté gaz
    
    % --- Calcul du coefficient de transfert thermique côté eau --- %
    Q_H2O = n16*M_H2O/rho_H2O; % m3/s 
    v_H2O = Q_H2O/(Ntubes*pi*Dti^2/4); % m/s Vitesse de l'eau dans les tubes
    Re_H2O = rho_H2O*v_H2O*Dti/mu_H2O; 
    Nu_H2O = 0.023*Re_H2O^(4/5)*Pr_H2O^0.3; % (Bergman, 2018, p.519)
    hi = k_H2O*Nu_H2O/Dti;
    
    % --- Calcul de l'aire interne --- %
    Ai = pi*Dti*Lt*Ntubes; % m2 Aire totale dans les tubes
    
    % --- Calcul de l'efficacité des ailettes --- %
    m = sqrt(2*ho/(k_al*ef));
    nf = tanh(m*hf)/(m*hf); % Efficacité d'une ailette
    no = 1 - (Af1m/Atot1m)*(1-nf); % Efficacité d'une surface
    
    % --- Recalcul du coefficient de transfert global côté gaz --- %
    Ucal = 1/((1/ho+Rfo)/no*(At1m/Af1m) + (Dto*log(Dto/Dti))/(2*k_al) + Dto/Dti*(Rfi + 1/hi)); % W/m2K (Sinnott, 2020, p.775 et p.891)
end

% --- Perte de charge dans les tubes --- %
Uo = Ucal;
Ui = Dto*Uo/Dti; % Wm2
T_wall = (Ui*Tm_gaz + Tm_H2O*(hi - Ui))/hi; % K Température estimée de la surface (Sinnott, 2020, p.807)
mu_wall = 721.2e-6; % Pas Viscosité à la surface (Bergman, 2018, p.919)
jf = 8.5e-3;
M = 0.25;
dP_tube = 8*jf*(Lt/Dti)*rho_H2O*v_H2O^2/2*(mu_H2O/mu_wall)^-M;

Dti2 = 0.0779; % Conduite d'entrée d'eau 3po cédule 40 (Crane)
fm = 0.042; % (Crane)
Kin = 0.78;
Kout = 1;
dP = dP_tube + rho_H2O*v_H2O^2/2*Ntubes*(Kin+Kout) + rho_H2O*(4*Q_H2O/(pi*Dti2^2))^2/2*(Kin+Kout); % Pa Perte de charge

fprintf('Hauteur : %.4f m\n', Hech)
fprintf('Largueur : %.4f m\n', Wech)
fprintf('Longueur : %.4f m\n', Lech)
fprintf('U : %.4f W/m2K \n', Ucal)
fprintf('Surface extérieur : %.4f m2\n', Ao)
fprintf('Coefficient de transfert : %.4f W/m2K\n', ho)
fprintf('Perte de charge côté eau/tube : %.4f Pa\n', dP)

