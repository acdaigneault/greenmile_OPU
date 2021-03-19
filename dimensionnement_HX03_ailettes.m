% -----------------------------------------------------------------------
% Nom du fichier : dimensionnement_HX03.m
% Description : Calculs concernant le dimensionnement de l'�changeur de
%               chaleur HX-03.
% Auteurs : GREEN MILE (Audrey Collard-Daigneault)
% Date : mars 2021
% -----------------------------------------------------------------------
clc
clear all

% --- Chaleur � transf�rer --- %
Q = -1.443371916708512e+05; % W Selon les bilans d'�nergie

% --- Informations sur les d�bits --- %
n15 = 70.576283371900020; % mol/s
y15_CH4 = 0.6243;
y15_CO2 = 1 - y15_CH4;
n16 = 44.092918817800125; % mol/s

% --- Variable pos�e --- %
U0 = 150; % W/m2K U approximatif moyen pour transfert gaz-eau (Hall, 2018, p.221)

% --- Temp�ratures d�termin�es et pression --- %
T15 = 343.15; % K T op�ration AD-01 (Okoro & Sun, 2019)
T18 = 288.15; % K T op�ration AB-01, choix de l'�quipe
T16 = 280.15; % K Moyenne temp�rature in de refroidissement (Hall, 2018, p.391)
T17 = 322.15; % K Temp�rature max vers tour de refroissement (Hall, 2018, p.394)
P = 800000; % Pa Pression d'op�ration de AB-01
R = 8.314; % m3*Pa/molK
Tm_gaz = (T15+T18)/2; % Temp�rature moyenne du gaz
Tm_H2O = (T16+T17)/2; % Temp�rature moyenne de l'eau de refroidissement

% --- Informations de conception pour les tubes --- %
N = 1; % Nombre de paire de tubes
Rfo = 0.0003; % m2C/W Facteur de Fouling max pour de l'eau de refroidissement (Sinnott, 2020, p. 780)
Rfi = 0.0002; % m2C/W Facteur de Fouling pour les hydrocarbures l�gers (Sinnott, 2020, p. 780)
Dto = 38/1000; % m Diam�tre sugg�r� (Frauss, 2015, p.23)
Dti = Dto - 2*1.7/1000; % m Diam�tre interne (choix arbitraire d'�paisseur selon litt�rature) (Sinnott, 2020, p.786)
Pf = 3/1000; % m Distance entre les ailettes circulaires et recommand� pour GN (Sinnott, 2020, p.890; Frauss, 2015, p.51)
hf = 14/1000; % m Hauteur d'une ailette (Sinnott, 2020, p.890)
Df = (Dto+2*hf); % m Diam�tre ailette
ef = 0.5/1000; % m Epaisseur d'une ailette (Frass, 2015, p.51)
Nf = floor(1/Pf); % Nombre d'ailette par 1m
At = pi*Dto; % m2 Aire d'un tube par 1m;
Ate = At - Nf*ef*pi*Dto; % m2 Aire tube sans ailette
Af = Nf*(2*pi*(Df^2 - Dto^2)/4 + ef*pi*Df); % m2 Aire d'une ailette
Atot = Ate + Af;
Laech = 1; % m Largeur pos�e � 1m
Nrangees = 8; % Nombre de rang�es (Frass, 2015, p.53)

% --- Propri�t�s du biogaz et calculs interm�diaires --- %
% Masses molaires
M_CH4 = 16.04/1000; % kg/mol
M_CO2 = 44.01/1000; % kg/mol
M_gaz = y15_CH4*M_CH4 + y15_CO2*M_CO2; % kg/mol

% Masse volumique (hypoth�se de fluide incompressible)
rho_gaz = P*M_gaz/(R*Tm_gaz); % kg/m3

% Viscosit� dynamique
mu_gaz = 1.03e-5; % Pa*s (Gallagher, 2006, p.51)

% Conductivit� (Yaws, 2014, table 84 et 85)
k_CO2 = -6.07829291578281E-03 + 7.53509863170161E-05*Tm_gaz + ...
    9.49276579202504E-09*Tm_gaz^2 + -1.12751601050646E-11*Tm_gaz^3; % W/mK
k_CH4 = 5.37671485384522E-03 + 5.15550830467679E-05*Tm_gaz + ...
    1.66548977272723E-07*Tm_gaz^2 + -5.71678321678305E-11*Tm_gaz^3; % W/mK
k_gaz = 1/(y15_CH4/k_CH4 + y15_CO2/k_CO2); % W/mK (Gilmore & Comings, 1966)

% Capacit� thermique (Yaws, 2014, table 165)
Cp_CH4 = 46.241914499032198 + -0.15641465372*Tm_gaz + 0.00063340927021*Tm_gaz^2 ...
    + -9.8938057655e-007*Tm_gaz^3 + 8.9552744542e-010*Tm_gaz^4; % J/molK
Cp_CO2 = 22.870906792800799 + 0.053885041357*Tm_gaz + -3.0277722822e-006*Tm_gaz^2 ...
    + -5.9414717025e-008*Tm_gaz^3 + 5.758653574e-011*Tm_gaz^4; % J/molK
Cp_gaz = y15_CH4*Cp_CH4 + y15_CO2*Cp_CO2; % J/molK

% Nombre de Prandtl
Pr_gaz = mu_gaz*Cp_gaz/(M_gaz*k_gaz);

% D�bit volum�trique
Q_gaz = n15*M_gaz/rho_gaz; % m3/s 

% --- Propri�t�s de l'eau de refroidissement et calculs interm�diaires --- %
% Masse molaire
M_H2O = 0.018015; % kg/mol

% Masse volumque 
rho_H2O = 0.997; % kg/m3 � Tf = 300K (Bergman, 2018, Table A6) 

% Cp de l'eau liquide 
Cp_H2O = 4179; % J/KgK � Tf = 300K (Bergman, 2018, Table A6) 

% Conductivit�
k_H2O = 613e-3; % W/mK � Tf = 300K (Bergman, 2018, Table A6) 

% Viscosit� dynamique
mu_H2O = 855e-6; % Pa*s � Tf = 300K (Bergman, 2018, Table A6) 

% Nombre de Prandtl
Pr_H2O = 5.83; % � Tf = 300K (Bergman, 2018, Table A6) 

% D�bit volum�trique de l'eau


% --- Calcul du facteur de correction courant crois� --- %
Tci = T15;
Tco = T18;
Thi = T16;
Tho = T17;
R_ = (Thi - Tho)/(Tco - Tci);
P_ = (Tco - Tci)/(Thi - Tci);
F = 0.565; % Facteur de correction (Bergman, 2018, fig 11S.3)

% --- Calcul de la diff�rence moyenne de temp�rature logarithmique --- %
dT1 = Thi-Tco; 
dT2 = Tho-Tci;
dTlm = (dT2-dT1)/log(dT2/dT1);

% --- D�but des int�rations avec le U --- %
U = U0;
Ucal = 0;
tol = 0.01;

A = Q/(F*U*dTlm); % m2 Surface d'�change entre les fluides


while abs(Ucal - U) > tol
    % --- Calcul de l'aire et de rarrangement des tubes --- %
    A = Q/(F*U*dTlm); % m2 Surface d'�change entre les fluides
    Ntubes = floor(A/Atot); % Nombre de tube de 1m
    Ntpr = Ntubes/Nrangees; % Nombre de tube par rang�e
    Hech = Ntpr*Df; % Hauteur de l'�changeur avec ailettes coll�es
    Loech = Nrangees*sqrt(3)/2*(Dto+2*hf); % Longueur de l'�changeur
    
    % --- Calcul de la vitesse transversale du gaz --- %
    Atr = Dto*Laech*Ntpr + Nf*2*hf*ef*Laech*Ntpr; % Aire d'une reng�e de tube transversale
    Alibre = Hech*Laech - Atr; % Aire de passage du gaz
    v_gaz = Q_gaz/Alibre; % m/s 
    
    % --- Calcul du nombre de Reynolds du gaz --- %
    L_gaz = Dto*Atot/At;
    Re_gaz = rho_gaz*v_gaz*L_gaz/mu_gaz;
    Nu_gaz = 0.45*Re_gaz^0.625*Pr_gaz^(1/3)*(Atot/At)^-0.375; %(Frauss, 2015, p,12)
    ho = k_gaz*Nu_gaz/L_gaz; % W/m2K h c�t� gaz
    
    % Nombre de Reynolds
    Q_H2
    v_H2O = (n16*M_H2O/rho_H2O)/(pi*Dti^2/4); % m/s Vitesse de l'eau dans un tube
    Re_H2O = rho_H2O*v_H2O*Dti/mu_H2O

    
    
    break
    
end
