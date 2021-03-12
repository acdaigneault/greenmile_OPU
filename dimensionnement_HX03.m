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
Q = -1.443371916708512e+05; % W Selon les bilans d'énergie

% --- Informations sur les débits --- %
n15 = 70.576283371900020; % mol/s
y15_CH4 = 0.6243; 
y15_CO2 = 1 - y15_CH4;
n16 = 44.092918817800125; % mol/s

% --- Variable posée --- %
U0 = 150; % W/m2K U approximatif moyen pour transfert gaz-eau (Hall, 2018, p.221)

% --- Températures déterminées et pression --- %
T15 = 343.15; % K T opération AD-01 (Okoro & Sun, 2019)
T18 = 288.15; % K T opération AB-01, choix de l'équipe
T16 = 280.15; % K Moyenne température in de refroidissement (Hall, 2018, p.391)
T17 = 223.15; % K Température max vers tour de refroissement (Hall, 2018, p.394)
P = 800000; % Pa Pression d'opération de AB-01
R = 8.314; % m3*Pa/molK
Tm_gaz = (T15+T18)/2; % Température moyenne 

% --- Informations de conception pour les tubes --- %
N = 1; % Nombre de paire de tubes
Npass = 2*N; % Nombre de passages
Rfo = 0.0003; % m2C/W Facteur de Fouling max pour de l'eau de refroidissement (Sinnott, 2020, p. 780) 
Rfi = 0.0002; % m2C/W Facteur de Fouling pour les hydrocarbures légers (Sinnott, 2020, p. 780) 
Dto = 16/1000; % m Diamètre minimum standard pour des tubes en acier (Sinnott, 2020, p.786)
Dti = (16 - 2*1.7)/1000; % m Diamètre interne (choix arbitraire d'épaisseur selon littérature) (Sinnott, 2020, p.786)
Ltube = 2.4384; % m Longueur posée selon la littérature (Sinnott, 2020, p.785)
Pratio = 1.25; % Pitch ratio minimal (ratio de distance entre les tubes) (Sinnott, 2020, p.787)
config = 1; % Type de configuration (1: triangular, 2 : square)
 
 
% --- Propriétés du gaz naturel et intermédiaires --- %
% Masses molaires
M_CH4 = 16.04/1000; % kg/mol
M_CO2 = 44.01/1000; % kg/mol
M_gaz = y15_CH4*M_CH4 + y15_CO2*M_CO2; % kg/mol

% Masse volumique (hypothèse de fluide incompressible)
rho_gaz = P*M_gaz/(R*Tm_gaz); % kg/m3 

% Viscosité dynamique
mu_gaz = 1.03e-5; % Pa*s (Gallagher, 2006, p.51)

% Conductivité
k_CO2 = -6.07829291578281E-03 + 7.53509863170161E-05*Tm_gaz + ...
    9.49276579202504E-09*Tm_gaz^2 + -1.12751601050646E-11*Tm_gaz^3; % W/mK (Yaws, 2014, table 84)
k_CH4 = 5.37671485384522E-03 + 5.15550830467679E-05*Tm_gaz + ...
    1.66548977272723E-07*Tm_gaz^2 + -5.71678321678305E-11*Tm_gaz^3; % W/mK (Yaws, 2014, table 85)
k_gaz = 1/(y15_CH4/k_CH4 + y15_CO2/k_CO2); % W/mK (Gilmore & Comings, 1966)

% Capacité thermique
% Cp du CH4 (Table 165 [112 - 612]K)
A_CH4 = 46.241914499032198; B_CH4 = -0.15641465372;
C_CH4 = 0.00063340927021; D_CH4 = -9.8938057655e-007;
E_CH4 = 8.9552744542e-010; F_CH4 = -3.8921664729e-013;
Cp_CH4 = A_CH4 + B_CH4*Tm_gaz + C_CH4*Tm_gaz^2 + D_CH4*Tm_gaz^3 + E_CH4*Tm_gaz^4 + F_CH4*Tm_gaz^5; % J/molK

% Cp du CO2 (Table 165 [250 - 1100]K)
A_CO2 = 22.870906792800799; B_CO2 = 0.053885041357;
C_CO2 = -3.0277722822e-006; D_CO2 = -5.9414717025e-008;
E_CO2 = 5.758653574e-011; F_CO2 = -1.7570478487e-014;
Cp_CO2 = A_CO2 + B_CO2*Tm_gaz + C_CO2*Tm_gaz^2 + D_CO2*Tm_gaz^3 + E_CO2*Tm_gaz^4 + F_CO2*Tm_gaz^5; % J/molK

Cp_gaz = y15_CH4*Cp_CH4 + y15_CH4*Cp_CO2; % J/molK

Q_gaz = n15*M_gaz/rho_gaz;

% --- Attribution des paramètres --- %
% Températures (K)
Tci = T15;
Tco = T18;
Thi = T16;
Tho = T17;

% Facteur de sécurité arbitraire
Fsec = 0.9; 

% --- Calcul du facteur de correction (échangeur à contre-courant)
R = (Thi - Tho)/(Tco - Tci);
P = (Tco - Tci)/(Thi - Tci);
W = ((1 - P*R)/(1 - P))^(1/N);
S = sqrt(R^2 + 1)/(R - 1);
F = S*log(W)/log((1+W-S+S*W)/(1+W+S-S*W));

% --- Calcul de la différence moyenne de température logarithmique --- %
dT1 = Thi-Tco;
dT2 = Tho-Tci;
dTlm = (dT2-dT1)/log(dT2/dT1);

% --- Début des intérations avec le U --- %
U = U0;
Ucal = 0;
tol = 0.01;

while abs(Ucal - U) > tol
    % --- Calcul de l'aire et de la géométrie des tubes --- % (Hall, 2008, p.210-213)
    A = Q/(F*U*dTlm); % m2 Surface d'échange entre les fluides
    Ltubes = A*Fsec/(pi*Dto); % m Longueur totale des tubes 
    Ntube = ceil((Ltubes/(Ltube*Npass) + 0.5)*Npass); % Nombre de tubes
    
    % Aire requise pour un tube selon la configuration (1: triangular, 2 : square)
    if config == 1
        Atube = 2*(Pratio*Dto)^2*(sqrt(3)/4);
    elseif config == 2
        Atube = (Pratio*Dto)^2;
    else
        fprintf('Erreur pour le type de configuration de tubes')
    end
    
    % --- Calcul de l'aire et de la géométrie de la calandre --- %  (Hall, 2008, p.213)
    Dtight = 2*sqrt(Ntube*Atube/pi); % Aire totale pour les tubes 
    Acorr = Dtight*Dto*(Npass-1)+(Ntube*Atube); % Facteur de correction pour l'aire
    Dmin = 2*sqrt(Acorr/pi) + 2*Dto; % m Diamètre minimum pour la calandre
    D = ceil(Dmin*39.37)/39.37; % m Arrondis à une valeur de po
    
    % --- Calcul du coefficient de transfert côté tubes --- %
    Stubes = Ntube*pi*Dti^2/4; % m2 Aire totale de section dans les tubes
    v_gaz = Q_gaz/Stubes; % m/s vitesse du gaz
    Re = rho_gaz*v_gaz*Dti/mu_gaz;
    Pr = mu_gaz*Cp_gaz/k_gaz; 
    Nu = 0.021*Re^0.8*Pr^0.33;
    hi = Nu*k_gaz/Dti;

    
    

    
    


        break

end


% vérifier hypothèse fluide incompressible
    Stubes = Ntube*pi*Dti^2/4;
    v_gaz = Q_gaz/Stubes;