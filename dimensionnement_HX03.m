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

% --- Variable posée --- %
U0 = 150; % W/m2K U approximatif moyen pour transfert gaz-eau (Hall, 2018, p.221)

% --- Températures déterminées --- %
T15 = 343.15; % K T opération AD-01 (Okoro & Sun, 2019)
T18 = 288.15; % K T opération AB-01, choix de l'équipe
T16 = 280.15; % K Moyenne température in de refroidissement (Hall, 2018, p.391)
T17 = 223.15; % K Température max vers tour de refroissement (Hall, 2018, p.394)

% --- Informations de conception pour les tubes --- %
N = 1; % Nombre de paire de tubes
Rfo = 0.0003; % m2C/W Facteur de Fouling max pour de l'eau de refroidissement (Sinnott, 2020, p. 780) 
Rfi = 0.0002; % m2C/W Facteur de Fouling pour les hydrocarbures légers (Sinnott, 2020, p. 780) 
Dto = 16/1000; % m Diamètre minimum standard pour des tubes en acier (Sinnott, 2020, p.786)
Dti = (16 - 2*1.7)/1000; % m Diamètre interne (choix arbitraire d'épaisseur selon littérature) (Sinnott, 2020, p.786)
Ltube = 2.4384; % m Longueur posée selon la littérature (Sinnott, 2020, p.785)
Pratio = 1.25; % Pitch ratio minimal (ratio de distance entre les tubes) (Sinnott, 2020, p.787)
config = 1; % Type de configuration (1: triangular, 2 : square)
k = 

T = [T15 T18 T16 T17];
Rf = [Rfi Rfo];
param_tube = [N Dto Dti Ltube Pratio config];

%%
% --- Attribution des paramètres --- %
% Températures (K)
Tci = T(1);
Tco = T(2);
Thi = T(3);
Tho = T(4);

% Facteurs d'encrassement (m2C/W)
Rfi = Rf(1); 
Rfo = Rf(2);

% Informations sur les tubes
N = param_tube(1); 
Dto = param_tube(2); % m
Dti = param_tube(3); % m
Ltube = param_tube(4); % m
Pratio = param_tube(5); 
config = param_tube(6); 
Npass = 2*N; % Nombre de passages

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
    
    

    
    


        break

end
