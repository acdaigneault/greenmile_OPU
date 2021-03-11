function [n16, Qech] = algorithme_tube_calandre(T,Q,U0,N)

% ------------------------------------------------------------------------
% Fonction 'bilan_energie_HX03' :
% Calcul du d�bit d'eau de refroidissement n�cessaire pour refroidir le gaz
% et calcul de l'�changeur de chaleur r�alis� entre les deux fluides.
% Argument d'entr�e :
% - n : D�bits entrants et sortants de l'�changeur de chaleur (mol/s)
% - T : Temp�ratures des courants reli�s � l'�changeur de chaleur (K)
% Argument de sortie :
% - n16 : D�bit d'eau de refroidissement (mol/s)
% - Qech : Transfert de chaleur (W)
% ------------------------------------------------------------------------

% --- Attribution des param�tres --- %
% Temp�ratures (K)
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

% Facteur de s�curit� arbitraire
Fsec = 0.9; 

% --- Calcul du facteur de correction (�changeur � contre-courant)
R = (Thi - Tho)/(Tco - Tci);
P = (Tco - Tci)/(Thi - Tci);
W = ((1 - P*R)/(1 - P))^(1/N);
S = sqrt(R^2 + 1)/(R - 1);
F = S*log(W)/log((1+W-S+S*W)/(1+W+S-S*W));

% --- Calcul de la diff�rence moyenne de temp�rature logarithmique --- %
dT1 = Thi-Tco;
dT2 = Tho-Tci;
dTlm = (dT2-dT1)/log(dT2/dT1);

% --- D�but des int�rations avec le U --- %
U = U0;
Ucal = 0;
tol = 0.01;

while abs(Ucal - U) < tol
    % --- Calcul de l'aire et de la g�om�trie des tubes --- % (Hall, 2008, p.210-213)
    A = Q/(F*U*dTlm); % m2 Surface d'�change entre les fluides
    Ltubes = A*Fsec/(pi*Dto); % m Longueur totale des tubes 
    Ntube = ceil((Ltubes/(Ltube*Npass) + 0.5)*Npass); % Nombre de tubes
    if strcmp(config,'triangular')
        Atube = 2*(Pratio*Dto)^2*(sqrt(3)/4); % Aire requise pour un tube 
    elseif strcmp(config,'square')
        Atube = (Pratio*Dto)^2;
    else
        fprintf('Erreur pour le type de configuration de tubes')
    end
    
    % --- Calcul de l'aire et de la g�om�trie de la calandre --- %  (Hall, 2008, p.213)
    Dtight = 2*sqrt(Ntube*Atube/pi); % Aire totale pour les tubes 
    Acorr = Dtight*Dto*(Npass-1)+(Ntube*Atube); % Facteur de correction pour l'aire
    Dmin = 2*sqrt(Acorr/pi) + 2*Dto; % Diam�tre minimum pour la calandre
    
    
    
    
    




end