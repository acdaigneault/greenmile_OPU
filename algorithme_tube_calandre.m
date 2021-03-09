function [n16, Qech] = algorithme_tube_calandre(T,Q,U0,N)

% ------------------------------------------------------------------------
% Fonction 'bilan_energie_HX03' :
% Calcul du débit d'eau de refroidissement nécessaire pour refroidir le gaz
% et calcul de l'échangeur de chaleur réalisé entre les deux fluides.
% Argument d'entrée :
% - n : Débits entrants et sortants de l'échangeur de chaleur (mol/s)
% - T : Températures des courants reliés à l'échangeur de chaleur (K)
% Argument de sortie :
% - n16 : Débit d'eau de refroidissement (mol/s)
% - Qech : Transfert de chaleur (W)
% ------------------------------------------------------------------------

% --- Attribution des paramètres --- %
Tci = T(1);
Tco = T(2);
Thi = T(3);
Tho = T(4);

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

% --- Calcul de l'aire totale d'échangeur --- %
A = Q/(F*U*dTlm);


end