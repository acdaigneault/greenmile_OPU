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
Tci = T(1);
Tco = T(2);
Thi = T(3);
Tho = T(4);

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

% --- Calcul de l'aire totale d'�changeur --- %
A = Q/(F*U*dTlm);


end