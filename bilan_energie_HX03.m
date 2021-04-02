function [n16, Qech] = bilan_energie_HX03(n,T)

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

% --- Attribution des paramètres en variables --- %
n15_CH4 = n(1); % débit CH4 (mol/s)
n15_CO2 = n(2); % débit CO2 (mol/s)

T1 = T(1); % Tin courant à refroidir
T2 = T(2); % Tout courant à refroidir
T3 = T(3); % Tin refroidissant
T4 = T(4); % Tout refroidissant

T0 = 298.15;

% --- Capacitances thermiques (Yaws sur Knovel) --- %
% Cp du CH4 (Table 165 [112 - 612]K)
A_CH4 = 46.241914499032198; B_CH4 = -0.15641465372;
C_CH4 = 0.00063340927021; D_CH4 = -9.8938057655e-007;
E_CH4 = 8.9552744542e-010; F_CH4 = -3.8921664729e-013;

% Cp du CO2 (Table 165 [250 - 1100]K)
A_CO2 = 22.870906792800799; B_CO2 = 0.053885041357;
C_CO2 = -3.0277722822e-006; D_CO2 = -5.9414717025e-008;
E_CO2 = 5.758653574e-011; F_CO2 = -1.7570478487e-014;

% Cp de l'eau liquide (Table 162 [273.16 - 620]K)
A_H2Oliq = -3385.59633521043; B_H2Oliq = 44.119996293;
C_H2Oliq = 	-0.22135593182; D_H2Oliq = 0.00054657156138;
E_H2Oliq = -6.6489932771e-007; F_H2Oliq = 3.1947490101e-010;

Cp_CH4 = @(T) A_CH4 + B_CH4*T + C_CH4*T.^2 + D_CH4*T.^3 + E_CH4*T.^4 + F_CH4*T.^5; % J/molK
Cp_CO2 = @(T) A_CO2 + B_CO2*T + C_CO2*T.^2 + D_CO2*T.^3 + E_CO2*T.^4 + F_CO2*T.^5; % J/molK
Cp_H2Oliq = @(T) A_H2Oliq + B_H2Oliq*T + C_H2Oliq*T.^2 + D_H2Oliq*T.^3 + E_H2Oliq*T.^4 + F_H2Oliq*T.^5; % J/molK

% --- Calcul des enthalpies (Ref CH4(g), CO2(g), H2O(liq) @ T6, 0.8Mpa) --- %
H1 = integral(Cp_CH4,T0,T1); % CH4 in (J/mol)
H2 = integral(Cp_CO2,T0,T1); % CO2 in (J/mol)
H3 = integral(Cp_H2Oliq,T0,T3); % H2O in (J/mol)
H4 = integral(Cp_CH4,T0,T2); % CH4 out (J/mol)
H5 = integral(Cp_CO2,T0,T2); % CO2 out (J/mol)
H6 = integral(Cp_H2Oliq,T0,T4); % H2O out (J/mol)

% --- Calcul du débit d'eau froide nécessaire (Échangeur de chaleur adiabatique) --- %
% Q = dH = 0, en isolant n3 (débit d'eau)
n16 = (n15_CH4*(H4 - H1) + n15_CO2*(H5 - H2))/(H3 - H6); % mol/s

% --- Calcul Q de l'échangeur de chaleur --- %
Qech = n15_CH4*(H4 - H1) + n15_CO2*(H5 - H2); % W


end