function [n23, Qech] = bilan_energie_F01(n,T)

% ------------------------------------------------------------------------
% Fonction 'bilan_energie_F01' :
% Calcul du débit d'eau à vaporiser et  calcul de l'échangeur de chaleur 
% réalisé entre les deux fluides dans la fournaise.
% Argument d'entrée :
% - n : Débits entrants et sortants de la fournaise (mol/s)
% - T : Températures des courants reliés à la fournaise (K)
% Argument de sortie :
% - n23 : Débit d'eau à vaporiser (mol/s)
% - Qech : Transfert de chaleur (W)
% ------------------------------------------------------------------------

% --- Attribution des paramètres en variables --- %
n21_CH4 = n(1); % débit CH4 in (mol/s)
n21_CO2 = n(2); % débit CO2 in (mol/s)
n25_CO2 = n(3); % débit CO2 out (mol/s)
n25_O2 = n(4); % débit O2 out (mol/s)
n25_H2O = n(5); % débit H2O out (mol/s)
n22_N2 = n(6); % débit N2 (mol/s)
n22_O2 = n(7); % débit O2 in (mol/s)

T21 = T(1); % Tin gaz
T25 = T(2); % Tout fumée
T23 = T(3); % Tin eau à vaporiser
T24 = T(4); % Tout vapeur saturée

eps = n21_CH4; % Combustion parfait => eps = débit CH4

T22 = 298.15; % Température posée de l'air sec (25C)
T0 = T22; % Température de référence (K)
% --- Capacitances thermiques (Yaws sur Knovel) --- %

% Cp du CH4 (Table 165 [112 - 612]K)
A_CH4 = 46.241914499032198; B_CH4 = -0.15641465372;
C_CH4 = 0.00063340927021; D_CH4 = -9.8938057655e-007;
E_CH4 = 8.9552744542e-010; F_CH4 = -3.8921664729e-013;

% Cp du CO2 (Table 165 [250 - 1100]K)
A_CO2 = 22.870906792800799; B_CO2 = 0.053885041357;
C_CO2 = -3.0277722822e-006; D_CO2 = -5.9414717025e-008;
E_CO2 = 5.758653574e-011; F_CO2 = -1.7570478487e-014;

% Cp de l'eau en vapeur saturée (Table 164 [273.16 - 610]K)
A_H2Ovap = -5199.3274153812099; B_H2Ovap = 67.020972634;
C_H2Ovap = 	-0.33821124538; D_H2Ovap = 0.00084089969775;
E_H2Ovap = -1.0311971019e-006; F_H2Ovap = 5.0051424738e-010;

% Cp de l'eau liquide (Table 162 [273.16 - 620]K)
A_H2Oliq = -3385.59633521043; B_H2Oliq = 44.119996293;
C_H2Oliq = 	-0.22135593182; D_H2Oliq = 0.00054657156138;
E_H2Oliq = -6.6489932771e-007; F_H2Oliq = 3.1947490101e-010;

% Cp du O2 gazeux (Table 164 [100 - 1000]K)
A_O2 = 31.901200000221301; B_O2 = -0.027541317952;
C_O2 = 	8.7685722629e-005; D_O2 = -8.6092599107e-008;
E_O2 = 2.8884615423e-011; F_O2 = 3.8461524933e-017;

% Cp du N2 gazeux (Table 164 [100 - 2000]K)
A_N2 = 31.328191279820899; B_N2 = -0.016958275935;
C_N2 = 	4.0270675931e-005; D_N2 = -3.1547389219e-008;
E_N2 = 1.1068014938e-011; F_N2 = -1.4759097237e-015;

Cp_CH4 = @(T) A_CH4 + B_CH4*T + C_CH4*T.^2 + D_CH4*T.^3 + E_CH4*T.^4 + F_CH4*T.^5;
Cp_CO2 = @(T) A_CO2 + B_CO2*T + C_CO2*T.^2 + D_CO2*T.^3 + E_CO2*T.^4 + F_CO2*T.^5;
Cp_H2Oliq = @(T) A_H2Oliq + B_H2Oliq*T + C_H2Oliq*T.^2 + D_H2Oliq*T.^3 + E_H2Oliq*T.^4 + F_H2Oliq*T.^5;
Cp_H2Ovap = @(T) A_H2Ovap + B_H2Ovap*T + C_H2Ovap*T.^2 + D_H2Ovap*T.^3 + E_H2Ovap*T.^4 + F_H2Ovap*T.^5;
Cp_O2 = @(T) A_O2 + B_O2*T + C_O2*T.^2 + D_O2*T.^3 + E_O2*T.^4 + F_O2*T.^5;
Cp_N2 = @(T) A_N2 + B_N2*T + C_N2*T.^2 + D_N2*T.^3 + E_N2*T.^4 + F_N2*T.^5;

% --- Informations sur les enthalpies --- %
Hvap_eau = 40.744*1000; % J/mol (table 63)
Hf_CH4 = -74.5*1000; % J/mol (table 57)
Hf_O2 = 0*1000; % J/mol (table 56)
Hf_H2O = -241.83*1000; % J/mol (table 56)
Hf_CO2 = -393.5*1000; % J/mol (table 57)

Hr = (2*Hf_H2O + Hf_CO2) - (Hf_CH4 + 2*Hf_O2); % J/mol

% --- Calcul des enthalpies (Ref CH4(g), CO2(g), H2O(liq), N2(g), O2(g) @25C, Pin) --- %
H1 = integral(Cp_CH4,T0,T21); % CH4 in (J/mol)
H2 = integral(Cp_CO2,T0,T21); % CO2 in (J/mol)
H3 = integral(Cp_N2,T0,T22); % N2 in (J/mol)
H4 = integral(Cp_O2,T0,T22); % O2 out (J/mol)
H5 = integral(Cp_H2Oliq,T0,T23); % H2O liq out (J/mol)

H6 = integral(Cp_CO2,T0,T25); % CO2 out (J/mol)
H7 = integral(Cp_H2Oliq,T0,100+273.15) + Hvap_eau + integral(Cp_H2Ovap,100+273.15,T25); % H2Orx in (J/mol)
H8 = integral(Cp_N2,T0,T25); % N2 out (J/mol)
H9 = integral(Cp_O2,T0,T25); % O2 out (J/mol)
H10 = integral(Cp_H2Oliq,T0,100+273.15) + Hvap_eau + integral(Cp_H2Ovap,100+273.15,T24); % H2Ovap out (J/mol)

% --- Calcul du débit d'eau froide nécessaire (Échangeur de chaleur adiabatique) --- %
% Q = dH = 0, en isolant n23 (débit d'eau)
n23 = (eps*Hr + n25_CO2*H6 + n25_H2O*H7 + n22_N2*(H8-H3) + n25_O2*H9 - n21_CH4*H1 - n21_CO2*H2 - n22_O2*H4)/(H10 - H5); % mol/s

% --- Calcul Q de l'échangeur de chaleur --- %
Qech =  n23*(H10 - H5);% W


end