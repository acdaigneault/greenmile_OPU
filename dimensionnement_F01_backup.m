% -----------------------------------------------------------------------
% Nom du fichier : dimensionnement_F01.m
% Description : Calculs concernant le dimensionnement de la fournaise F-01
% Auteurs : GREEN MILE (Audrey Collard-Daigneault)
% Date : avril 2021
% -----------------------------------------------------------------------

clc
clear all

%% --- Informations obtenues par le bilan de matière --- %%
nin_CH4 = 1.2061; % mol/s Débit molaire de CH4 entrant dans F-01
nin_CO2 = 0.0373; % mol/s Débit molaire de CO2 entrant dans F-01

nout_CO2 = 1.2434; % mol/s Débit molaire de CO2 sortant de F-01
nout_H2O = 2.4123; % mol/s Débit molaire de H2O sortant de F-01
nout_N2 = 9.9823; % mol/s Débit molaire de N2 sortant de F-01
nout_O2 = 0.2412; % mol/s Débit molaire de O2 sortant de F-01
nout = nout_CO2 + nout_H2O + nout_N2 + nout_O2; % mol/s Débit molaire sortant

Tin_gaz = 15 + 273.15; %K (T opération AB-01)
Tout_gaz = 150 + 273.15; %K (T gaz sortant 150C)(Song & al. 2004)
Tin_H2O = 25 + 273.15; %K (Eau à 35C minimum (T opération RX-01))
Tout_H2O = 100 + 273.15; %K (Vapeur saturée à 1 atm)

%% --- Puissance introduite --- %%
Ting = 905; % K température minimum d'ignition du CH4 (Perry & Green, 2008, p. 24-13)
Hc_CH4 = 892.68*1000; % J/mol Enthalpie de combustion(table 60)
Hc_O2 = (-8.4165335332647153 + 0.026684484417*Ting + 5.5210633049e-006*Ting^2 +...
    -1.17125478e-009*Ting^3 + 1.026e-013*Ting^4)*1000; % J/mol Enthalpie de combustion (table 76)
Hc_H2O = (-9.6702911285054416 + 0.030947811819*Ting + 4.9273659307e-006*Ting^2 +...
    2.2142365e-010*Ting^3 + -8.579e-014*Ting^4)*1000; % J/mol Enthalpie de combustion (table 76)
Hc_CO2 = 0.02*1000; % J/mol Enthalpie de combustion (table 60)

Hr = abs((2*Hc_H2O + Hc_CO2) - (Hc_CH4 + 2*Hc_O2)); % J/mol % Loi de Hess
Pi = nin_CH4*Hr; % W Puissance introduite

%% --- Rendement et puissance utile --- %%
n = 0.8; % rendement de la fournaise (Sinnott, 2020, p.906)
Pu = n*Pi; % W Puissance utile transmise aux tubes

%% --- Section radiante --- %%
sigma = 5.67e-8; % W/m2K4
transfert_rad = 0.7; % 70% transféré dans la partie radiante (Sinnott, 2020, p.609
qrad = transfert_rad*Pu; % Énergie dans les tubes radiants

% Propriétés des tubes
Dto = 0.075; % m Diamètre externe (Sinnott, 2020, p.901)
Dti = Dto - 2*0.002; % m Diamètre interne
k_ss = 14.9; % W/mK Conductivité du stainless steel (Bergman, 2018)
esp = 0.1*Dto; % m Espace entre 2 tubes

Tgaz = 1200 + 273.15; % K (Song & al., 2004)
Tsat = 450; % K % vapeur saturé à 130.5 psi (Bergman, 2018)

Ts0 = Tsat;
Ts = 0;
while abs(Ts0-Ts) > 0.1
    Ts = Ts0;
    
    % Aire de la section radiante (alpha = 1 et F = 1)
    Arad = qrad /(sigma*(Tgaz^4 - Ts^4));
    
    % Calcul des dimensions avec des intervalles.
    W = 1;
    L = 2*W;
    H = (Arad - 2*W*L)/(2*L+W);
    
    while (abs(H/W) < 1 || H/W > 1.5)
        W = W - 0.01;
        L = 2*W;
        H = (Arad - 2*W*L)/(2*L+W);
    end
    
    Vrad = W*L*H; % m3 Volume radiant
    
    Ltrad = W*L; % m Longueur d'un tube qui contourne la fournaise
    nbtrad = floor((2*H+W)/(Dto+esp)); % Nombre de tubes dans le foyer radiant
    
    % --- Calcul du débit massique d'eau à réchauffer --- %
    Cpliq = 4400; % J/kgK Cp du liquide à 450K
    hfg = 2024000; % J/kg Enthalpie de vaporisation
    Xrad = 0.7; % Fraction de gaz
    m_eau = qrad/(Cpliq*(Tout_H2O-Tin_gaz) + hfg*Xrad);
    
    % Coefficient de transfert par convection du liquide
    mu_liq = 152e-6; % Pa*s Viscosité dynamique
    k_liq = 678e-3; % W/mK
    rho_liq = 1/1.123e-3;
    Pr_liq = 0.99;
    Re_liq = 4*m_eau/(mu_liq*Dti*pi);
    f_liq = (0.790*log(Re_liq)-1.64)^-2; % (Bergman, 2018, eq8.21)
    hsp = k_liq/Dti * ((f_liq/8)*(Re_liq-1000)*Pr_liq)/(1+12.7*sqrt(f_liq/8)*(Pr_liq^(2/3)-1));
    
    % Coefficient de transfert par convection du fluide en vaporisation
    mpp = (m_eau/nbtrad)/(pi*Dti^2/4); % kg/m2s mass flow rate per unit cross-sectional area
    %v_vap = 4*m_eau/(rho_liq*Dti^2*pi); % m/s Vitesse (vitesse du gaz = vitesse du liquide...)
    rho_vap = 1/0.208; % kg/m3
    Fr = (mpp/rho_liq)^2/(9.81*Dti); % Nombre de Froude (vitesse du gaz = vitesse du liquide)
    fFr = 2.63*Fr^0.3; % f(Fr)
    Gsf = 1; % Because stainless steel (Bergman, 2018, p.673)
    qconvpp = qconv/Aconv;
    
    X = Xrad;
    hi1 = hsp * 1.136 * (rho_liq/rho_vap)^0.45 * X^0.72 * (1-X)^0.08 * fFr + 667.2 * (qradpp/(mpp*hfg))^0.7*(1-X)^0.8*Gsf;
    hi2 = hsp * 0.6683 * (rho_liq/rho_vap)^0.1 * X^0.16 * (1-X)^0.64 * fFr + 1058 * (qradpp/(mpp*hfg))^0.7*(1-X)^0.8*Gsf;
    hi = max(hi1,hi2); % W/m2K
    
    % Calcul des résistances de conduction et convection
    Rcond = log(Dto/Dti)/(2*pi*nbtrad*Ltrad*k_ss); % Résistance thermique conduction
    Rconvi = 1/(hi1*(nbtrad*pi*Dti*Ltrad)); % Réstance thermique convection dans les tubes
    
    % Recalcul de la température à la surface extérieure
    Ts0 = (qrad*(Rcond+Rconvi)+Tsat);
end

%% --- Section convective --- %%
qconv_tubes = Pu - (qrad + 0.02*Pi); % W Énergie restante après transfert radiant

% Coefficient des Cp des gaz de combustion
% Cp du CO2 (Table 165 [250 - 1100]K)
A_CO2 = 22.870906792800799; B_CO2 = 0.053885041357;
C_CO2 = -3.0277722822e-006; D_CO2 = -5.9414717025e-008;
E_CO2 = 5.758653574e-011; F_CO2 = -1.7570478487e-014;

% Cp de l'eau en vapeur saturée (Table 39 [150 - 1500]K)
A_H2Ovap = 33.1743818702509; B_H2Ovap = -3.24633355028777E-03;
C_H2Ovap = 	1.74365275243944E-05; D_H2Ovap = -5.97957584435752E-09;

% Cp du O2 gazeux (Table 164 [100 - 1000]K)
A_O2 = 31.901200000221301; B_O2 = -0.027541317952;
C_O2 = 	8.7685722629e-005; D_O2 = -8.6092599107e-008;
E_O2 = 2.8884615423e-011; F_O2 = 3.8461524933e-017;

% Cp du N2 gazeux (Table 164 [100 - 2000]K)
A_N2 = 31.328191279820899; B_N2 = -0.016958275935;
C_N2 = 	4.0270675931e-005; D_N2 = -3.1547389219e-008;
E_N2 = 1.1068014938e-011; F_N2 = -1.4759097237e-015;

Cp_CO2 = @(T) A_CO2 + B_CO2*T + C_CO2*T.^2 + D_CO2*T.^3 + E_CO2*T.^4 + F_CO2*T.^5;
Cp_H2Ovap = @(T) A_H2Ovap + B_H2Ovap*T + C_H2Ovap*T.^2 + D_H2Ovap*T.^3;
Cp_O2 = @(T) A_O2 + B_O2*T + C_O2*T.^2 + D_O2*T.^3 + E_O2*T.^4 + F_O2*T.^5;
Cp_N2 = @(T) A_N2 + B_N2*T + C_N2*T.^2 + D_N2*T.^3 + E_N2*T.^4 + F_N2*T.^5;

% Itérations pour trouver la température de sortie du foyer radiant/entrée de la section convective
Tflamme = 2191; % K Température de flamme du méthane (Perry & Green, 2008, p.24-23)
Tconv_in = 275;
qrad_cal = 0;
while abs((qrad + 0.02*Pi) - qrad_cal)/qrad > 0.01  && Tconv_in < Tgaz
    Tconv_in = Tconv_in + 5; % Trouve la température Tconv de 275 K par bonds de 5 K
    qrad_cal = abs(nout_CO2*integral(Cp_CO2,Tflamme,Tconv_in) + ...
        nout_H2O*integral(Cp_H2Ovap,Tflamme,Tconv_in) + ...
        nout_O2*integral(Cp_O2,Tflamme,Tconv_in) +...
        nout_N2*integral(Cp_N2,Tflamme,Tconv_in));
end

% Aire de la section convective (vérification des X)
Xconv = 1; % La vapeur est 100% gaz à la sortie de la section convective
qconv = m_eau*hfg*(Xconv-Xrad);
qsurc = qconv_tubes - qconv; % W Énergie restante pour le surchauffeur

% Itérations pour la température de sortie de la section convective
Tconv_out = 275;
qconv_cal = 0;
while abs(qconv - qconv_cal)/qrad > 0.01  && Tconv_out < Tconv_in
    Tconv_out = Tconv_out + 5; % Trouve la température Tconv out de 275 K par bonds de 5 K
    qconv_cal = abs(nout_CO2*integral(Cp_CO2,Tconv_in,Tconv_out) + ...
        nout_H2O*integral(Cp_H2Ovap,Tconv_in,Tconv_out) + ...
        nout_O2*integral(Cp_O2,Tconv_in,Tconv_out) +...
        nout_N2*integral(Cp_N2,Tconv_in,Tconv_out));
end

% Calcul de l'efficacité du transfert de chaleur par la méthode NTU
Tmconv_gaz = (Tconv_in + Tconv_out)/2; % K % Temperature de gaz moyen
Cconv_min = nout_CO2*Cp_CO2(Tmconv_gaz) +  nout_H2O*Cp_H2Ovap(Tmconv_gaz) + ...
    nout_O2*Cp_O2(Tmconv_gaz) + nout_N2*Cp_N2(Tmconv_gaz);

econv = qconv/(Cconv_min*(Tconv_in-Tsat)); % Efficacité de la section convective
NTUconv = -log(1-econv); % Calcul NTU pour Cr = 0 (vaporisation) (Bergman, 2018, p725)

% Propriété des fluides 
% Flue gas @ Tmconv_gaz (Increase Performance, 2020)
M_gaz = 27.78/1000; % kg/mol
kconv_gaz = 0.054654; % W/mK 
Cpconv_gaz = 1233; % J/kgK
muconv_gaz = 0.0000337; % Pa*s
rhoconv_gaz = 0.4389; % kg/m3
Prconv_gaz = muconv_gaz*Cpconv_gaz/kconv_gaz;

% Eau en vapeur (moyenne)
Xm = (Xconv + Xrad)/2; % moyenne des fraction liq-gaz
M_H2O = 18.02/1000; % kg/mol
kconv_H2O = 1/((1-Xm)/678e-3 + Xm/33.1e-3); % W/m2K
Cpconv_H2O = (1-Xm)*4400 + Xm*2560; % J/kgK
muconv_H2O = 1/((1-Xm)/152e-6 + Xm/14.85e-6); % Pa*s
rhoconv_H2O = (1-Xm)*(1/1.123) + Xm*(1/0.208);  % kg/m3
Prconv_vap = 2560*14.85e-6/33.1e-3;
Prconv_liq = 4400*152e-6/678e-3;
Prconv_H2O = Xm*Prconv_vap + (1-Xm)*Prconv_liq;

% Débit volumétrique et vitesse d'entrée dans la section radiante
P = 101325 + 25; % Pa (pression dans fournaise posée)
Vconv_gaz = nout*8.314*Tconv_in/P;





 

    % --- Calcul de l'aire et de rarrangement des tubes --- %
    Aconvo = Cconv_min*NTUconv/Uconv; % Aire des tubes
    Ltconv = L; % m Longueur des tubes = longeueur fournaise
    nbtconv = floor(Aconvo/(pi*Dto*Ltconv)); % Nopmbre de tubes
    Wconv = W/2;
    nbtconv_rangee = 3; % Nombre de tube par rangée
    
    % --- Calcul du coefficient de transfert thermique côté gaz --- %
    aire_libre = Ltconv*(Wconv - nbtconv_rangee*Dto); %m2
    vconv_max = Vconv_gaz/aire_libre; % m/s 
    Reconv_gaz = rhoconv_gaz*vconv_max*Dto/muconv_gaz; % Nombre de Reynolds
    Tconvs = (Tmconv_gaz + Tsat)/2; % K T surface évaluée à la moyenne de Tsat et T moyenne de gaz
    kconvs_gaz = 0.054654; % W/mK 
    Cpconvs_gaz = 1186; % J/kgK Cp surface
    muconvs_gaz = 0.0000287; % Pa*s viscosité surface
    kconvs_gaz = 0.045315; % W/m2K Conductivité surface
    Prconvs_gaz = muconvs_gaz*Cpconvs_gaz/kconvs_gaz;
    C1 = 0.35; % 
    C2 = 0.95; % Nb rangée = 7
    m = 0.6;
    Nuconv_gaz = C2*C1*Reconv_gaz^m*Prconv_gaz^0.36*(Prconv_gaz/Prconvs_gaz)^(1/4);
    hconv_o = Nuconv_gaz*kconv_gaz/Dto; % W/m2K 
    
    % --- Calcul du coefficient de transfert thermique côté eau --- %
    Vconv_H2O = m_eau/rhoconv_H2O;
    vconv = Vconv_H2O/(nbtconv*pi*Dti^2/4); 
    Reconv_H2O = rhoconv_H2O*vconv*Dti/muconv_H2O; 
    Nuconv_H2O = 0.023*Reconv_H2O^(4/5)*Prconv_H2O^0.3; % (Bergman, 2018, p.519)
    hconv_i = kconv_H2O*Nuconv_H2O/Dti;
    
    X = Xconv - Xrad;
    hi1 = hsp * 1.136 * (rho_liq/rho_vap)^0.45 * X^0.72 * (1-X)^0.08 * fFr + 667.2 * (qradpp/(mpp*hfg))^0.7*(1-X)^0.8*Gsf;
    hi2 = hsp * 0.6683 * (rho_liq/rho_vap)^0.1 * X^0.16 * (1-X)^0.64 * fFr + 1058 * (qradpp/(mpp*hfg))^0.7*(1-X)^0.8*Gsf;
    hi = max(hi1,hi2); % W/m2K
    
    Rconv_cond = log(Dto/Dti)/(2*pi*nbtrad*Ltconv*k_ss); % Résistance thermique conduction
    Rconv_convi = 1/(hi1*(nbtrad*pi*Dti*Ltrad)); % Réstance thermique convection dans les tubes
    
                                                                                                                                                                                                           
    
%     Nu_gaz = 0.134*Re_gaz^0.681*Pr_gaz^0.33*((pf-ef)/hf)^0.2*(pf/ef)^0.1134; % (Bermann
%     ho = k_gaz*Nu_gaz/L_gaz; % W/m2K h côté gaz
    
    
    

%end


