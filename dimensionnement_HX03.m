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
Q = -8.458080117842977e+04; % A changer!

% --- Températures déterminées --- %
T15 = 343.15; % K T opération AD-01 (Okoro & Sun, 2019)
T18 = 298.15; % K A VALIDER (T opération AB-01)
T16 = 280.15; % K Moyenne température in de refroidissement (Hall, 2018, p.391)
T17 = 285.15; % K < que la température max vers tour de refroissement (Hall, 2018, p.394)

% --- Autres informations de conception --- %
N = 1; % Nombre de paire de tubes
Rfo = 0.0003; % m2C/W Facteur de Fouling max pour de l'eau de refroidissement (Sinnott, 2020, p. 780) 
Rfi = 0.0002; % m2C/W Facteur de Fouling pour les hydrocarbures légers (Sinnott, 2020, p. 780) 
Dto = 16/1000; % m Diamètre minimum s+ standard pour des tubes en acier (Sinnott, 2020, p.786)
Dti = (16 - 2*1.7)/1000; % m Diamètre interne (choix arbitraire d'épaisseur selon littérature) (Sinnott, 2020, p.786)



% --- Variable posée --- %
U0 = 150; % W/m2K U approximatif moyen pour transfert gaz-eau (Hall, 2018, p.221)

