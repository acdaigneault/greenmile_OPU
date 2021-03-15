% -----------------------------------------------------------------------
% Nom du fichier : reseau_AD01_AB01.m
% Description : Calculs concernant le dimensionnement du réseau entre le
%               AD-01 et AB-01 avec les équipements hydrauliques reliés.
% Auteurs : GREEN MILE (Audrey Collard-Daigneault)
% Date : mars 2021
% -----------------------------------------------------------------------
clc
clear all

% Informations sur le fluides
n = 70.576283371900020; % mol/s Débit molaire de biogaz
y_CH4 = 0.6243; % Fraction molaire de CH4
y_CO2 = 1 - y_CH4; % Fraction molaire de CO2