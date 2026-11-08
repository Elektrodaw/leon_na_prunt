% clear all; close all; clc;
% 
% %% Dane 
% m_s = 950;                          % [kg] masa wlasna
% dw = (16*2.54)/100 + (205*0.55*2)/1000; % [m] srednica kola
% m_v = m_s + 70;                     % [kg] masa z kierowca
% 
% P_n;                                % [kW] moc silnika
% i;                                  % przekladnia
% 
% jtm = 0.05/100;                     % [kgm^2/kW] wspolczynnik z (?)
% m_1 = P_n*jtm*(4*i^2)/(dw^2);       % [kg] masa od momentu 
% m_2 = m_v + (P_n/1000)*0.5;         % [kg] masa z kierowca i silnikiem
% m_3 = 100;                          % [kg] masa osprzetu kół
% 
% aku;                                % [kWh] pojemnosc aku
% m_aku = 8*aku;                      % [kg] masa aku
% 
% P_ac = 0.03 * aku + 0.05 * P_n + 1; % [kW] moc osprzetu
% 
% m = m_1+ m_2+ m_3+ m_aku;           % [kg] masa calkowita

function [J4]=leon(x)

% wartości zadane
t_100_opt = 5; % [s]
t_1_4_opt = 13; % [s]
V_max_opt = 190; % [km/h]
V_max_opt = V_max_opt/3.6; % [m/s]
S_max = 300; % [km]

%zmienne wejściowe
Pn = x(1);
i = x(2);

d_k = (16*2.54)/100 + (205*0.55*2)/1000; % [m] srednica kola
m1 = JJ1=Pn/100000*0.05; %mom. bezwładnosci: 0.05kgm2 na każde 100kW

%masy [kg]
m1=JJ1*4*i^2/(d_k^2);
m2 = mv + Pn/1000*0.5 + 8*aku;
m3 = 100;
