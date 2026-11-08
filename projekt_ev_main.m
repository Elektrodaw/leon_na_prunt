clear all; close all; clc;

%%  PARAMETRY DOCELOWE
p.t100opt = 5.0;         % Docelowy czas przyspieszenia 0-100 km/h [s]
p.t14opt  = 13.0;        % Docelowy czas wyścigu na 1/4 mili [s]
p.voptmax = 190 / 3.6;   % Docelowa prędkość maksymalna [m/s] (190 km/h)
p.sp      = 300 * 1000;  % Wymagany zasięg autostradowy [m] (300 km)
p.vp      = 150 / 3.6;   % Prędkość autostradowa dla zasięgu [m/s] (150 km/h)

%%  PARAMETRY BAZOWE LEONA
p.masa_ice_z_kierowca = 950 + 80; % Masa bazowa auta spalinowego z kierowcą [kg]
p.d_k   = (16*25.4 + 2*205*0.55) / 1000; % Średnica koła [m] (dla opon 205/45R17)
p.rho   = 1.225;              % Gęstość powietrza [kg/m^3]
p.Cd    = 0.32;                % Współczynnik oporu powietrza
p.A     = 2.1;                % Powierzchnia czołowa pojazdu [m^2]
p.g     = 9.81;               % Przyspieszenie ziemskie [m/s^2]

%%  PARAMETRY SILNIKA ELEKTRYCZNEGO (PMSM) I INNE STAŁE
p.wI    = 400;                % Koniec strefy stałego momentu obrotowego [rad/s]
p.wIII  = 1200;               % Maksymalna prędkość kątowa silnika [rad/s]
p.m3    = 100;                % Stała masa wirujących elementów kół [kg]
p.P_aux_stale = 1.0;          % Moc systemów niezwiązanych z napędem [kW]

%%  WARUNKI POCZĄTKOWE DLA ALGORYTMU OPTYMALIZACJI (x0)
Mmax0 = 350;        % Startowy maksymalny moment obrotowy [Nm]
przekladnia0 = 8.0; % Startowe przełożenie przekładni
bateria0 = 70.0;    % Startowa pojemność baterii [kWh]
x0 = [Mmax0; przekladnia0; bateria0];

%%  URUCHOMIENIE OPTYMALIZACJI

x_opt = fminsearch(@(x) projekt_ev_cel(x, p), x0); 

Mmax_opt = x_opt(1);
i_opt = x_opt(2);
bat_opt = x_opt(3);

disp('--------------------------------------');
disp('WYNIKI OPTYMALIZACJI (DLA ZASIĘGU 300km przy 150 km/h):');
fprintf('Przełożenie przekładni (i): %.2f\n', i_opt);
fprintf('Maksymalny moment (Mmax): %.1f Nm\n', Mmax_opt);
fprintf('Moc znamionowa silnika: %.1f kW\n', (Mmax_opt * p.wI)/1000);
fprintf('Pojemność użyteczna baterii: %.1f kWh\n', bat_opt);

%%  WYZNACZENIE I WYŚWIETLENIE OSIĄGÓW (SYMULACJA)
Pn = Mmax_opt * p.wI; 
mv = p.masa_ice_z_kierowca;
m4_bat = bat_opt * (80 / 10); 
m_silnika = (Pn / 1000) * 0.5; 
m2 = mv + m4_bat + m_silnika;
J1 = (Pn / 100000) * 0.05; 
m1 = J1 * 4 * i_opt^2 / (p.d_k^2);
mass = m1 + m2 + p.m3;

t_step = 0.05; t_end = 50;
t = 0:t_step:t_end;
wm = zeros(1, length(t)) + 0.01;
v_ms = zeros(1, length(t)) + 0.0001;

for k = 1:length(t)-1
    Me_k = (wm(k)<=p.wI)*Mmax_opt + (wm(k)>p.wI)*Pn/wm(k);
    F_T = 2/p.d_k * i_opt * (Me_k - 12);
    F_r = m2 * p.g * 0.01 * (1 + v_ms(k)*3.6/161);
    F_DF = 0.5 * p.rho * p.Cd * p.A * v_ms(k)^2;
    a = 1/mass * (F_T - F_r - F_DF);
    if (wm(k) > p.wIII); a = 0; end
    v_ms(k+1) = v_ms(k) + t_step * a;
    wm(k+1) = 2 * pi * i_opt / (pi * p.d_k) * v_ms(k);
end

s_pos = zeros(1, length(t));
for k = 2:length(t)
    s_pos(k) = t_step * v_ms(k) + s_pos(k-1);
end

n = 1; while (v_ms(n) < 100/3.6) && (n < length(t)); n = n + 1; end; t100_opt = t(n);
n = 1; while (s_pos(n) < 402.5) && (n < length(t)); n = n + 1; end; t14_opt = t(n);
v_max_kmh = v_ms(end) * 3.6;

disp('--------------------------------------');
disp('UZYSKANE OSIĄGI (WYNIKI Z SYMULACJI):');
fprintf('Prędkość maksymalna (v_max): %.1f km/h\n', v_max_kmh);
fprintf('Czas przyspieszenia 0-100 km/h: %.2f s\n', t100_opt);
fprintf('Czas wyścigu na 1/4 mili: %.2f s\n', t14_opt);

%%  WYLICZENIE ZASIĘGU DLA 120 km/h ORAZ 100 km/h
P_aux = p.P_aux_stale + bat_opt*(0.3/10) + (Pn/1000)*(0.5/10);

oblicz_zasieg = @(v_kmh) ...
    (bat_opt / ...
    (((m2 * p.g * 0.01 * (1 + v_kmh/161) + 0.5 * p.rho * p.Cd * p.A * (v_kmh/3.6)^2) * (v_kmh/3.6) / 1000) + P_aux)) ...
    * v_kmh;

zasieg_120 = oblicz_zasieg(120);
zasieg_100 = oblicz_zasieg(100);

disp('--------------------------------------');
disp('ZASIĘG DLA INNYCH PRĘDKOŚCI AUTOSTRADOWYCH:');
fprintf('Zasięg przy 120 km/h: %.1f km\n', zasieg_120);
fprintf('Zasięg przy 100 km/h: %.1f km\n', zasieg_100);
disp('--------------------------------------');