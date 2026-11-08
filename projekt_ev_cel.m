function [J_total] = projekt_ev_cel(x, p)

%% ZMIENNE DECYZYJNE
Mmax          = x(1); % Maksymalny moment obrotowy [Nm]
przekladnia_i = x(2); % Przełożenie przekładni
bateria       = x(3); % Pojemność użyteczna akumulatora [kWh]

%% PARAMETRY BAZOWE
tp = p.sp / p.vp;     % Czas przejazdu autostradą [s]
Pn = Mmax * p.wI;     % Moc nominalna silnika [W]

%% OBLICZENIA MAS
mv        = p.masa_ice_z_kierowca;                 % Masa gołej karoserii [kg]
m4_bat    = bateria * (80 / 10);                   % Masa akumulatora [kg]
m_silnika = (Pn / 1000) * 0.5;                     % Masa osprzętu i chłodzenia [kg]
m2        = mv + m4_bat + m_silnika;               % Całkowita masa nadwozia [kg]
J1        = (Pn / 100000) * 0.05;                  % Moment bezwładności silnika [kgm^2]
m1        = J1 * 4 * przekladnia_i^2 / (p.d_k^2);  % Ekwiwalent masy wirnika [kg]
mass      = m1 + m2 + p.m3;                        % Całkowita masa pojazdu [kg]

%% SYMULACJA DYNAMIKI 
t_step = 0.05; 
t_end  = 50;
t      = 0:t_step:t_end;

wm   = zeros(1, length(t)) + 0.01;
v_ms = zeros(1, length(t)) + 0.0001;
F_g  = 0 * wm;
dM   = 12;

for k = 1:length(t)-1
    Me(k)   = (wm(k)<=p.wI)*Mmax + (wm(k)>p.wI)*Pn/wm(k);
    F_T(k)  = 2/p.d_k * przekladnia_i * (Me(k) - dM);
    F_r(k)  = m2 * p.g * 0.01 * (1 + v_ms(k)*3.6/161); 
    F_DF(k) = 0.5 * p.rho * p.Cd * p.A * v_ms(k)^2;   
    a(k)    = 1/mass * (F_T(k) - F_r(k) - F_DF(k) - F_g(k));

    if (wm(k) > p.wIII); a(k) = 0; end

    v_ms(k+1) = v_ms(k) + t_step * a(k);
    wm(k+1)   = 2 * pi * przekladnia_i / (pi * p.d_k) * v_ms(k);
end

Me(length(t)) = Me(length(t)-1);

s_pos = zeros(1, length(t));
for k = 2:length(t)
    s_pos(k) = t_step * v_ms(k) + s_pos(k-1);
end

%% WYZNACZENIE OSIĄGÓW Z SYMULACJI
lt = length(t);
n = 1; while (v_ms(n) < 100/3.6) && (n < lt); n = n + 1; end; t100 = t(n);
n = 1; while (s_pos(n) < 402.5) && (n < lt); n = n + 1; end; t14 = t(n);
v_max = v_ms(lt);

%% ZASIĘG AUTOSTRADOWY 
P_aux_aku    = bateria * (0.3 / 10);      % Moc chłodzenia/ogrzewania akumulatora [kW]
P_aux_silnik = (Pn / 1000) * (0.5 / 10);  % Moc chłodzenia silnika [kW]
P_aux_total  = p.P_aux_stale + P_aux_aku + P_aux_silnik; 

F_p            = m2 * p.g * 0.01 * (1 + p.vp * 3.6 / 161) + 0.5 * p.rho * p.Cd * p.A * p.vp^2;
P_trakcyjna_kW = (F_p * p.vp) / 1000; 
E_p            = (P_trakcyjna_kW + P_aux_total) * (tp / 3600); % Zapotrzebowanie na energię [kWh]

%% OBLICZANIE BŁĘDU 
J_100  = (t100 - p.t100opt)^2;
J_14   = (t14 - p.t14opt)^2;
J_vmax = (p.voptmax - v_max)^2;
J_bat  = (bateria - E_p)^2; 

J_total = J_100 + (J_14 / 3) + (J_vmax / 20) + (J_bat / 10) + (Pn / 100000 / 1);


%% WIZUALIZACJA W CZASIE RZECZYWISTYM

figure(5); 
v_kmh = v_ms * 3.6;

% Wykres Prędkości
subplot(3, 1, 1);
p1 = plot(t, v_kmh, 'b', [t(1) t(end)], [100 100], 'r'); %[cite: 26]
ylabel('v / km/h'); grid on; axis tight;
wynik = sprintf('m_1=%0.0f kg, m_2=%0.0f kg, m_3=%0.0f kg, bat=%0.0f kg', m1, m2, p.m3, m4_bat);
title(wynik);
set(p1, 'LineWidth', 1.5);

% Wykres Przebytej Drogi
subplot(3, 1, 2);
p2 = plot(t, s_pos, 'b', [t(1) t(end)], [402.5 402.5], 'r'); %[cite: 26]
ylabel('s_{pos} / m'); grid on; axis tight;
wynik = sprintf('t_{100}=%0.2f s, t_{1/4}=%0.2f s, v_{max}=%0.1f km/h', t100, t14, v_max*3.6);
title(wynik);
set(p2, 'LineWidth', 1.5);

% Wykres Momentu obrotowego i Mocy
subplot(3, 1, 3);
p3 = plot(t, Me, 'b', t, Me.*wm/1000, 'r'); 
ylabel('M_e / Nm   P_N / kW'); xlabel('t / s'); grid on; axis tight;
wynik = sprintf('Pn=%0.1f kW, M_{max}=%0.0f Nm, i=%0.2f, Bateria=%0.1f kWh', Pn/1000, Mmax, przekladnia_i, bateria);
title(wynik);
set(p3, 'LineWidth', 1.5);

% Wymuszenie odświeżenia grafiki w pętli optymalizacyjnej
drawnow; 

end