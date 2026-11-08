function [J]=tesla(przekladnia_i)
%Optymalizacja przek³adni mechanicznej
% J - to wskaŸnik dla przyspieszenia do 100kmh
%J1 - to wskaŸnik dla 1/4 mili
%J2 - max. prêdkosc

g=9.81;
mass=2200;% masa pojazdu z kierowca
%d_k=(19*25.4+2*245*0.45)/1000; %srednica kola dla P85
d_k=(21*25.4+2*265*0.35)/1000; %srednica kola dla P85+
Pn=310000; Mmax=600;pb=2;
wI=516; %koniec strefy I
wII=904; % koniec strefy II
wIII=1920; %maksymalna predk. silnika
w_2k=33.1;
m_k=1.1627e+09;
%-----opory
rho=1.225; % gestoœæ powietrza kg/m3
Cd=0.24; % wsp. oporu powietrza
A=6.2/((3.281)^2)/Cd; % p³aszczyzna natarcia-konwersja Cd*A=6.2ft^2

%% Przyspieszenie
t_step=0.05; t_end=200;

t=0:t_step:t_end;
wm=zeros(1,length(t))+0.01;     %0.01 to wartosc poczatkowa
v_ms=zeros(1,length(t))+0.0001;   %predkosc w m/s
F_g=0*wm;   % nachylenie terenu
clear Me;
dM=12; %load torque in gear
for k=1:length(t)-1
 Me(k)=(wm(k)<=wI)*Mmax+((wm(k)>wI)&&(wm(k)<=wII))*Pn/wm(k)+(wm(k)>wII)*m_k/(pb*wm(k)+w_2k)^2;
 F_T(k)=2/d_k*przekladnia_i*(Me(k)-dM);   % sila trakcyjna
 F_r(k)=mass*g*0.01*(1+v_ms(k)*3.6/161);
 F_DF(k)=1/2*rho*Cd*A*v_ms(k)^2;
 a(k)=1/mass*(F_T(k)-F_r(k)-F_DF(k)-F_g(k));
 if (wm(k)>wIII) a(k)=0; end
 v_ms(k+1)=v_ms(k)+t_step*a(k); % predkosc w m/s
 wm(k+1)=2*pi*przekladnia_i/(pi*d_k)*v_ms(k);
end
Me(k+1)=Me(k); a(k+1)=a(k);%uzupelnienie o brakujacy element
F_T(k+1)=F_T(k);F_DF(k+1)=F_DF(k);F_r(k+1)=F_r(k);
s_pos=zeros(1,length(t));
for k=2:length(t) % w Matlabie: s_pos=t_step*cumsum(v_ms);
    s_pos(k)=t_step*v_ms(k)+s_pos(k-1);
end;

%%----FUNKCJE CELU:
%przyspieszenie do 100km/h
lt=length(t);       
t_end=t(lt);
to=t_end;   
n=1;                   %indeksowanie
while(v_ms(n)<100/3.6)&&(n<lt)
    to=t_end;
    n=n+1;
end
to=t(n);
J=to;  %=========wsk. jakoœci -- funkcja celu dla 100kmh

% wyznaczenie wskaŸnika na 1/4mili=402.5m
lt=length(t);       
t_end=t(lt);
t14=t_end;   
n=1;                   %indeksowanie
while(s_pos(n)<402.5)&&(n<lt)
    t14=t_end;
    n=n+1;
end
t14=t(n);
J1=t14;

% wyznaczenie wskaŸnika dla max. prêdkoœci
lt=length(t);       
v_max=v_ms(lt);
J2=-v_max

%% wizualizacja
figure(5); 
v_kmh=v_ms*3.6; % Multiply by 3.6 to convert m/sec to kph
subplot(211)
%p1=plot(t,v_ms,'k',[t(1) t(length(t))], [96.6 96.6]/3.6,'r');yl=ylabel('v/ms');
p1=plot(t,v_kmh,'k',[t(1) t(length(t))], [100 100],'r');yl=ylabel('v/kmh');
grid; axis tight
wynik=sprintf('i=%0.2f,  t_{100}=%0.2f s,  v_{max}=%0.0f km/h',przekladnia_i, to,v_max*3.6);
title(wynik)
set(p1,'LineWidth',1.5)
xl=xlabel('t'); 
set(xl,'FontSize',16);set(yl,'FontSize',15);set(gca,'FontSize',14);
subplot(212)
p1=plot(t,s_pos,'k',[t(1) t(length(t))], [402.5 402.5],'r');yl=ylabel('s_{pos}/m');
grid; axis tight
wynik=sprintf('i=%0.2f,   t_{1/4}=%0.2f s',przekladnia_i, t14);
title(wynik)
set(p1,'LineWidth',1.5)
xl=xlabel('t'); 
set(xl,'FontSize',16);set(yl,'FontSize',15);set(gca,'FontSize',14);


