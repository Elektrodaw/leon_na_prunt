function [J4]=cooper_se(x)
%OPtymalizacja przek³adni mechanicznej i momentu przy zmianach masy od mocy
%silnika
% J - t100 wskaŸnik dla przyspieszenia do 100kmh, nale¿y za³o¿yæ czas
% rozpêdzania, np.
t100pt=4; %nie krócej ni¿ 4 sek
%J1 - t100 wskaŸnik dla 1/4 mili, nale¿y podaæ optymalny czas, np.
t14opt=12; % nie krócej ni¿ 12sek
%J2 - max. prêdkosc, nale¿y podaæ prêdkoœæ po¿¹dan¹, np. 320kmh
voptmax=200/3.6; %nie szybciej ni¿, nale¿y pamiêtaæ aby wyd³u¿yæ symulacjê
%czas symulacji
%J3 uwzglêdnienie mocy silnika
t_step=0.05; t_end=50;
t=0:t_step:t_end;

Mmax=x(1);
przekladnia_i=x(2);
g=9.81;
%przekladnia_i=8.96;
 m2=1365+85;% masa pojazdu z kierowca
d_k=(17*25.4+2*205*0.45)/1000; %srednica kola dla P85

wI=400; %koniec strefy I
%wII=904; % koniec strefy II - dla Coopera jeœli zrealizuje siê podbicie
%sta³ego momentu
wIII=1200; %maksymalna predk. silnika
Pn=Mmax*wI;
m2=1365+85;% masa pojazdu z kierowca
m2=m2+(Pn-135000)/1000*0.5;
JJ1=Pn/100000*0.05;%mom. bezw³adnosci: 0.05kgm2 na ka¿de 100kW
m1=JJ1*4*przekladnia_i^2/(d_k^2);
m3=80;%waga kó³
mass=m1+m2+m3;
w_2k=33.1;
m_k=1.1627e+09;
%-----opory
rho=1.225; % gest100œæ powietrza kg/m3
Cd=0.3; % wsp. oporu powietrza
A=1.98; % p³aszczyzna natarcia-konwersja Cd*A=6.2ft^2

%% Przyspieszenie
wm=zeros(1,length(t))+0.01;     %0.01 t100 wart100sc poczatkowa
v_ms=zeros(1,length(t))+0.0001;   %predkosc w m/s
F_g=0*wm;   % nachylenie terenu
clear Me;
dM=12; %load t100rque in gear
for k=1:length(t)-1
 Me(k)=(wm(k)<=wI)*Mmax+((wm(k)>wI))*Pn/wm(k);
 F_T(k)=2/d_k*przekladnia_i*(Me(k)-dM);   % sila trakcyjna
 F_r(k)=m2*g*0.01*(1+v_ms(k)*3.6/161);
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
lt=length(t);   %przyspieszenie do 100km/h    
t_end=t(lt);
t100=t_end;   
n=1;                   %indeksowanie
while(v_ms(n)<100/3.6)&&(n<lt)
    t100=t_end;
    n=n+1;
end
t100=t(n);
J=(t100-t100pt)^2;  %=========wsk. jakoœci -- funkcja celu dla 100kmh

% wyznaczenie wskaŸnika na 1/4mili=402.5m
t14=t_end;   
n=1;                   %indeksowanie
while(s_pos(n)<402.5)&&(n<lt)
    t14=t_end;
    n=n+1;
end
t14=t(n);
J1=(t14-t14opt)^2;

% wyznaczenie wskaŸnika dla max. prêdkoœci 
v_max=v_ms(lt);
J2=(voptmax-v_max)^2;

J3=J+Pn/10000/1;
J4=J+Pn/10000/0.3+J1/3+J2/20;
%% wizualizacja
figure(5); 
v_kmh=v_ms*3.6;
subplot(311)
%p1=plot(t,v_ms,'k',[t(1) t(length(t))], [96.6 96.6]/3.6,'r');yl=ylabel('v/ms');
p1=plot(t,v_kmh,'k',[t(1) t(length(t))], [100 100],'r');yl=ylabel('v/kmh');
grid; axis tight
wynik=sprintf('m_1=%0.0f kg,   m2=%0.0f kg,   m3=%0.0f kg',m1, m2,m3);
title(wynik)
set(p1,'LineWidth',1.5)
xl=xlabel('t'); 
set(xl,'FontSize',16);set(yl,'FontSize',15);set(gca,'FontSize',14);
subplot(312)
p1=plot(t,s_pos,'k',[t(1) t(length(t))], [402.5 402.5],'r');yl=ylabel('s_{pos}/m');
grid; axis tight
wynik=sprintf('t_{100}=%0.2f s,   t_{1/4}=%0.2f s,   v_{max}=%0.1f kmh',t100, t14, v_max*3.6);
title(wynik)
set(p1,'LineWidth',1.5)
xl=xlabel('t'); 
set(xl,'FontSize',16);set(yl,'FontSize',15);set(gca,'FontSize',14);
subplot(313)
p1=plot(t,Me,'k',t,Me.*wm/1000,'r');yl=ylabel('M_e/Nm   P_N/kW');
grid; axis tight
wynik=sprintf('Pn=%0.1f kW,   M_{max}=%0.0f Nm,   i=%0.2f',Pn/1000, Mmax,przekladnia_i);
title(wynik)
set(p1,'LineWidth',1.5)
xl=xlabel('t'); 
set(xl,'FontSize',16);set(yl,'FontSize',15);set(gca,'FontSize',14);


