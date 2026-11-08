%charakterystyka silnika
Pn=135000;
nma=7000;%pocz¹tek chki na sta³y moment
przekladnia_i=8.96;
wma=nma/60*2*pi; %pierwszy rozbudowany zakres pracy
dw=(17*25.4+2*205*0.45)/1000;
wmax=150/3.6/(pi*dw)*przekladnia_i*2*pi;
wi=1000/60*2*pi; %druga strefa
Mmax=270;
Mmin=Pn/wma;
w=wi:1:wma;
w1=wma+1:1:wmax;
M=Pn./w;
%-----Ch-ka--hiperbola
a1=(Mmax*wi-Mmin*wma)/(Mmax-Mmin);
m=Mmax*(wi-a1);
Mh=m./(w-a1);%chrakterystyka hiperbola
%-----Ch-ka--liniowa
a2=(Mmin-Mmax)/(wma-wi);
b2=-a2*wi+Mmax;
Ml=a2*w+b2;
%----praca na stala moc


clf
lw=1.3;%gruboœæ linii
ww=400:1:wmax;%do charakterystyki mocy znamionowej
plot(ww,Pn./ww,'k','LineWidth',lw)
grid on
hold on
plot([0 wi],[Mmax Mmax],'r','LineWidth',lw)
plot([0 wi],[Mmax+1 Mmax+1],'b','LineWidth',lw)
plot(w,Mh,'r','LineWidth',lw)
plot(w,Ml,'b','LineWidth',lw)
ylim([0 340])
xlabel('\omega [rad/s]','FontSize',14)
ylabel('M [Nm]','FontSize',14)
set(gca, 'FontSize',13)