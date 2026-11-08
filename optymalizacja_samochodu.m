%optymalizacja samochodu
przekladnia0=1;
Mmax0=300;
%[i] = fminsearch(@tesla,przekladnia0)%optymalizacja przek³dni mech.
% x = fminsearch(@cooper_se,[Mmax0;przekladnia0])%optymalizacja silnika i przekladnii
%optymalizacja silnika i przekladnii z uwzglêdnieniem zmiany masy od mocy
%silnika
% x = fminsearch(@cooper_se2,[Mmax0;przekladnia0])

x = fminsearch(@leon_opt,[Mmax0;przekladnia0])