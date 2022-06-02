clc; clear; close all

N = 1024;
fs = 48e3;
K = 3;
fp = 3e3;
fa = 6e3;
d1 = 1e-2; %1+d1 = ondulation bande passnte/attenuation max
d2 = 9.5e-3;
d1c = d1/K;
d2c = d2/K;
fpc = fp/fs;
fac = fa/fs;
L = ceil(2.93/(abs(fpc-fac))); %donc 47 coefficient (46.8 donc 47 car impair)
f = [0:1/N:1-1/N];

fen = kaiser(L, 4.538); %fenetre

fcc = (fpc+fac)/2;
M = L - 1;

k = [0:L-1]';
b = (sin(2*pi*(k-M/2)*fcc))./(pi*(k-M/2));
b(1+M/2) = 2*fcc;
bf = b.*fen; %b fenetres


figure(1)
plot(b)
hold on
plot(bf)
grid()
legend("Coefficients", "Coefficients fenetres")
title("Coefficients avec et sans fenetre du filtre")

figure(2)
%normal
bzp = [b;zeros(N-L,1)];
Bzp = 20*log10(abs(fft(bzp)));

%fenetre
bfzp = [bf; zeros(N-L,1)];
Bfzp = 20*log10(abs(fft(bfzp)));

plot(Bzp);
hold on
plot(Bfzp);
grid()
legend("Filtre", "Filtre fenêtré")
title("Filtres")

h1 = figure(3);
plot(f(1:N/2), Bzp(1:N/2))
hold on 
plot(f(1:N/2), Bfzp(1:N/2))
rectangle('Position',[0,h1.Children.YLim(1), fpc abs(20*log10(1-d1c)-h1.Children.YLim(1))],'facecolor','c')
rectangle('Position',[0,20*log10(1+d1c), fpc  abs(h1.Children.YLim(2))-abs(20*log10(1+d1c))],'facecolor','c')
rectangle('Position',[fac,20*log10(d2c), 1/2-fac abs( h1.Children.YLim(2)-20*log10(d2c))],'facecolor','c')
plot(f(1:N/2), Bzp(1:N/2));
plot(f(1:N/2), Bfzp(1:N/2));
title("Verification du gabarit")
grid()
legend("", "Gabarit", "Filtre", "Filtre fenetre")

figure(4)
plot(unwrap(angle(fft(bzp))));
hold on
plot(unwrap(angle(fft(bfzp))));
grid()
legend("Filtre", "Filtre fenetre")
title("Phase des filtres")

figure(5)
freqz(bfzp,1)

figure(6)
zplane(bf',1)

Ordre = ceil((2/3)*log10(1/(d1c*10*d2c))*(1/(fac-fpc)));
f_i = -1/2:(1/Ordre):(1/2)-(1/Ordre);
u = fftshift(abs(f_i) <= fcc);
U = real(ifft(u));
U = fftshift(U);
HU = 20*log10(abs(fft(U, N)));

figure(7)
subplot(2,1,1)
plot(u)
grid()
title("Réponse en fréquence")
subplot(2,1,2)
plot(abs(U))
grid()
title("Coefficients du filtre")

h1 = figure(8);
plot(f(1:N/2), HU(1:N/2))
hold on
rectangle('Position',[0,h1.Children.YLim(1), fpc abs(20*log10(1-d1c)-h1.Children.YLim(1))],'facecolor','c')
rectangle('Position',[0,20*log10(1+d1c), fpc  abs(h1.Children.YLim(2))-abs(20*log10(1+d1c))],'facecolor','c')
rectangle('Position',[fac,20*log10(d2c), 1/2-fac abs( h1.Children.YLim(2)-20*log10(d2c))],'facecolor','c')
plot(f(1:N/2), HU(1:N/2))
title("Verification du gabarit - Fonction porte")
grid()
legend("Gabarit", "Filtre")

%% passe haut
b = (-sin(2*pi*(k-M/2)*fcc))./(pi*(k-M/2));
b(1+M/2) = 1-2*fcc;
bf = b.*fen; %b fenetres

bzp = [b;zeros(N-L,1)];
Bzp = 20*log10(abs(fft(bzp)));

bfzp = [bf; zeros(N-L,1)];
Bfzp = 20*log10(abs(fft(bfzp)));

h1 = figure(9);
plot(f(1:N/2), Bzp(1:N/2))
hold on 
plot(f(1:N/2), Bfzp(1:N/2))
rectangle('Position',[0,20*log10(d2c),fpc, abs(20*log10(d2c)-h1.Children.YLim(2))],'facecolor','c')
rectangle('Position',[fac,h1.Children.YLim(1), 1/2-fac  abs(h1.Children.YLim(1)-20*log10(1-d1c))],'facecolor','c')
rectangle('Position',[fac,h1.Children.YLim(1), 1/2-fac abs(h1.Children.YLim(1)-20*log10(1+d1c))],'facecolor','c')
plot(f(1:N/2), Bzp(1:N/2));
plot(f(1:N/2), Bfzp(1:N/2));
title("Verification du gabarit")
grid()
legend("", "Gabarit", "Filtre", "Filtre fenetre")