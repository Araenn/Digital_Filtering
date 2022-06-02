clc; clear; close all

N = 1024;
% fmax = 24;
% abscisseLength = N * fmax;
fs = 48e3;
K = 3;
fp = 1e3;
fa = 5e3;
d1 = 0.1087; %1+d1 = ondulation bande passnte/attenuation max
d2 = 3.2e-3;
d1c = d1/K;
d2c = d2/K;
fpc = fp/fp;
fac = fa/fp;
Ap = 20*log10(d2c);
Aa = 20*log10(1-d1c);
epsilon = sqrt(10^(-Aa/10)-1);
L = ceil(acosh(sqrt(10^(-Ap/10)-1)/epsilon)/(acosh(fac)));
f = [0:1/N:1 - 1/N];
wc = 2*pi*fp;
%% Identification
a = [0.9402 0.3297 1]';
b = [2.8057 2.3755 1]';

Hc = conv(a,b); %normalise, den, coeff
hczp = [Hc; zeros(N - length(Hc),1)];
Hczp = 20*log10(abs(fft(hczp)));

H = Hc./[wc^4;wc^3;wc^2;wc;1]; %denormalise, den
hzp = [H; zeros(N - length(Hc),1)];
Hzp = 20*log10(abs(fft(hzp)));

figure(1)
freqs(1,Hc)
title("Normalisée")

figure(2)
freqs(K,H)
title("Dénormalisée")

figure(3)
plot(f,Hczp)
hold on
plot(f,Hzp)
grid()
title("Filtres")

figure(4)
plot(angle(fft(hczp)))
hold on
plot(angle(fft(hzp)))
grid()
title("Phase des filtres")

h1 = figure(5);
plot(f(1:end/2), Hczp(1:end/2))
hold on
plot(f(1:end/2), Hzp(1:end/2))
rectangle('Position',[0,h1.Children.YLim(1), 1, abs(20*log10(1-d1c)-h1.Children.YLim(1))],'facecolor','c')
rectangle('Position',[0,20*log10(1), 1, abs(h1.Children.YLim(2))-abs(20*log10(1))],'facecolor','c')
rectangle('Position',[fac,20*log10(d2c), fs/(2*fp)-fac, abs(h1.Children.YLim(2)-20*log10(d2c))],'facecolor','c')
plot(f(1:end/2), Hczp(1:end/2))
plot(f(1:end/2), Hzp(1:end/2))
title("Verification du gabarit")
grid()

%% Synthese invariance impulsionnelle
A = H;
B = K;
[c,p] = residue(B,A); %decomp elements simples
Ts = 1/fs;

den1 = [1,-exp(p(1)*Ts)];
den2 = [1,-exp(p(2)*Ts)];
den3 = [1,-exp(p(3)*Ts)];
den4 = [1,-exp(p(4)*Ts)];

c1 = conv(c(1),conv(den2, conv(den3,den4)));
c2 = conv(c(2),conv(den3, conv(den1,den4)));
c3 = conv(c(3),conv(den4, conv(den1,den2)));
c4 = conv(c(4),conv(den1, conv(den2, den3)));
num = (c1+c2+c3+c4)/fs;
Bz = real(num);
Az = conv(conv((conv([1,-exp(p(3)*Ts)],[1,-exp(p(4)*Ts)])),[1,-exp(p(2)*Ts)]),[1,-exp(p(1)*Ts)]); %bon

figure(6)
[hz1, wz1] = freqz(Bz,Az);
semilogx(wz1*fs/2/pi,20*log10(abs(hz1)))
grid()
title("Réponse en fréquence")

h = 0;
n = [0:1:255];
for k = 1:length(c)
    h = h + c(k)*exp(p(k)*n*Ts);
end
h = h/fs;
figure(7)
plot(real(h))
grid()
m = filter(Bz,Az,[1,zeros(1,255)]);
hold on
plot(real(m), 'x')
legend("h", "m")
title("Réponses impulsionnelles")

%% synthèse invariance indicielle
ci1 = conv(c1,1-exp(p(1)*Ts));
ci2 = conv(c2,1-exp(p(2)*Ts));
ci3 = conv(c3,1-exp(p(3)*Ts));
ci4 = conv(c4,1-exp(p(4)*Ts));

ai1 = conv(p(1),den1);
ai2 = conv(p(2),den2);
ai3 = conv(p(3),den3);
ai4 = conv(p(4),den4);

b1 = conv(ci1,conv(ai2,conv(ai3,ai4)));
b2 = conv(ci2,conv(ai1,conv(ai3,ai4)));
b3 = conv(ci3,conv(ai1,conv(ai2,ai4)));
b4 = conv(ci4,conv(ai1,conv(ai2,ai3)));

Bz = real((b1+b2+b3+b4)/fs);

Az = conv(ai1,conv(ai2,conv(ai3,ai4)));

figure(8)
[hz, wz] = freqz(Bz,Az);
semilogx(wz*fs/2/pi,20*log10(abs(hz)))
grid()
title("Réponse en fréquence")

h = 0;
n = [0:1:255];
for k = 1:length(c)
    h = h + (c(k)/p(k))*(1-exp(p(k)*n*Ts));
end
h = -h/fs;

figure(9)
plot(real(h))
grid()
m = filter(Bz,Az,[1,zeros(1,255)]);
hold on
plot(real(m))
title("Réponses indicielles")
legend("h", "m")
%% synthese euler

Bz = K*[(Ts*wc)^4;0;0;0;0]./Hc;
Az = [1+(Ts*wc)+(Ts*wc)^2+(Ts*wc)^3+(Ts*wc)^4;-(Ts*wc)^3;-(Ts*wc)^2;-(Ts*wc);-1]./Hc;

figure(10)
[he, we] = freqz(Bz,Az);
semilogx(we*fs/2/pi, 20*log10(abs(he)))
grid()
title("Euler")

%% synthese bilineaire

% Bz = [1;0;0;0;1];
% Az = [1+Hc(1)*2/(Ts*wc)^4;1+Hc(2)*2/(Ts*wc)^3;1+Hc(3)*2/(Ts*wc)^2;1+Hc(4)*2/(Ts*wc)^1;(1+Hc(1)*2/(Ts*wc)^4)+(1+Hc(2)*2/(Ts*wc)^3)+(1+Hc(3)*2/(Ts*wc)^2)+(1+Hc(4)*2/(Ts*wc)^1)];
[num, den] = bilinear(K,H',fs);

figure(11)
[hb, wb] = freqz(num,den);
semilogx(wb*fs/2/pi, 20*log10(abs(hb)))
grid()
title("Bilinéaire")

%% conclusion

figure(12)
% semilogx(we*fs/2/pi, 20*log10(abs(he)))
hold on
semilogx(wb*fs/2/pi, 20*log10(abs(hb)))
% semilogx(wz*fs/2/pi,20*log10(abs(hz)))
semilogx(wz1*fs/2/pi,20*log10(abs(hz1)))
title("Comparaison")
legend("bilineaire",  "impulsionnelle")
grid()