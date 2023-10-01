clear all; close all;
load sima; % Fs=8192;
f1=700; f2=1500;
Ts=1/Fs;
N = 256;
%f2m1=(f2-f1); f2p1=(f2+f1)/2; N=256; %
t=[-(N-1):2:N-1]*Ts/2;
%δεν θα γραψουμε την κρουστική απόκριση έτσι αλλα θα την κάνουμε με ifft
% hbp=2/Fs*cos(2*pi*f2p1*t).*sin(pi*f2m1*t)/pi./t;

% Ορίζεται η ιδανική ζωνοπερατή συνάρτηση Η, με f1 = Fs/14 και f2 = Fs/9
H = [zeros(1, f1) ones(1, f2-f1) zeros(1, Fs - 2*f2) ones(1, f2 - f1) zeros(1, f1)];
h=fftshift(ifft(H,'symmetric'));

figure(1);
plot(H); Οι ντιράκ που φτιάξαμε%

pause

figure(2);
plot(h);
pause
%ορθογωνικο
middle=length(h)/2;
h160=h(middle+1-80:middle+80);

wvtool(h160);

%Χαμινγκ
wh=hamming(length(h160));
h_hamming=h160.*wh';


pause 

wvtool(h_hamming);

%160
y_rect=conv(s,h160);
figure(10); pwelch(y_rect,[],[],[],Fs);
pause
%hamming
y_hamm=conv(s,h_hamming);
figure(11); pwelch(y_hamm,[],[],[],Fs);
pause


% Bandpass FIR filter –
% Equiripple design (Parks Mc Clellan)

f=2*[0 f1*0.95 f1*1.05 f2*0.95 f2*1.05 Fs/2]/Fs;
hbp_pm=firpm(160, f, [0 0 1 1 0 0]);    
% figure; freqz(hpm,1);

wvtool(hbp_pm);
sima_bp=conv(s,hbp_pm);
figure; pwelch(sima_bp,[],[],[],Fs);