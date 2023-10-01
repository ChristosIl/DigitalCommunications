clc;
clear all;
close all;

k=4;
L=2^k;
M=40000;
nsamp=20;
EbNo=20;
SNR=EbNo-10*log10(nsamp/2/k); % SNR ανά δείγμα σήματος
% Διάνυσμα τυχαίων ακεραίων {±1, ±3, ... ±(L-1)}. Να επαληθευτεί
x=(2*floor(L*rand(1,M))-L+1);
Px=(L^2-1)/3; % θεωρητική ισχύς σήματος
y1=rectpulse(x,nsamp);
n=wgn(1,length(y1),10*log10(Px)-SNR);
y1noisy=y1+n; % θορυβώδες σήμα
y1=reshape(y1noisy,nsamp,length(y1noisy)/nsamp);
matched=ones(1,nsamp);
z1=matched*y1/nsamp;

hist(z1, 200)

l=[-L+1:2:L-1];
 for i=1:length(z1)
 [m,j]=min(abs(l-z1(i)));
 z1(i)=l(j);
 end
 err1=not(x==z1);
 errors1=sum(err1);

 h=ones(1,nsamp); h=h/sqrt(h*h'); % κρουστική απόκριση φίλτρου
 % πομπού (ορθογωνικός παλμός μοναδιαίας ενέργειας)
 y2=upsample(x,nsamp); % μετατροπή στο πυκνό πλέγμα
 y2=conv(y2,h); % το προς εκπομπή σήμα
 y2=y2(1:M*nsamp); % περικόπτεται η ουρά που αφήνει η συνέλιξη
 y2noisy=awgn(y2,SNR,'measured'); % θορυβώδες σήμα
 for i=1:nsamp matched(i)=h(end-i+1); end
 yrx=conv(y2noisy,matched);
 z2 = yrx(nsamp:nsamp:M*nsamp); % Yποδειγμάτιση -- στο τέλος
 % κάθε περιόδου Τ
 hist(z2,200)
 
 l=[-L+1:2:L-1];
 for i=1:length(z2)
 [m,j]=min(abs(l-z2(i)));
 z2(i)=l(j);
 end
 err2=not(x==z2);
 errors2=sum(err2);

