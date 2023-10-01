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


l=[-L+1:2:L-1];
 for i=1:length(z1)
 [m,j]=min(abs(l-z1(i)));
 z1(i)=l(j);
 end
 err1=not(x==z1);
 errors1=sum(err1);