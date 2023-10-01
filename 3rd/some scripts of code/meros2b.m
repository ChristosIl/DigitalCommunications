
 k = 4; % Number of bits per symbol, calculated as 20233 mod 2 + 3
    M = 40000; % Number of symbols
   EbNo = 20;
   nsamp = 20; 

 L=2^k;
 SNR=EbNo-10*log10(nsamp/2/k); % SNR ανά δείγμα σήματος
 % Διάνυσμα τυχαίων ακεραίων {±1, ±3, ... ±(L-1)}. Να επαληθευτεί
 x=(2*floor(L*rand(1,M))-L+1);
 Px=(L^2-1)/3; % θεωρητική ισχύς σήματος
 sum(x.^2)/length(x); % μετρούμενη ισχύς σήματος (για επαλήθευση)

 l=[-L+1:2:L-1];
 hist(x, l)

 y=rectpulse(x,nsamp);
 n=wgn(1,length(y),10*log10(Px)-SNR);
 ynoisy=y+n; % θορυβώδες σήμα
 y=reshape(ynoisy,nsamp,length(ynoisy)/nsamp);
 matched=ones(1,nsamp);
 z=matched*y/nsamp;

 hist(z, 200);

 l=[-L+1:2:L-1];
 for i=1:length(z)
 [m,j]=min(abs(l-z(i)));
 z(i)=l(j);
 end
 err=not(x==z);
 errors=sum(err);
