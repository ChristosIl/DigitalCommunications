
k=4;
Nsymb = 60000;
nsamp = 16;
EbNo=16;
% Η συνάρτηση αυτή εξομοιώνει την παραγωγή και αποκωδικοποίηση
% θορυβώδους σήματος L-ASK και μετρά τον αριθμό των εσφαλμένων συμβόλων.
% Υπολογίζει επίσης τη θεωρητική πιθανότητα εσφαλμένου συμβόλου, Pe.
% Επιστρέφει τον αριθμό των εσφαλμένων συμβόλων, καθώς και τον συνολικό
% αριθμό των συμβόλων που παρήχθησαν.
% k είναι ο αριθμός των bits/σύμβολο, ώστε L=2^k,
% M ο αριθμός των παραγόμενων συμβόλων (μήκος ακολουθίας L-ASK)
% nsamp ο αριθμός των δειγμάτων ανά σύμβολο (oversampling ratio)
% EbNo είναι ο λόγος Eb/No, σε db
L=2^k;
SNR=EbNo-10*log10(nsamp/2/k); % SNR ανά δείγμα σήματος
% Διάνυσμα τυχαίων ακεραίων {±1, ±3, ... ±(L-1)}. Να επαληθευτεί
x=2*floor(L*rand(1,Nsymb))-L+1;
Px=(L^2-1)/3; % θεωρητική ισχύς σήματος
sum(x.^2)/length(x); % μετρούμενη ισχύς σήματος (για επαλήθευση)

%h=ones(1,nsamp); h=h/sqrt(h*h'); % κρουστική απόκριση φίλτρου πομπού (ορθογωνικός παλμός μοναδιαίας ενέργειας)
%h=cos(2*pi*(1:nsamp)/nsamp); h=h/sqrt(h*h');
%y=upsample(x,nsamp); % μετατροπή στο πυκνό πλέγμα
%y=conv(y,h); % το προς εκπομπή σήμα
%y=y(1:Nsymb*nsamp); % περικόπτεται η ουρά που αφήνει η συνέλιξη
%n=wgn(1,length(y),10*log10(Px)-SNR);
%ynoisy=awgn(y,SNR,'measured'); % θορυβώδες σήμα
%ynoisy=y;
%for i=1:nsamp matched(i)=h(end-i+1); end
%yrx=conv(ynoisy,matched);

%z = yrx(nsamp:nsamp:Nsymb*nsamp); % Yποδειγμάτιση -- στο τέλος κάθε περιόδου Τ
y=rectpulse(x,nsamp);
stem(y(1:100));
n=wgn(1,length(y),10*log10(Px)-SNR);
ynoisy=y+n;          % θορυβώδες σήμα
y=reshape(ynoisy,nsamp,length(ynoisy)/nsamp);
h=cos(2*pi*(1:nsamp)/nsamp); h=h/sqrt(h*h');
matched=h;
for i=1:nsamp matched(i)=h(end-i+1); end
yrx=conv(ynoisy,matched);
z = yrx(nsamp:nsamp:M*nsamp); % Yποδειγμάτιση -- στο τέλος
% κάθε περιόδου Τ   

%figure; stem(x(1:20));
%figure; stem(y(1:20*nsamp));
%figure; stem(yrx(1:20*nsamp));
l=[-L+1:2:L-1];
for i=1:length(z)
    [m,j]=min(abs(l-z(i)));
    z(i)=l(j);
end
err=not(x==z);
errors=sum(err);
