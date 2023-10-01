function errors=ask_errors(k,M,nsamp,EbNo)
 % Η συνάρτηση αυτή εξομοιώνει την παραγωγή και αποκωδικοποίηση
 % θορυβώδους σήματος L-ASK και μετρά τον αριθμό των εσφαλμένων
%συμβόλων.
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
 x=(2*floor(L*rand(1,M))-L+1);

 l=[-L+1:2:L-1];
 hist(x, l)
 
 Px=(L^2-1)/3; % θεωρητική ισχύς σήματος
 sum(x.^2)/length(x); % μετρούμενη ισχύς σήματος (για επαλήθευση)
 
 h=ones(1,nsamp); h=h/sqrt(h*h'); % κρουστική απόκριση φίλτρου
 % πομπού (ορθογωνικός παλμός μοναδιαίας ενέργειας)
 y=upsample(x,nsamp); % μετατροπή στο πυκνό πλέγμα
 y=conv(y,h); % το προς εκπομπή σήμα
 y=y(1:M*nsamp); % περικόπτεται η ουρά που αφήνει η συνέλιξη

 %n=wgn(1,length(y),10*log10(Px)-SNR);
 ynoisy= y     %awgn(y,SNR,'measured'); % θορυβώδες σήμα1
 
 
 matched=ones(1,nsamp);
 
 for i=1:nsamp matched(i)=h(end-i+1); end
 yrx=conv(ynoisy,matched);
 z = yrx(nsamp:nsamp:M*nsamp); % Yποδειγμάτιση -- στο τέλος
 % κάθε περιόδου Τ

 hist(z, 200) 

 l=[-L+1:2:L-1];
 for i=1:length(z)
 [m,j]=min(abs(l-z(i)));
 z(i)=l(j);
 
 end
 err=not(x==z);
 errors=sum(err);
 end