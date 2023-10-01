 clear all; close all;
 % Το αρχείο "sima.mat" περιέχει το σήμα s και τη συχνότητα
 % δειγματοληψίας Fs. Το φάσμα του σήματος εκτείνεται σχεδόν σε όλη την
 % περιοχή συχνοτήτων μέχρι 4 KHz. Πάνω από 1 KHz, όμως, είναι θόρυβος
 % και πρέπει να φιλτραριστεί.
 load sima;
 %Fs=8192Hz
 %t = 0:1/Fs:1-1/Fs Διάρκεια 1sec
 %s = sin(2*pi*700*t) +sin(2*pi*900*t)+sin(2*pi*1400*t)+sin(2*pi*2500*t);

 figure; pwelch(s,[],[],[],Fs);
 pause
 % Ορίζεται η ιδανική βαθυπερατή συνάρτηση Η, με συχνότ. αποκοπ. Fs/8
 H=[ones(1,Fs/8) zeros(1,Fs-Fs/4) ones(1,Fs/8)];
 % Υπολογίζεται η κρουστική απόκριση με αντίστροφο μετασχ. Fourier
 % Εναλλακτικά, μπορεί να χρησιμοποιηθεί η αναλυτική σχέση Sa(x)

 h=ifft(H,'symmetric');
 h = fftshift(h);
   %παμε απο φασμα συχνοτήτων μονο στα θετικα. Το μετατρέπουμε στο χρονο και 
   %μετα το shiftαρουμε για να παρουμε αμφιπλευρο σήμα χρόνου.

 middle=length(h)/2;
 
 h32=h(middle+1-16:middle+17);
 h64=h(middle+1-32:middle+33);
 h128=h(middle+1-64:middle+65);
 h160=h(middle+1-80:middle+81);
 % figure; stem([0:length(h64)-1],h64); grid;
 % figure; freqz(h64,1); % σχεδιάζουμε την απόκριση συχνότητας της h64
 wvtool(h32,h128,h160);
 pause
 % αποκρίσεις συχνότητας των περικομμένων h
 % Οι πλευρικοί λοβοί είναι υψηλοί!
 % Πολλαπλασιάζουμε την περικομμένη κρουστική απόκριση με κατάλληλο
 % παράθυρο. Χρησιμοποιούμε την h64 και παράθυρα hamming και kaiser
 wh=hamming(length(h160));


 h_hamming=h160.*wh';
 % figure; stem([0:length(h64)-1],h_hamming); grid;
 % figure; freqz(h_hamming,1);
 wvtool(h160,h_hamming);
 % Φιλτράρουμε το σήμα μας με καθένα από τα τρία φίλτρα
 y_rect=conv(s,h160);
 figure; pwelch(y_rect,[],[],[],Fs);
 y_hamm=conv(s,h_hamming);
 figure; pwelch(y_hamm,[],[],[],Fs);
 %
 % Βαθυπερατό Parks-MacClellan
 hpm=firpm(160 , [0 0.10 0.15 0.5]*2, [1 1 0 0]);
 % figure; freqz(hpm,1);
 s_pm=conv(s,hpm);
 figure; pwelch(s_pm,[],[],[],Fs);
 sound(20*s); % ακούμε το αρχικό σήμα, s
 %sound(20*s_lp); % ακούμε το φιλτραρισμένο σήμα, s_lp

 % Βαθυπερατό Parks-MacClellan
 hpm=firpm(160 , [0 0.11 0.12 0.5]*2, [1 1 0 0]);
 % figure; freqz(hpm,1);
 s_pm1=conv(s,hpm);
 figure; pwelch(s_pm1,[],[],[],Fs);
 sound(20*s); % ακούμε το αρχικό σήμα, s
 %sound(20*s_lp); % ακούμε το φιλτραρισμένο σήμα, s_lp