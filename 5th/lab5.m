clc;
clear all;
close all;

k = 4;
l = k/2;
L = 2^l;
M = L^2;
Nsymb = 30000;
nsamp = 16;

fc = 10;
EbNo = 16;
SNR = EbNo-10*log10(nsamp/k/2);

W = (11.25-8.75)*10^6;
R = 12*10^6;
a = 0.25;

core=[1+i;1-i;-1+i;-1-i];
mapping = core;

if(l>1)
    for j=1:l-1
        mapping=mapping+j*2*core(1);
        mapping=[mapping;conj(mapping)];
        mapping=[mapping;-conj(mapping)];
    end
end

x=floor(2*rand(k*Nsymb,1)); % τυχαία δυαδική ακολουθία
xsym=bi2de(reshape(x,k,length(x)/k).','left-msb')';
y=[];

for n=1:length(xsym)
        y=[y mapping(xsym(n)+1)];
end


%% Ορισμός φίλτρου μορφοποίησης

delay = 8; % Group delay (# περιόδων Τ)
filtorder = delay*nsamp*2;
rolloff = a; % συντελεστής εξάπλωσης φίλτρου
rNyquist = rcosine(1,nsamp,'fir/sqrt',rolloff,delay);


%% Εκπεμπόμενο σήμα
ytx=upsample(y,nsamp);
ytx = conv(ytx,rNyquist);

figure(1); 
pwelch(ytx,[],[],[],nsamp*2);

% quadrature modulation
m=(1:length(ytx));
s=real(ytx.*exp(1j*2*pi*fc*m/(2*nsamp)));

figure(2); 
pwelch(s,[],[],[],2*nsamp);

% ---------------------
% Διάγραμμα οφθαλμού, όταν ορθ. παλμός.
% if (pulse_type==1) eyediagram(ytx(1:2000),nsamp*2);
% Προσθήκη λευκού γκαουσιανού θορύβου
%Ps=10*log10(s*s'/length(s)); % ισχύς σήματος, σε db
%Pn=Ps-SNR; % αντίστοιχη ισχύς θορύβου, σε db
%n=sqrt(10^(Pn/10))*randn(1,length(ytx));
snoisy=awgn(s, SNR, 'measured'); % θορυβώδες ζωνοπερατό σήμα 
clear s;
clear ytx;

%%%%%%%%% ΔΕΚΤΗΣ %%%%%%%%%%%%%


% Αποδιαμόρφωση
yrx = 2*snoisy.*exp(-1j*2*pi*fc*m/(2*nsamp)); %clear s;

figure(3); 
pwelch(yrx,[],[],[],2*nsamp);

yrx = conv(yrx,rNyquist);
%yrx = downsample(yrx,nsamp); % Υποδειγμάτιση στο πλέγμα nT.


figure(4); 
pwelch(real(yrx),[],[],[],2*nsamp);

yrx = downsample(yrx,nsamp); % Υποδειγμάτιση στο πλέγμα nT.

%yrx = downsample(yrx,nsamp); % Υποδειγμάτιση στο πλέγμα nT.
yrx = yrx(2*delay:length(y)-2*delay); % περικοπή άκρων συνέλιξης.

yi=real(yrx); yq=imag(yrx); % συμφασική και εγκάρσια συνιστώσα
xrx=[]; % διάνυσμα δυαδικής ακολουθίας εξόδου -- αρχικά κενό
q=[-L+1:2:L-1];
for n=1:length(yrx) % επιλογή πλησιέστερου σημείου
    [m,j]=min(abs(q-yi(n)));
    yi(n)=q(j);
    [m,j]=min(abs(q-yq(n)));
    yq(n)=q(j);
    m=1;
    while (mapping(m)~=yi(n)+i*yq(n)) 
        m=m+1; 
    end
    xrx=[xrx; de2bi(m-1,k,'left-msb')'];
end

scatterplot(yrx); 
title("Scatter plot, Eb/No=16dB");

ber1=sum(not(y==(yi+1i*yq)))/length(x);
ber2=sum(not(xrx==x))/length(x);