clc 
clear all;
close all;

%variables
k=3;
L=2^k; %filter 8-ASK
num_bits=9999;
Nsymb = num_bits/k; %mhkos ths eksomioumenhs akolouthias
nsamp = 32; %syntelesths uperdigmatishs
delay = 4   %group delay
filterorder = delay*nsamp*2; %filter's order = 256, nsamp =32 and delay =4T
rolloff=0.5; %rolloff factor, suntelesths ptwshs
rNyquist = rcosine(1, nsamp, 'fir/sqrt', rolloff, delay); %

%DHMIOURGIA 9999 bits anamesa se 0 kai 1
x = round(rand(1, num_bits));

%kwdikopoihsh gray
step = 2; %apostash geitwnikwn shmeiwn L-ASK
mapping = [step/2; -step/2];
if(k>1)
    for j=2:k
        mapping = [mapping+2^(j-1)*step/2; -mapping-2^(j-1)*step/2];
    end
end;
xsym = bi2de(reshape(x, k, length(x)/k).','left-msb');
y=[];
for i=1:length(xsym)
        y=[y mapping(xsym(i)+1)];
end

%shma poy ekpempoume
y1 = upsample(y, nsamp);

%figure(5)
%stem(y1(1:10*nsamp));

ytx = conv(y1, rNyquist);


%eyediagram(ytx(1:2000), nsamp*2);
%an theloume thorubo xreiazomaste snr
%clear y1;

%shma pou lambanoume

yrx = conv(ytx, rNyquist);%filtrarisma shmatos

yrx = yrx(2*delay*nsamp+1:end - 2*delay*nsamp); %perikoph logw kathisterhshs



figure(1)
grid; %enwsh twn duo grafikwn
plot(yrx(1:10*nsamp)); hold;
stem([1:nsamp:nsamp*10], y(1:10), 'filled');

figure(2)
pwelch(yrx, [],[],[], nsamp);

figure(3)
stem(yrx(1:30*nsamp));

figure(4)
plot(yrx(1:10*nsamp));




