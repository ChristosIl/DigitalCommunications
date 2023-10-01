clear all;
close all;
clc

k = 4; %20233 mod 2 +3
M = 40000; %Nsymb 
L = 2^k;
nsamp = 16;
numBits = k * M;

EbNo = 1:20;
Pe = ((L-1)/L)*erfc(sqrt(3*k/(L^2-1)*(10.^(EbNo/10))));
BER = Pe/k;

for i = 1:20
        errors(i) = ask_errors(k, M, nsamp, i);
        simBER(i) = errors(i)/numBits;
end 


figure(3);
hold on;
set(gca, 'yscale', 'log');
semilogy(EbNo, simBER, 'r+');
semilogy(EbNo, BER, 'b-');
title('BER of 16-ASK');
xlabel('EbNo dB');
ylabel('BER');
hold off;