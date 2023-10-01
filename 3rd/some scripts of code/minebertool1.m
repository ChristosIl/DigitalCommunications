    % Clear all variables, close all open figures, and clear command window
    clear all;
    close all;
    clc
    
    % Define system parameters
    k = 4; % Number of bits per symbol, calculated as 20233 mod 2 + 3
    M = 40000; % Number of symbols
    L = 2^k; % Number of amplitude levels
    nsamp = 16; % Number of samples per symbol
    numBits = k * M; % Total number of bits
    
    % Define Eb/No values and calculate theoretical bit error rate (BER)
    EbNo = 1:20;
    % also changing EbNo from dB to Watt.
    Pe = ((L-1)/L) * erfc(sqrt(3*k/(L^2-1) * (10.^(EbNo/10)))); 
    theoreticalBER = Pe / k;
    
    % Simulate BER for each Eb/No value
    for i = 1:20
        symbolErrors(i) = ask_errors(k, M, nsamp, i);
        simulatedBER(i) = symbolErrors(i) / numBits;
    end
    
    % Plot the simulated and theoretical BER performance
    figure(3);
    hold on;
    set(gca, 'yscale', 'log');
    semilogy(EbNo, simulatedBER, 'r+');
    semilogy(EbNo, theoreticalBER, 'b-');
    title('BER of 16-ASK');
    xlabel('EbNo dB');
    ylabel('BER');
    hold off;
