function [ber,numBits] = msk_ber_func(EbNo, maxNumErrs, maxNumBits)
% Import Java class for BERTool.
import com.mathworks.toolbox.comm.BERTool;
% Initialize variables related to exit criteria.
totErr = 0; % Number of errors observed
numBits = 0; % Number of bits processed
% �. --- Set up parameters. ---
% --- INSERT YOUR CODE HERE.
k=1; % number of bits per symbol
Nbits=2000; % number of symbols in each run
nsamp=128; % oversampling,i.e. number of samples per T
% Simulate until number of errors exceeds maxNumErrs
% or number of bits processed exceeds maxNumBits.
while((totErr < maxNumErrs) && (numBits < maxNumBits))
% Check if the user clicked the Stop button of BERTool.
if (BERTool.getSimulationStop)
break;
end
% �. --- INSERT YOUR CODE HERE.
errors=msk_errors(Nbits,nsamp,EbNo);
% Assume Gray coding: 1 symbol error ==> 1 bit error
totErr=totErr+errors;
numBits=numBits + k*Nbits;
end % End of loop
% Compute the BER
ber = totErr/numBits;