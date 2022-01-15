% 19uec023 - Hitesh Goyal
% Experiment 8 - Observation 2

% This code will perform the convolutinal encoding for two types of rates -
% 1/2 and 1/3. It will also perform Viterbi Decoding for HARD DECISION
% DECODER.

% Further we will find BER values for different ebn values and plot the
% corresponding graphs

clc;
clear;
close all;

% generating random sequence of 0's and 1's
n = 1000000;
R1 = 1/2;
R2 = 1/3;
data = randi([0,1],1,n);

% generating trellis structures for 1/2 and 1/3 convolutional encoders
trellis1 = poly2trellis(3,[7 5]);
trellis2 = poly2trellis(3,[4 5 7]);

% generating codewords
codeword1 = convenc(data,trellis1);
codeword2 = convenc(data,trellis2);

% initialization of arrays for calculation of BER
maxebn = 10;
ber_coded_using2 = zeros(1,maxebn);
ber_coded_using3 = zeros(1,maxebn);

% array to store BER of uncoded system
ber_uncoded = zeros(1,maxebn);
ebndb_values = zeros(1,maxebn);

% Loop to calculate the BPSK modulation and calculate the corresponding BER
% values for different SNR values
for ebndb = 1:maxebn
    ebndb_values(ebndb) = ebndb;
    ebn = 10^(ebndb/10);
    sigma1 = sqrt(1/(2*R1*ebn));
    sigma2 = sqrt(1/(2*R2*ebn));
    noise1 = sigma1*randn(1,n/R1);
    noise2 = sigma1*randn(1,n/R2);
    
    % BPSK Modulation
    tx1= 2*codeword1 -1;
    rx1 = tx1 + noise1;
    
    tx2 = 2*codeword2 -1;
    rx2 = tx2 + noise2;
    
    % Hard Decision Decoder
    recovered_codeword1 = (rx1>0);
    recovered_codeword2 = (rx2>0);
    recovered_data1 = vitdec(recovered_codeword1,trellis1,20,'trunc','hard');
    recovered_data2 = vitdec(recovered_codeword2,trellis2,20,'trunc','hard');
    ner_hard1 = sum(data~=recovered_data1);
    ner_hard2 = sum(data~=recovered_data2);
    ber_coded_using2(ebndb) = ner_hard1/n;
    ber_coded_using3(ebndb) = ner_hard2/n;
    ber_uncoded(ebndb) = qfunc(sqrt(ebn));
end

% Plotting of graph between BER vs. SNR (dB)
semilogy(ebndb_values,ber_coded_using2);
hold on;
semilogy(ebndb_values,ber_coded_using3);
hold on;
semilogy(ebndb_values,ber_uncoded);
legend('BER using 1/2 convolution encoder','BER using 1/3 convolution encoder','BER without convolution coded');
xlabel('SNR (dB) ->');
ylabel('Bit error rate (BER) ->');
title('19uec023 - Hitesh Goyal','BER vs. SNR (dB) for 1/2 and 1/3 convolutional encoder and for uncoded system');
grid on;
hold off;