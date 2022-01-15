% 19ucc023
% Mohit Akhouri
% Observation 2 - Calculation and Plotting of Probability of Error vs.
% Modulation order

% This code will first calculate the Probability of error for different
% modulation order and Plot the graphs between them

clc;
clear all;
close all;

val1 = 5; % SNR = 5 dB stored in val1
val2 = 10; % SNR = 10 dB stored in val2
M=[2 4 8 16 32]; % Initializing the modulation order (M) array

PBE_for_SNR_5 = zeros(1,5); % Initializing array 1 for PBE for SNR = 5 dB
PBE_for_SNR_10 = zeros(1,5); % Initializing array 2 for PBE for SNR = 10 dB

% Calculation of Probability of Error for SNR = 5 dB

SNR_1 = 10.^(val1/10); % Calculating SNR 1

for i=1:5
    PBE_for_SNR_5(i) = 2 * qfunc(sqrt(SNR_1))*sin(pi/M(i));
end

% Calculation of Probability of Error for SNR = 10 dB

SNR_2 = 10.^(val2/10); % Calculating SNR 2

for i=1:5
    PBE_for_SNR_10(i) = 2*qfunc(sqrt(SNR_2))*sin(pi/M(i));
end

% Displaying the values of probability of error for different values of M
display('Probability of error (for SNR = 5dB) values for M=2,4,8,16,32 :');
display(PBE_for_SNR_5);

display('Probability of error (for SNR = 10dB) values for M=2,4,8,16,32 :');
display(PBE_for_SNR_10);

% Plots of Probability of error vs. Modulation order M for SNR = 5 dB
figure;
plot(PBE_for_SNR_5);
xlabel('Modulation order (M) ->');
ylabel('Probability of error ->');
title('19ucc023 - Mohit Akhouri','Probability of error vs. Modulation order (M) for SNR = 5 dB');
grid on;

% Plots of Probability of error vs. Modulation order M for SNR = 10 dB
figure;
plot(PBE_for_SNR_10);
xlabel('Modulation order (M) ->');
ylabel('Probability of error ->');
title('19ucc023 - Mohit Akhouri','Probability of error vs. Modulation order (M) for SNR = 10 dB');
grid on;