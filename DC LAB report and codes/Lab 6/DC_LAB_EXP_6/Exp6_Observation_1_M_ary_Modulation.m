% 19ucc023
% Mohit Akhouri
% Observation 1 - Modulation order and Probability of Error for M-ary PSK

clc;
clear all;
close all;

size1 = 10; % intializing the size for BER and SER
size2 = 10000; % intializing the size for signal x[n]

SER_Practical = zeros(1,size1); % Initializing Practical SER array
SER_Theoretical = zeros(1,size1); % Initializing Theoretical SER array

BER_Real_part = zeros(1,size1); % Initializing Real part of BER_Practical
BER_Imag_part = zeros(1,size1); % Initializing Imaginary part of BER_Practical

BER_Practical = zeros(1,size1); % Initializing Practical BER array
BER_Theoretical = zeros(1,size1); % Initializing Theoretical BER array

x_real = zeros(1,size2); % Initializing Real part of x[n]
x_img = zeros(1,size2); % Initializing Imaginary part of x[n]
x = zeros(1,size2); % Intializing input signal x[n]
y = zeros(1,size2); % Initializing Output signal y[n]

% Main loop algorithm for calculation of transmitted signal x[n]
for i=1:size2
    
    rnd1 = rand(); % random value 1 generation
    rnd2 = rand(); % random value 2 generation
    
    % Constructing the real part of signal on the basis of decision
    if(rnd1 > 0.5)
        x_real(i)=1;
    else
        x_real(i)=-1;
    end
    
    % Constructing the Imaginary part of signal on the basis of decision
    if(rnd2 > 0.5)
        x_img(i)=1;
    else
        x_img(i)=-1;
    end
end

x = x_real + (1j * x_img); % Overall transmitted signal x[n]

SNR_dB = 0:9; % defining the range of Signal to Noise Ratio ( Measured in dB )

% Main loop algorithm for calculation of x[n],y[n], noise "n"
% and calculation of theoretical and practical BER and SER
for i=1:length(SNR_dB)
    
    SNR=10^((i-1)/10);
    N = 1/SNR;
    M = sqrt(N/2);
        
    n=zeros(1,size2); % Initializing noise signal n
    
    % loop for calculation of noise signal
    for j=1:size2
        n(j) = M*randn() + ( 1j * M * randn());
    end
    
    % loop for calculation of received AWGN + x[n] signal
    for j=1:size2
        y(j) = x(j) + n(j);
    end
    
    yn = zeros(1,size2);
    y_real = zeros(1,size2);
    y_img = zeros(1,size2);
    
    % Main Loop algorithm for ML-Detection of M-ary and QPSK Modulation
    for j=1:size2
        if(real(y(j)) >= 0)
            y_real(j) = 1;
        else
            y_real(j) = -1;
        end
        
        if(imag(y(j))>=0)
            y_img(j) = 1;
        else
            y_img(j) = -1;
        end
        
        yn = y_real + (1j * y_img);
    end
    
    % Comparing the transmitted and received message signal 
    % and calculating the Practical BER and SER
    for j=1:size2
        if(x(j)~=yn(j))
            SER_Practical(i) = SER_Practical(i) + 1;
        end
        
        if(real(x(j)) ~= real(yn(j)))
            BER_Real_part(i) = BER_Real_part(i) + 1;
        end
        
        if(imag(x(j)) ~= imag(yn(j)))
            BER_Imag_part(i) = BER_Imag_part(i) + 1;
        end
    end
    
    BER_Practical(i) = ((BER_Real_part(i)/size2) + (BER_Imag_part(i)/size2))/2;
    BER_Theoretical(i) = qfunc(sqrt(2/N));
    
    SER_Practical(i) = SER_Practical(i)/size2;
    SER_Theoretical(i) = 2 * qfunc(sqrt(2/N));
end

SER_matrix = zeros(10,5); % Matrix for storing SER values for different modulation orders
M = [2,4,8,16,32]; % array of modulation orders M

% Loop for calculation of SER for different modulation order M
for i=1:length(SNR_dB)
    SNR = 10^((i-1)/10);
    N = 1/SNR;
    
    for m = 1:length(M)
        SER_matrix(i,m) = 2 * qfunc(sqrt(2/N) * sin(pi/M(m)));
    end
    
end

% Displaying the SER matrix
disp('SER vs. Modulation order matrix is given as:');
disp(SER_matrix); 

% Plot of SER for different modulation order
figure;
semilogy(SNR_dB,SER_matrix(:,1),'color','blue');
hold on;
semilogy(SNR_dB,SER_matrix(:,2),'color','black');
semilogy(SNR_dB,SER_matrix(:,3),'color','red');
semilogy(SNR_dB,SER_matrix(:,4),'color','magenta');
semilogy(SNR_dB,SER_matrix(:,5),'color','cyan');
ylabel('SER ->');
xlabel('SNR(dB) ->');
title('19ucc023 - Mohit Akhouri','Plots of SER for different values of Modulation order (M) for M-ary Phase Shift Keying');
legend('SER for M = 2','SER for M = 4','SER for M = 8','SER for M = 16','SER for M = 32');
grid on;
hold off;

% Plots of practical and theoretical SER vs. SNR ( in dB )
figure;
semilogy(SNR_dB,SER_Practical,'Color','blue');
hold on;
semilogy(SNR_dB,SER_Theoretical,'Color','red');
xlabel('SNR(dB) ->');
ylabel("SER ->");
title('19ucc023 - Mohit Akhouri','Plots of Practical and Theoretical SER vs. SNR (dB) for 4-QPSK Modulation');
legend('Practical SER','Theoretical SER');
grid on;
hold off;

% Plots of practical and theoretical BER vs. SNR ( in dB )
figure;
semilogy(SNR_dB,BER_Practical,'Color','blue');
hold on;
semilogy(SNR_dB,BER_Theoretical,'Color','red');
xlabel('SNR(dB) ->');
ylabel("BER ->");
title('19ucc023 - Mohit Akhouri','Plots of Practical and Theoretical BER vs. BER (dB) for 4-QPSK Modulation');
legend('Practical BER','Theoretical SER');
grid on;
hold off;

