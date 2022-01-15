% 19uec023 - Hitesh Goyal
% Experiment 8 - Observation 1

% This code will perform the Convolutional Encoding via both INBUILT
% FUNCTION and via MATLAB coding

clc;
clear;
close all;

% Rates for encoding
R1 = 1/2;
R2 = 1/3;

% Generation of random sequence of 0's and 1's
n = 5;
data = randi([0,1],1,n);

% using inbuilt function poly2trellis for convolutional encoding
trellis1 = poly2trellis(3,[7 5]);
codeword1_inbuilt = convenc(data,trellis1);
disp("Original Data =");
disp(data);
disp("Codeword using inbuilt functions =");
disp(codeword1_inbuilt);

% using matlab coding poly2trellis for convolutional encoding
codeword1_matlab = zeros(1,2*n);
current_data = zeros(1,3);
for i = 1:n
    current_data(2:3) = current_data(1:2);
    current_data(1) = data(i);
    codeword1_matlab(2*i-1) = mod(current_data(1)+current_data(2)+current_data(3) ,2);
    codeword1_matlab(2*i) = mod(current_data(1) + current_data(3) ,2);
end
disp("Codeword using MATLAB = ");
disp(codeword1_matlab);