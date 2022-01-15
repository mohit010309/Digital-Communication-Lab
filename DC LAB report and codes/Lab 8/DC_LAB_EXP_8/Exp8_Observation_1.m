% 19ucc023
% Mohit Akhouri
% Observation 1 - Generating Convolutional Encoded symbols ( via both INBUILT
% FUNCTION and MATLAB CODING )

% In this code, we will generate convolutionally encoded codewords using
% inbuilt function poly2trellis and also via MATLAB coding with help of
% arrays and loops.

clc;
clear all;
close all;

n = 5; % Number of random numbers to be generated
data = randi([0,1],1,n); % Generates a random sequence of 0's and 1's of length 1xn

% Display of random sequence of binary digits 0's and 1's generated
disp('The Random data generated is : ');
disp(data);

% Using Inbuilt function poly2trellis to generate convolutional codes
codes_trellis = poly2trellis(3,[7 5]); % Generating Trellis structure for convolutional code
codeword_inbuilt = convenc(data,codes_trellis); % Generating convolutional codeword with the help of trellis structure

% Display of convolutional codeword generated via INBUILT FUNCTION
% poly2trellis and convenc
disp('The Convolutional Codeword generated via INBUILT FUNCTION is : ');
disp(codeword_inbuilt);

% Using MATLAB coding to generate convolutional codes
codeword_matlab_coding = zeros(1,2*n); % To store the codeword generated
curr_data = zeros(1,3); % Temporary array to store the binary digits 0 and 1

% Main Loop algorithm for the calculation of convolutional codewords
for i = 1:n
    curr_data(2:3) = curr_data(1:2); % Replacing bits 2 and 3 with bits 1 and 2
    curr_data(1) = data(i); % Starting point of codeword
    
    % Calculation of remaining convolutional codewords with the help of
    % 'mod' function
    codeword_matlab_coding(2*i-1) = mod(curr_data(1) + curr_data(2) + curr_data(3) ,2);
    codeword_matlab_coding(2*i) = mod(curr_data(1) + curr_data(3) ,2);
end

% Display of convolutional codewords generated via MATLAB coding
disp('The Convolutional Codeword generated via MATLAB CODING is : ');
disp(codeword_matlab_coding);