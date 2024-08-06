function C_mimo = RVQ_MIMO_QRforK(M,N,B,K)
% Generation for the codebook
%This code refers to the following scientific article:
%
% Wentao Zhou, Di Zhang, Mérouane Debbah, and Inkyu Lee,
% "Robust Precoding Designs for Multiuser MIMO Systems with Limited Feedback,
%" IEEE Transactions on Wireless Communications, To appear.
% 
% This is version 1.0 (last edited: 2024-04-08)
% 
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
% Input:
% M: # of transmit antennas
% N: # of receive antennas
% B: feedback bits
% K: # of users
% Output:
% C_mimo: codebook

% C: codebook for K users
% D: M dimension
% B: feedback bits
% K: the number of users

len = 2^B;
C_mimo = zeros(M,N,len,K);
for idx1 = 1:1:K
    for idx2 = 1:1:len
        C = randn(M,N)+1i*randn(M,N);
        [Q, R] = qr(C);
        C_mimo(:,:,idx2,idx1) = Q(1:N,:)';
    end
end