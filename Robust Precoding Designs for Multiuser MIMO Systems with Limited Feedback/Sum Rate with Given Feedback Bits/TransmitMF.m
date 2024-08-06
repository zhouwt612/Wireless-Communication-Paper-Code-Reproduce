function P = TransmitMF(H,pow)
% Maximum ratio transmission precoding
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
% H: channel realizations
% pow: transmit power
% Output:
% P: precoding matrices
[M,N,K] = size(H);
P = zeros(M,N,K);
for idx = 1:1:K
    P(:,:,idx) = sqrt(pow/K/trace(H(:,:,idx)*H(:,:,idx)'))*H(:,:,idx);
end