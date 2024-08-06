function H = channel(M,N,K)
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
% M: the number of transmit antennas
% N: the number of receive antennas
% K: the number of users
% Output:
% H: channel realizations
H = (randn(M,N,K) + 1i*randn(M,N,K))/sqrt(2);