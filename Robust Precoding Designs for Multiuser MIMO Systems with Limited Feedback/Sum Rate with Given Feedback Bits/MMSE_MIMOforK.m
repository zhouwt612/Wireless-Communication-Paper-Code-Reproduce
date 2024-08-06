function PMMSE = MMSE_MIMOforK(H,pow)
% MMSE precoding
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
% PMMSE: precoding matrices
[M,N,K] = size(H);
sigma2 = M/pow;
Hdedicated = [];
for idx1 = 1:1:K
    Hdedicated = [Hdedicated H(:,:,idx1)];
end
H = Hdedicated;
P= inv(M*H*H'+sigma2*eye(K*N))*H;
PMMSE = zeros(M,N,K);
for idx = 1:1:K
    precoder = P(:,(idx-1)*N+1:N*idx);
    precoder = Hbasis(precoder);
    PMMSE(:,:,idx) = sqrt(pow/K/trace(precoder*precoder'))*precoder;
end