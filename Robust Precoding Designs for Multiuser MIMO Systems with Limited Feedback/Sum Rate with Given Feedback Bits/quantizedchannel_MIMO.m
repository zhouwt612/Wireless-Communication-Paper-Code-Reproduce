function F = quantizedchannel_MIMO(Htilde,C)
% Calculation for the quantized channel
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
% Htilde: channel subspaces
% C: codebook
% Output:
% F: quantized channels
[M,N,L,K] = size(C);
F = zeros(M,N,K);
for idx1 = 1:1:K
    dsquart = zeros(1,L);
    for idx2 = 1:1:L
        dsquart(1,idx2) = N - trace(Htilde(:,:,idx1)'*C(:,:,idx2,idx1)*C(:,:,idx2,idx1)'*Htilde(:,:,idx1));
    end
    [~,index] = min(dsquart);
    F(:,:,idx1) = C(:,:,index,idx1);
end