function Htilde = Hbasis(H)
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
% H: channel realization for K users (size: (M,N,K))
% Output:
% Htilde: channel subspaces
[M,N,K] = size(H);
Htilde = zeros(M,N,K);
for idx1 = 1:1:K
%     Htilde(:,:,idx1) = orth(H(:,:,idx1));
    [V,~] = eig(H(:,:,idx1)*H(:,:,idx1)');
    Htilde(:,:,idx1) = V(:,end-N+1:end);
end