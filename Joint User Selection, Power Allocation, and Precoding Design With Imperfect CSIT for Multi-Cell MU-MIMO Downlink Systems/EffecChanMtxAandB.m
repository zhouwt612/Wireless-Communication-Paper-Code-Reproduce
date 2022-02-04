function [A,B] = EffecChanMtxAandB(QC,Phi,sigma2,P)
% QC: Quantized Channel (N, K) (same format as real channel)
% Phi: (N, N), the error covariance matrix for all users
% sigma2: the sum of the inter-cell-interference and noise power
% P: power constraint
% A and B: the effective channel matrix (NK, NK, K)
[K,N] = size(QC);
Ainit = zeros(N*K,N*K,K);
for idx1 = 1:1:K
    for idx2 = 1:1:K
        Ainit((N*(idx2-1)+1):(N*idx2),(N*(idx2-1)+1):(N*idx2),idx1) = QC(idx1,:)'*QC(idx1,:) + Phi(:,:,idx1);
    end
end
A = Ainit + (sigma2/P)*eye(N*K,N*K);

MinusTerm = zeros(N*K,N*K,K);
for idx3 = 1:1:K
    MinusTerm((N*(idx3-1)+1):(N*idx3),(N*(idx3-1)+1):(N*idx3),idx3) = QC(idx3,:)'*QC(idx3,:);
end
B = A - MinusTerm;