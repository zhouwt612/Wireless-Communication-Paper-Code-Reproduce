function F = quantizedchannel_MIMO(Htilde,C)
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