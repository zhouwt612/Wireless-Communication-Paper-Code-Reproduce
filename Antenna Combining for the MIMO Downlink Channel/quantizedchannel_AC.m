function F = quantizedchannel_AC(spanH,C)
% spanH: the space spanned by the channel realization (size: M, N, K)
  % spanH = [q1,...,qN]
% M: the number of transmit antennas
% N: the number of receive antennas
% C: codebook (size: M, L=(2^B), K), therefore k is for the k-th user
% K: the number of user equipments
% F: the quantized version of real channel H (size: M,K)

[M,N,K] = size(spanH);
[~,L,~] = size(C);
F = zeros(M,1,K);
for idx1 = 1:1:K
    summax = zeros(1,L);
    for idx2 = 1:1:L
        for idx3 = 1:1:N
            summax(idx2) = summax(idx2)+ abs(C(:,idx2,idx1)'*spanH(:,idx3,idx1))^2;
        end
    end
    [~,index] = max(summax);
    F(:,1,idx1) = C(:,index,idx1);
end