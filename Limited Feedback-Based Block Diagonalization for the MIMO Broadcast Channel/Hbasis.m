function Htilde = Hbasis(H)
% H: channel realization for K users (size: (M,N,K))

[M,N,K] = size(H);
Htilde = zeros(M,N,K);
for idx1 = 1:1:K
    Htilde(:,:,idx1) = orth(H(:,:,idx1));
end