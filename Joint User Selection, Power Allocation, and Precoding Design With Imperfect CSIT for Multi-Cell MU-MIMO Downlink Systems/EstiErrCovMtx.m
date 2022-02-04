function Phi = EstiErrCovMtx(N,K,phi)
% phi: [phi_1, ... , phi_K]
Phi = zeros(N,N,K);
for idx1 = 1:1:K
    Phi(:,:,idx1) = phi(idx1)*eye(N);
end