function PMMSE = MMSE_MIMOforK(H,pow)
[M,N,K] = size(H);
sigma2 = M/pow;
Hdedicated = [];
for idx1 = 1:1:K
    Hdedicated = [Hdedicated H(:,:,idx1)];
end
P = Hdedicated*inv(Hdedicated'*Hdedicated + sigma2*eye(K*N));
PMMSE = zeros(M,N,K);
for idx = 1:1:K
     Pbasis = P(:,(idx-1)*N+1:N*idx);
     Pbasis = orth(Pbasis);
     PMMSE(:,:,idx) = sqrt(pow/K/N)*Pbasis;
end