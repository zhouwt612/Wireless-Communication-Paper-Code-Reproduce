function G_SQBC = SQBCcombiner_forK(H,Hhat)
% This function is for reproduce of "Antenna combiners for block-diagonalization
% based multi-user MIMO with limited feedback"
[~,N,K] = size(H);
G_SQBC = zeros(N,N,K);
for idx1 = 1:1:K
    [Uk,Sk,Vk] = svd(H(:,:,idx1));
    Lambdak = Sk(1:N,1:N)^(-2);
    Wk = Hhat(:,:,idx1)'*Uk(:,1:N);
    G_SQBC(:,:,idx1) = Vk*Sk(1:N,1:N)^(-1)*Uk(:,1:N)'*Hhat(:,:,idx1)*Wk*Lambdak^(-1/2);
end