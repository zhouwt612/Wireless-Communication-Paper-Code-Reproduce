function [H_sort,Hhat_sort] = Channel_Sort_Imp(H,Hhat,pow,sigma2)
[M,N,K] = size(Hhat);
sig_inv = 1/sigma2;
Rk = zeros(1,K);
for idx1 = 1:1:K
cvx_begin quiet
    variable Xk(M,M) complex semidefinite
    maximize(log_det(eye(N)+ sig_inv*Hhat(:,:,idx1)'*Xk*Hhat(:,:,idx1)))
    subject to
        %real(trace(Xk)) <= pow;
        norm(Xk,'fro') <= sqrt(pow);
cvx_end
Rk(idx1) = cvx_optval;
end
[~,I] = sort(Rk,'descend');
H_sort = zeros(M,N,K);
Hhat_sort = zeros(M,N,K);
for idx2 = 1:1:K
    H_sort(:,:,idx2) = H(:,:,I(idx2));
    Hhat_sort(:,:,idx2) = Hhat(:,:,I(idx2));
end