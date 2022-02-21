function [G_SQBC,Heff] = J_SQBCcombiner_forK(H,L)
% This function is for reproduce of "Subspace Quantization Based Combining
%for Limited Feedback Block-Diagonalization"
[M,N,K] = size(H);
G_SQBC = zeros(N,N,K);
Heff = zeros(M,N,K);
Wk = [eye(L);zeros(N-L,L)];
for idx = 1:1:K
    [Uk,~,~] = svd(H(:,:,idx));
    Bk = Uk(:,1:L);
    Kk = (Wk'*inv(Bk'*H(:,:,idx)*H(:,:,idx)'*Bk)*Wk)^(-1/2);
    Dk = Bk'*H(:,:,idx);
    G_SQBC(:,:,idx) = inv(Dk)*Wk*Kk;
    Heff(:,:,idx) = Bk*Wk*Kk;
end