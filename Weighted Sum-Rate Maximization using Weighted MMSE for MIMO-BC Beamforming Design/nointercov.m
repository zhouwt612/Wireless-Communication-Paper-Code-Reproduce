function Rvv = nointercov(H,B,P,Q,K)
% Rn: Effective noise covariance matrix at user k
Rvv = zeros(Q,Q,K);
B_ite = zeros(P,Q,K);
H_ite = zeros(Q,P,K);

for j = 1:1:K
    B_ite(:,:,j) = B(:,(j-1)*Q+1:j*Q);
    H_ite(:,:,j) = H((j-1)*Q+1:j*Q,:);
end

for k = 1:1:K
    for i = 1:1:K
        if i ~= k
            Rvv(:,:,k) = Rvv(:,:,k) + H_ite(:,:,k)*B_ite(:,:,i)*B_ite(:,:,i)'*H_ite(:,:,k)';
        else
            Rvv(:,:,k) = Rvv(:,:,k) + zeros(Q);
        end
    end 
end
Rvv = Rvv + eye(Q);