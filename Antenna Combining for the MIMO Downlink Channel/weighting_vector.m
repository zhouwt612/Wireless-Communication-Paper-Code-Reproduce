function gamma = weighting_vector(H,s_proj)
% H: the channel realization (size: M,N,K)
% s_proj: the projection version of the quantization vector onto spanH
% gamma: the weighting vector (size: M,1,K)

[M,N,K] = size(H);
gamma = zeros(N,1,K);
for idx1 = 1:1:K
    HHHs = inv(H(:,:,idx1)'*H(:,:,idx1))*H(:,:,idx1)'*s_proj(:,1,idx1);
    gamma(:,1,idx1) = HHHs/norm(HHHs);
end