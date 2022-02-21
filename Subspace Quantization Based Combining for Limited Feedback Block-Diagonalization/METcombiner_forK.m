function [H_MET,G_MET] = METcombiner_forK(H,L)
[M,N,K] = size(H);
H_MET = zeros(M,N,K);
G_MET = zeros(N,N,K);
for idx = 1:1:K
    [Uk, ~, Vk] = svd(H(:,:,idx));
    H_MET(:,:,idx) = Uk(:,1:L);
    G_MET(:,:,idx) = Vk(:,1:L);
end