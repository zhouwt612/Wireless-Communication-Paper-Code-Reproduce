function Heff = effective_channel(H,gamma)
% H: channel realization (size: M,N,K)
% gamma: weighting vector for all users
% Heff: the effective channel

[M,N,K] = size(H);
Heff = zeros(M,1,K);
for idx = 1:1:K
    Heff(:,1,idx) = H(:,:,idx)*gamma(:,:,idx);
end