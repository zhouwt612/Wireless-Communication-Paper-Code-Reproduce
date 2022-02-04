function F = quantizedchannel(H,C)
% H: channel
% C: codebook
% C = {w1,...,wi}: codeword
% F: quantized channel

% H: K x M  C: M x N, where N = 2^B
[K,M] = size(H);
[~,N] = size(C);
F = zeros(K,M);
hwmulti = zeros(1,N);
for idx1 = 1:1:K
    for idx2 = 1:1:N
        hwmulti(1,idx2) = abs(H(idx1,:)*C(:,idx2));
    end
    [~,index] = max(hwmulti);
    F(idx1,:) = C(:,index)';
end