function V = ACZeroForcing(F)
% F: the quantized version of the channel realization
[M,~,K] = size(F);
Fcom = zeros(M,K);
for idx = 1:1:K
    Fcom(:,idx) = F(:,1,idx);
end
Vb = inv(Fcom');
V = zeros(size(Vb));
for idx2 = 1:1:K
    V(:,idx2) = Vb(:,idx2)/norm(Vb(:,idx2));
end