function s_proj = quantizedvector_proj(spanH,F)
% spanH: the space spanned by the channel realization (size: M, N, K)
  % spanH = [q1,...,qN]
% F: the quantized version of real channel H (size: M,K)
% s_proj: the projection version of the quantization vector onto spanH

[M,N,K] = size(spanH);
s_proj = zeros(M,1,K);
for idx1 = 1:1:K
    sk = zeros(M,1);
    for idx2 = 1:1:N
        sk = sk + spanH(:,idx2,idx1)*(F(:,1,idx1)'*spanH(:,idx2,idx1));
    end
    s_proj(:,1,idx1) = sk/norm(sk);
end