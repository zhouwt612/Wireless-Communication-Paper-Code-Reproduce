function spanH = Space_SpanH(H)
% H: the channel realization (size: M,N,K)
% spanH: the space spanned by the channel realization (size: M,N,K)

[M,N,K] = size(H);
spanH = [];
for idx = 1:1:K
    spanH(:,:,idx) = orth(H(:,:,idx));
end