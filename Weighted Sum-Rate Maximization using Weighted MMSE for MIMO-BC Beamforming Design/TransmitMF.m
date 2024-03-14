function P = TransmitMF(H,pow)
[M,N,K] = size(H);
P = zeros(M,N,K);
for idx = 1:1:K
    P(:,:,idx) = sqrt(pow/K/trace(H(:,:,idx)*H(:,:,idx)'))*H(:,:,idx);
end