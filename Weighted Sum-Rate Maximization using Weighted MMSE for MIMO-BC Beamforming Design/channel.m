function H = channel(M,N,K)

H = zeros(M,N,K);
for idx = 1:1:K
    H(:,:,idx) = (randn(M,N)+1i*randn(M,N))/sqrt(2);
end