function XXX = Xcombiner(X)
[M,N,K] = size(X);
XXX = [];
for idx = 1:1:K
    XXX = [XXX, X(:,:,idx)];
end