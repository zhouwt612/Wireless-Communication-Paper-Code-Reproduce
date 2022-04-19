function C = RVQforK(M,B,K)
% C: codebook for K users
% D: M dimension
% B: feedback bits
% K: the number of users

len = 2^B;
C = zeros(M,len,K);
for idx1 = 1:1:K
    Codeword = randn(M,len) + 1i*randn(M,len);
    for idx2 = 1:1:len
        C(:,idx2,idx1) = Codeword(:,idx2)/norm(Codeword(:,idx2));
    end
end