function C = RVQ_MIMOforK(M,N,B,K)
% C: codebook for K users
% D: M dimension
% B: feedback bits
% K: the number of users

len = 2^B;
C = zeros(M,N,len,K);
for idx1 = 1:1:K
    for idx2 = 1:1:len
        Codeword = randn(M,N) + 1i*randn(M,N);
        C(:,:,idx2,idx1) = Codeword/norm(Codeword);
    end
end