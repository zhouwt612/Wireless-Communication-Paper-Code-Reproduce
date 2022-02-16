function C = RVQ_MIMO(M,N,B)
% C: codebook for K users
% D: M dimension
% B: feedback bits
% K: the number of users

len = 2^B;
C = zeros(M,N,len);
for idx1 = 1:1:len
    Codeword = randn(M,N) + 1i*randn(M,N);
    C(:,:,idx1) = Codeword/norm(Codeword);
end