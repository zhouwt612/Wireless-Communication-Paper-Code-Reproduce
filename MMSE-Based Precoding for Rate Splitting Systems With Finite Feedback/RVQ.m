function C = RVQ(D,B)
% C: codebook
% D: M dimension
% B: feedback bits
len = 2^B;
Codeword = randn(D,len) + 1i*randn(D,len);
C = zeros(D,len);
for idx = 1:1:len
    C(:,idx) = Codeword(:,idx)/norm(Codeword(:,idx));
end