function C = RVQcodebook(D,B)
% C: codebook
% D: M dimension
% B: feedback bits

len = 2^B;
realpart = randn(D,len);
imagpart = randn(D,len);

RP = zeros(D,len);
IP = zeros(D,len);
for idx = 1:1:len
    RP(:,idx) = realpart(:,idx)/norm(realpart(:,idx));
    IP(:,idx) = imagpart(:,idx)/norm(imagpart(:,idx));
end
C = (RP+1i*IP)/sqrt(2);