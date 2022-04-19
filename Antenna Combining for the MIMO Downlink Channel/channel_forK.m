function H = channel_forK(M,N,K)
% M: the number of transmit antennas
% N: the number of receive antennas
% K: the number of users

H = zeros(M,N,K);
for idx1 = 1:1:K
    H(:,:,idx1) = (randn(M,N) + 1i*randn(M,N))/sqrt(2);
end