function H = channel(M,N,K)
% M: the number of transmit antennas
% N: the number of receive antennas
% K: the number of users

H = (randn(M,N,K) + 1i*randn(M,N,K))/sqrt(2);