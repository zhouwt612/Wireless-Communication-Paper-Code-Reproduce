function H = channel(K,M)
% Channel function
% K: the number of UE (assuming every UE with 1 antenna)
% M: the number of transmit antennas

H=(randn(K,M)+1i*randn(K,M))/sqrt(2);