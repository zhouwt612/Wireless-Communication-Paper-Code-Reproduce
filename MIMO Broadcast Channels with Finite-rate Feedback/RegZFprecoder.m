function precoder = RegZFprecoder(H,P)
% G: regularized zero-forcing precoder
% H: channel
% M: the number of transmit antennas
% P: the transmit power
[K,M] = size(H);
precoder = H'/(H*H' + (M/P)*eye(K));