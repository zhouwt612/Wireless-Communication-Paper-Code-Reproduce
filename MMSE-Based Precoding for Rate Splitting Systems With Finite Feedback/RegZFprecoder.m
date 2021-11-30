function precoder = RegZFprecoder(H,pow)
% G: regularized zero-forcing precoder
% H: channel
% M: the number of transmit antennas
% P: the transmit power
[K,M] = size(H);
W = H'/(H*H' + (M/pow)*eye(K));
precoder = sqrt(pow/trace(W*W'))*W;