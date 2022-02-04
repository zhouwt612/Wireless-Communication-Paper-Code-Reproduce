function precoder = RegZFprecoder(H,P)
% G: regularized zero-forcing precoder
% H: channel
% M: the number of transmit antennas
% P: the transmit power
[K,M] = size(H);
PRE = H'/(H*H' + (M/P)*eye(K));
b = sqrt(P/trace(PRE*PRE'));
precoder = b*PRE;