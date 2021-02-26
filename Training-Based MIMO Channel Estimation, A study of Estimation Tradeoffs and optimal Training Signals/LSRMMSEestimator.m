function Hlsrmmse = LSRMMSEestimator(S,P,r,t,sigma2,pow)
% LS-RMMSE estimator
% S: received signal
% P: pilot matrix
% r,t: the number of transmit and receive antennas
% sigma2: the power of receiver noise (sigma^2)
% pow: transmit power
% Additional: noise power, transmit power

Hlsrmmse = t*trace(S*S')/(pow*(trace(S*S')+sigma2*r*t))*S*P';