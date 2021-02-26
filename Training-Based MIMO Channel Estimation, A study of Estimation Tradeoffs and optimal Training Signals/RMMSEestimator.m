function Hrmmse = RMMSEestimator(S,P,trRh,r,t,sigma2)
% Relaxed minimum mean square error estimator
% S: received signal
% P: pilot matrix
% trRh: the trace of channel correlation matrix
% r,t: the number of transmit and receive antennas
% sigma2: the power of receiver noise (sigma^2)
% Additional: the trace of channel correlation matrix, noise power

Hrmmse = S/(P'*P+sigma2*r*t/trRh*eye(size(P'*P)))*P';