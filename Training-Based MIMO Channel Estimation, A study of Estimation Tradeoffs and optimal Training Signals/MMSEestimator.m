function Hmmse = MMSEestimator(S,P,Rh,r,sigma2)
% Minimum mean square error estimator
% S: received signal
% P: pilot matrix
% Rh: Channel correlation matrix
% sigma2: the power of receiver noise (sigma^2)
% Additional: channel correlation matrix, receiver noise power

Ao = inv(P'*Rh*P+sigma2*r*eye(size(P'*P)))*P'*Rh;
Hmmse = S*Ao;