function Hsls = SLSestimator(S,P,Rh,r,sigma2)
% Scaled least square estimator
% S: received signal
% P: pilot matrix
% Rh: Channel correlation matrix
% sigma2: the power of noise (sigma^2)
% Additional: channel correlation matrix, noise power

Hsls = (trace(Rh)/(sigma2*r*trace(inv(P*P'))+trace(Rh)))*S*P'/(P*P');