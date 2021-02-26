function [Jmmse,Jmmse_eigdecom] = MMSEmse(r,P,sigma2,Rh)
% Equation 28 & 34
% r: the number of receive antennas
% P: pilot sequences
% sigma2: the power of noise
% Rh: channel correlation matrix
% Jmmse: MMSE estimation error
% Jmmse_eigdecom: MMSE estimation error with eigenvalue decomposition

% Eq. 28
Jmmse = trace(inv(inv(Rh) + (1/(sigma2*r))*(P*P')));
% Eq. 34
[Q, lambda] = eig(Rh);
P_tilde = (1/sqrt((sigma2*r)))*Q'*P;
Jmmse_eigdecom = trace(inv(inv(lambda) + P_tilde*P_tilde'));
