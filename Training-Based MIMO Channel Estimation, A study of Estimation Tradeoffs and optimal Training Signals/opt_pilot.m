function P = opt_pilot(N,r,t,sigma2,Rh,mu_o,U)
% Equation 44
% N: the length of pilots
% r: the number of receive antennas
% t: the number of transmit antennas
% pow: power constraint
% sigma2: the power of noise
% Rh: channel correlation matrix
% U: an arbitrary NxN unitary matrix (Identity matrix is fine)
if nargin < 7
    U = eye(N);
end

[Q,lambda] = eig(Rh);
Premtx = mu_o*eye(t) - inv(lambda);
Premtx(Premtx<0) = 0;
Premtx2 = [sqrt(Premtx) zeros(t,N-t)];
P = sqrt(sigma2*r)*Q*Premtx2*U;