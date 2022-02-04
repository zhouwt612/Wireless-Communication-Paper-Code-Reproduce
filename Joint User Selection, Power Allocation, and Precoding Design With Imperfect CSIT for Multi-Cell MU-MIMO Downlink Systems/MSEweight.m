function W = MSEweight(H,B,Rn,u)
% E: MMSE matrix
% H (QxP): Channel realization
% B (PxQ): Transmitter filter (Beamformer)
% Rn (QxQ): Noise covariance matrix
% u: parameter
if nargin < 4
    u = 1;
end
E = inv(eye(size(H*H')) + B'*H'/Rn*H*B);
W = u*inv(E);