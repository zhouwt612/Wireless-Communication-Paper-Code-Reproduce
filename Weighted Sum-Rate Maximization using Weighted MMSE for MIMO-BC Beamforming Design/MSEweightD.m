function WD = MSEweightD(H,B,Rn,u)
% H (QxP): Channel realization
% B (PxQ): Transmitter filter (Beamformer)
% Rn (QxQ): Noise covariance matrix
% u = parameter
if nargin < 4
    u = 1;
end

VAV = B'*H'/Rn*H*B;
[V,D] = eig(VAV);
% B_tilde = B*V;
% ED_inv = eye(size(D))+B_tilde'*H'*Rn*H*B_tilde;
ED_inv = eye(size(D)) + D;
WD = u*ED_inv;