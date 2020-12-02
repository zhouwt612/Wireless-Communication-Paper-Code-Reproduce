function A_mmse = RxMMSE(H,B,Rn)
% H (QxP): Channel realization
% B (PxQ): Transmitter filter (Beamformer)
% Rn (QxQ): Noise covariance matrix
A_mmse = B'*H'/(H*B*B'*H' + Rn);
