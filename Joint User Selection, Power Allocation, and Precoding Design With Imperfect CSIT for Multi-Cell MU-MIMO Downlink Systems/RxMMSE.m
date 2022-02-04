function A_wmmse = RxMMSE(H,B,Rn)
% H (QxP): Channel realization
% B (PxQ): Transmitter filter (Beamformer)
% Rn (QxQ): Noise covariance matrix
A_wmmse = B'*H'/(H*B*B'*H' + Rn);