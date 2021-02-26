function Hls = LSestimator(S,P)
% Least square estimator
% S: received signal
% P: pilot matrix
% Additional: None

Hls = S*P'/(P*P');