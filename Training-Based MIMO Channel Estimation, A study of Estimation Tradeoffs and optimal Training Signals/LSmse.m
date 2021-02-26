function Jls = LSmse(snr_linear,t,r)
% LS channel estimation error (Mean square error)
% r,t: the number of transmit and receive antennas
% snr: signal-to-noise ratio (pow/sigma^2) (linear)

Jls = (1/snr_linear)*t^2*r;