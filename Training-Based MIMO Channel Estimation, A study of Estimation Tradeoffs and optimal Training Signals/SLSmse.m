function Jsls = SLSmse(trRh,t,r,sigma2,pow)
% SLS channel estimation error (Mean square error)
% trRh: the trace of the channel correlation matrix
% r,t: the number of transmit and receive antennas
% sigma2: the power of receiver noise (sigma^2)
% pow: the transmit pow

Jsls = sigma2*t^2*r*trRh/(sigma2*t^2*r+pow*trRh);