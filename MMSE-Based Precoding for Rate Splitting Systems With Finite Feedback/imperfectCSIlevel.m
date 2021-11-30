function tau2 = imperfectCSIlevel(B,M)
% B: the number of feedback bits
% M: the number of transmit antennas

tau2 = (2^B + M/(2*(M-1)))^(-1/(M-1));