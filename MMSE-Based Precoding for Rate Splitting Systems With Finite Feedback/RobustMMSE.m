function W = RobustMMSE(H,pow,M,B)
% W: MMSE precoder
% H: Channel
% P: transmit power constraint
% M: the number of transmit antennas
% B: feedback bits

delta = 2^(-B/(M-1));
Wb = inv(H'*H+(1+pow*delta)/(pow*(1-delta))*eye(M))*H';
eta = sqrt(pow/(norm(Wb,'fro')^2));
W = eta*Wb;