function varepsilonWF = REWFmse(H,P,R_s,R_n)
% Compute the MSEs of RxWF

Jrx = (R_s^(1/2))'*P'*H'/R_n*H*P*R_s^(1/2);
varepsilonWF = trace(R_s) - trace(inv(Jrx+eye(size(R_s)))*Jrx*R_s);
