function varepsilonZF = REZFmse(H,P,R_s,R_n)
% Compute the MSEs of RxZF

Jrx = (R_s^(1/2))'*P'*H'/R_n*H*P*R_s^(1/2);
varepsilonZF = trace(R_s) - (trace(R_s)^2)/(trace(R_s)+trace(inv(Jrx)*R_s));
