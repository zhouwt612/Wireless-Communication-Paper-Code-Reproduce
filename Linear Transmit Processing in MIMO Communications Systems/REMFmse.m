function varepsilonMF = REMFmse(H,P,R_s,R_n)
% Compute the MSE of RxMF

Jrx = (R_s^(1/2))'*P'*H'/R_n*H*P*R_s^(1/2);
varepsilonMF = trace(R_s) - (trace(Jrx*R_s)^2)/trace((Jrx^2+Jrx)*R_s);
