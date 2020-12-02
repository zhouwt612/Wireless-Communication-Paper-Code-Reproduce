function varepsilonMF = TXMFmse(H,G,R_s,R_n,Etr)
% Compute the MSE of TxMF

Jtx = G*H*H'*G'*(Etr/trace(G*R_n*G'));
varepsilonMF = trace(R_s) - (trace(Jtx*R_s)^2)/trace((Jtx^2+Jtx)*R_s);
