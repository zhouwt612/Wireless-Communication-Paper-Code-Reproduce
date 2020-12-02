function varepsilonWF = TXWFmse(H,G,R_s,R_n,Etr)
% Compute the MSE of TxWF

Jtx = G*H*H'*G'*(Etr/trace(G*R_n*G'));
varepsilonWF = trace(R_s) - trace(inv(Jtx+eye(size(R_s)))*Jtx*R_s);
