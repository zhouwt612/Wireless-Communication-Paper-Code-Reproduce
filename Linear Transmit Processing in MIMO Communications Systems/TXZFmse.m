function varepsilonZF = TXZFmse(H,G,R_s,R_n,Etr)
% Compute the MSE of TxZF

Jtx = G*H*H'*G'*(Etr/trace(G*R_n*G'));
varepsilonZF = trace(R_s) - (trace(R_s)^2)/(trace(R_s)+trace(inv(Jtx)*R_s));
