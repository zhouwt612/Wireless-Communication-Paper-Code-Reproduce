function P = TxWF(H,G,R_s,R_n,Etr)
% Transmit Wiener Filter
% Precoder

F = H'*G'*G*H + (trace(G*R_n*G')/Etr)*eye(size(H'*G'*G*H));
beta_txwf = sqrt(Etr/trace(F^(-2)*H'*G'*R_s*G*H));
P = beta_txwf*inv(F)*H'*G';
