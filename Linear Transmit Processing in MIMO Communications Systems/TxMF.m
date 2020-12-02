function P = TxMF(H,G,R_s,Etr)
% Tranmit Matched Filter
% Precoder

beta_txmf = sqrt(Etr/trace(H'*G'*R_s*G*H));
P = beta_txmf*H'*G';