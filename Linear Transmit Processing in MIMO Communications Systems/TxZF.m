function [P, beta_txzf] = TxZF(H,G,R_s,Etr)
% Transmit Zero-forcing Filter
% Precoder

beta_txzf = sqrt(Etr/trace((G*H*H'*G')^(-1)*R_s));
P = beta_txzf*H'*G'/(G*H*H'*G');
