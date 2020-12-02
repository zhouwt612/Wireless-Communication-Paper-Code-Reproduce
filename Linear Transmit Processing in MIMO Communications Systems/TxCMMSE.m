function P = TxCMMSE(H,G,R_s,Etr,xi)
% Transmit Constrained MMSE transmit filter
% Precoder

P_tilde = inv(H'*G'*G*H + xi*eye(2))*H'*G';
beta = sqrt(Etr/trace(P_tilde*R_s*P_tilde'));
P = beta * P_tilde;
