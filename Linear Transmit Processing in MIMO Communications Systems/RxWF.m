function G = RxWF(P,H,R_s,R_n)
% Receive Wiener Filter
% Decoder

G = ((P'*H'/R_n*H*P + inv(R_s))^(-1))*P'*H'/R_n;