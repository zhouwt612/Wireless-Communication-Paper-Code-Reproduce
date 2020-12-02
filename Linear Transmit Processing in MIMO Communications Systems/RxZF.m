function G = RxZF(P,H,R_n)
% Receive Zero-Forcing Filter
% Decoder

G = ((P'*H'/R_n*H*P)^(-1))*P'*H'/R_n;